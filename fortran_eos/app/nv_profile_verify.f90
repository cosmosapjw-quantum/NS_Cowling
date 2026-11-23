! Read named-column CSV projections on stdin and independently check them against
! the physical CompOSE slice on disk. The front end handles CSV dialect and JSON.
program nv_profile_verify
 use compose_eos_types, only: wp, eos_data_t
 use compose_eos_io, only: read_compose, read_metadata
 use, intrinsic :: iso_fortran_env, only: input_unit, error_unit
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 type(eos_data_t) :: eos, refined_eos
 character(4096) :: primary, refined, path
 character(12), parameter :: models(3)=[character(12)::'printed','d03x100','uniform_beta']
 character(12) :: model
 integer :: k, mesh, fine_mesh, nodes, fine_nodes
 logical :: primary_declared, refined_declared
 real(wp),allocatable :: table(:,:), fine(:,:), anchors(:,:), ye(:)
 character(:),allocatable :: mesh_prefix
 character(40) :: tag
 if(command_argument_count()/=2) error stop 'expected primary and refined profile directories'
 model='mesh'
 call get_command_argument(1,primary); call get_command_argument(2,refined)
 call mesh_info(primary,128,mesh,nodes,primary_declared)
 call mesh_info(refined,256,fine_mesh,fine_nodes,refined_declared)
 call require(fine_mesh>mesh,'refined radial mesh must be finer')
 write(tag,'(a,i0,a,i0,a)') 'mesh_',mesh,'_vs_',fine_mesh,'_'
 mesh_prefix=trim(tag)
 model='root'
 call integer_value('fine_radial_mesh',fine_mesh)
 call integer_value('published_control_radial_mesh',mesh)
 write(*,'(a,1x,l1,1x,l1)') 'mesh_metadata_present',primary_declared,refined_declared
 do k=1,size(models)
    model=models(k)
    path=trim(primary)//'/'//trim(model)
    call read_slice(path,k==3,eos)
    if(k==3)then
       ! n, epsilon, P, Ye, neutron/proton/electron mu, beta residual,
       ! local and sampled charge curvature, FD pressure error and sound speed.
       call input_table(12,eos%n,table)
       ye=table(4,:)
       call uniform_checks()
    else
       ! n, epsilon, P, N, Z, E/A-mn, max density, particle error, charge
       ! affinity, FD affinity error, sector and finite-dilation pressure error.
       call input_table(12,eos%n,table)
       call input_table(3,eos%n,fine)   ! n, P, E/A-mn
       call read_slice(trim(refined)//'/'//trim(model),.false.,refined_eos)
       call require(refined_eos%n==eos%n,'refined CompOSE row count')
       call input_table(3,11,anchors)    ! P, E/A-mn, published E/A-mn
       if(nodes>0)call require(eos%n==nodes,'primary density node count')
       if(fine_nodes>0)call require(eos%n==fine_nodes,'refined density node count')
       call require(all(table(4,:)==950).and.all(table(5,:)==32),'cell composition changed')
       ye=table(5,:)/(table(4,:)+table(5,:))
       call profile_checks()
    endif
    call common_checks()
 enddo
contains
 subroutine require(ok,msg)
   logical,intent(in)::ok
   character(*),intent(in)::msg
   if(ok)return
   write(error_unit,'(a)')trim(model)//': '//msg
   error stop 1
 end subroutine
 subroutine mesh_info(dir,legacy,radial,node_count,declared)
   character(*),intent(in)::dir
   integer,intent(in)::legacy
   integer,intent(out)::radial,node_count
   logical,intent(out)::declared
   integer::u,ios,radial_cells,density_nodes
   character(:),allocatable::name
   namelist /profile_run/ radial_cells,density_nodes
   name=trim(dir)//'/profile_run.nml'
   inquire(file=name,exist=declared)
   radial=legacy;node_count=0
   if(.not.declared)return
   radial_cells=-1;density_nodes=-1
   open(newunit=u,file=name,status='old',action='read',iostat=ios)
   if(ios/=0)error stop 'cannot open profile_run.nml'
   read(u,nml=profile_run,iostat=ios);close(u)
   if(ios/=0.or.radial_cells<32.or.density_nodes<4)error stop 'invalid profile_run.nml'
   radial=radial_cells;node_count=density_nodes
 end subroutine
 subroutine read_slice(dir,beta,x)
   character(*),intent(in)::dir
   logical,intent(in)::beta
   type(eos_data_t),intent(out)::x
   type(eos_data_t)::meta
   character(512)::status
   call read_compose(trim(dir),x,status);call require(status=='ok',trim(status))
   call read_metadata(trim(dir)//'/metadata.nml',meta,status);call require(status=='ok',trim(status))
   call require(x%n>=2.and.x%temperature_mev==0,'cold slice or row count')
   call require((x%beta_equilibrated.eqv.beta).and.(meta%beta_equilibrated.eqv.beta),'beta convention')
   call require(meta%mn==x%mn.and.meta%mp==x%mp.and.meta%charge_fraction_grid==x%charge_fraction_grid,'slice metadata')
   call require(all(ieee_is_finite(x%q)).and.all(ieee_is_finite(x%nb)),'nonfinite CompOSE state')
 end subroutine
 subroutine input_table(columns,rows,values)
   integer,intent(in)::columns,rows
   real(wp),allocatable,intent(out)::values(:,:)
   integer::count_rows,i,ios
   read(input_unit,*,iostat=ios) count_rows
   call require(ios==0,'missing CSV row count')
   call require(count_rows==rows,'unexpected CSV row count')
   allocate(values(columns,rows))
   do i=1,rows
      read(input_unit,*,iostat=ios)values(:,i)
      call require(ios==0,'invalid numeric CSV field')
      call require(all(ieee_is_finite(values(:,i))),'nonfinite numeric CSV field')
   enddo
 end subroutine
 subroutine common_checks()
   real(wp)::n,e,p,last_e,last_p,slope,euler_error,minimum,maximum
   integer::i
   minimum=huge(1._wp);maximum=0;euler_error=0
   last_e=0;last_p=0
   call require(all(eos%ye_present),'missing electrons in composition')
   call require(all(abs(eos%ye-ye)<1e-12_wp),'CSV electron fraction differs from CompOSE')
   do i=1,eos%n
      n=eos%nb(i)
      call require(n>0.and.abs(n-table(1,i))<=1e-14_wp*n,'density mismatch')
      e=(1+eos%q(7,i))*eos%mn*n
      p=eos%q(1,i)*n
      call require(abs(e-table(2,i))<=2e-13_wp*abs(e),'energy mismatch')
      call require(abs(p-table(3,i))<=2e-13_wp*abs(p),'pressure mismatch')
      if(i>1)then
         call require(n>eos%nb(i-1).and.e>last_e,'nonmonotone density or energy')
         slope=(p-last_p)/(e-last_e)
         call require(slope>0.and.slope<1,'pressure secant outside (0,1)')
         minimum=min(minimum,slope);maximum=max(maximum,slope)
      endif
      euler_error=max(euler_error,abs((e+p)/n-eos%mn*(1+eos%q(3,i)+ye(i)*eos%q(5,i))))
      last_e=e;last_p=p
   enddo
   call require(euler_error<1e-8_wp,'generalized Euler identity')
   call integer_value('rows',eos%n)
   call value('n_min_fm-3',eos%nb(1)); call value('n_max_fm-3',eos%nb(eos%n))
   call value('cs2_over_c2_secant_range',minimum,maximum)
   call value('generalized_euler_max_abs_MeV',euler_error)
 end subroutine
 subroutine profile_checks()
   integer::split
   real(wp)::dp
   call require(all(table(1,:)==fine(1,:)),'mesh density grid mismatch')
   call require(all(abs(fine(1,:)-refined_eos%nb)<=1e-14_wp*refined_eos%nb),'refined CompOSE density')
   call require(all(abs(fine(2,:)-refined_eos%q(1,:)*refined_eos%nb)<=2e-13_wp*abs(fine(2,:))),'refined CompOSE pressure')
   call require(all(abs(fine(3,:)-refined_eos%mn*refined_eos%q(7,:))<1e-10_wp),'refined CompOSE energy per baryon')
   call require(all(abs(table(6,:)-eos%mn*eos%q(7,:))<1e-10_wp),'primary CompOSE energy per baryon')
   call require(all(table(7,:)>=0.and.table(7,:)<=.15_wp),'DME density domain')
   call require(all(table(8,:)>=0.and.table(8,:)<1.e-8_wp),'particle-number residual')
   call require(all(abs(table(9,:)-eos%mn*eos%q(5,:))<1e-9_wp),'charge affinity differs from CompOSE')
   call require(all(table(11,:)==eos%iphase),'sector label mismatch')
   split=count(eos%iphase==1)
   call require(split>0.and.split<eos%n,'missing shape sector')
   call require(all(eos%iphase(:split)==1).and.all(eos%iphase(split+1:)==2),'shape partition')
   dp=maxval(abs(table(3,:)-fine(2,:)))
   call require(dp<5e-5_wp,'radial refinement pressure')
   call require(all(table(12,:)>=0.and.table(12,:)<2e-7_wp),'dilation pressure derivative')
   call require(all(table(10,:)>=0.and.table(10,:)<2e-5_wp),'composition derivative')
   call integer_value('clustered_rows',split)
   call value('last_clustered_fm-3',table(1,split));call value('first_homogeneous_fm-3',table(1,split+1))
   call value('max_density_any_species_fm-3',maxval(table(7,:)))
   call value('mu_L_range_MeV',minval(table(9,:)),maxval(table(9,:)))
   call value('max_particle_number_residual',maxval(table(8,:)))
   call value('max_dilation_pressure_relative_error',maxval(table(12,:)))
   call value('max_FD_charge_derivative_error_MeV',maxval(table(10,:)))
   call value(mesh_prefix//'max_abs_pressure_MeV_fm-3',dp)
   call value(mesh_prefix//'max_abs_E_per_baryon_MeV',maxval(abs(table(6,:)-fine(3,:))))
   call value('anchor_E_per_baryon_minus_table3_range_MeV',minval(anchors(2,:)-anchors(3,:)),maxval(anchors(2,:)-anchors(3,:)))
   call integer_value('negative_pressure_anchor_rows',count(anchors(1,:)<0))
 end subroutine
 subroutine uniform_checks()
   call require(all(eos%yn_present).and.all(eos%yp_present),'missing uniform nucleon fractions')
   call require(all(abs(eos%yn+eos%yp-1)<1e-12_wp).and.all(abs(eos%yp-ye)<1e-12_wp),'uniform neutrality or nucleon closure')
   call require(all(abs(table(5,:)-eos%mn*(1+eos%q(3,:)))<1e-9_wp),'uniform neutron chemical potential')
   call require(all(abs(table(7,:)+eos%mn*eos%q(4,:))<1e-9_wp),'uniform electron chemical potential')
   call require(all(abs(table(6,:)+table(7,:)-table(5,:)-table(8,:))<1e-9_wp),'CSV beta affinity inconsistent')
   call require(all(abs(table(8,:)-eos%mn*eos%q(5,:))<1e-9_wp).and.maxval(abs(table(8,:)))<1e-8_wp,'uniform beta equilibrium')
   call require(all(table(9,:)>0).and.all(table(12,:)>0),'uniform charge curvature')
   call require(all(table(10,:)>=0.and.table(10,:)<1e-4_wp),'uniform first-law FD pressure')
   call require(all(table(11,:)>0.and.table(11,:)<1),'uniform local sound speed')
   call value('Ye_range',minval(ye),maxval(ye))
   call value('max_abs_beta_residual_MeV',maxval(abs(table(8,:))))
   call value('max_first_law_P_rel_error',maxval(table(10,:)))
   call value('analytic_charge_curvature_min_sampled_MeV_fm3',minval(table(12,:)))
   call value('cs2_over_c2_local_range',minval(table(11,:)),maxval(table(11,:)))
 end subroutine
 subroutine integer_value(key,i)
   character(*),intent(in)::key
   integer,intent(in)::i
   write(*,'(a,1x,a,1x,i0)')trim(model),key,i
 end subroutine
 subroutine value(key,x,y)
   character(*),intent(in)::key
   real(wp),intent(in)::x
   real(wp),intent(in),optional::y
   if(present(y))then
      write(*,'(2(a,1x),2(es25.16e3,1x))')trim(model),key,x,y
   else
      write(*,'(2(a,1x),es25.16e3)')trim(model),key,x
   endif
 end subroutine
end program
