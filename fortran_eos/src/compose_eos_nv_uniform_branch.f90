! Beta-equilibrated homogeneous comparison for the NV1973 Table 1 functional.
module compose_eos_nv_uniform_branch
 use compose_eos_types, only: wp,eos_data_t,allocate_eos,validate_eos
 use compose_eos_nv_original, only: nv_original_source_t,nv_uniform_state_t, &
      nv1973_uniform_beta_state,nv1973_uniform_composition
 use compose_eos_io, only: write_compose,write_metadata,read_compose
 implicit none
 private
 public :: write_nv_uniform_branch
contains
 subroutine write_nv_uniform_branch(src,out,status)
  type(nv_original_source_t),intent(in)::src
  character(*),intent(in)::out
  character(*),intent(out)::status
  integer,parameter::nodes=257
  real(wp),parameter::delta=1e-5_wp
  integer::i,j,u
  real(wp)::n,n1,n2,yp,cs2,fd,euler,energy_old,pressure_old,ytest,ymax,emincheck,mntol
  type(nv_uniform_state_t)::s,probe,lo,hi
  type(eos_data_t)::eos,back
  logical::ok
  call allocate_eos(eos,nodes)
  eos%mn=src%fit%mn
  eos%mp=src%fit%mp
  eos%beta_equilibrated=.true.
  eos%composition_file_present=.true.
  eos%composition_available=.true.
  eos%nT=1
  eos%nYq=1
  eos%t_first=1
  eos%t_last=1
  eos%yq_first=1
  eos%yq_last=1
  eos%selected_iT=1
  eos%selected_iYq=1
  eos%yn_present=.true.
  eos%yp_present=.true.
  eos%ye_present=.true.
  eos%ymu_present=.true.
  ! The uniform phase contains no nuclei. Do not fabricate a zero-A tuple.
  eos%row_status=1
  open(newunit=u,file=trim(out)//'/uniform_beta.csv',status='replace')
  write(u,'(a)')'nb_fm3,epsilon_MeV_fm3,pressure_MeV_fm3,Ye,mu_n_MeV,mu_p_MeV,mu_e_MeV,'// &
     'beta_residual_MeV,beta_curvature_MeV_fm3,first_law_P_rel_error,euler_MeV,cs2_over_c2_FD,'// &
     'min_sampled_energy_excess_MeV,minimum_sampled_charge_curvature_MeV_fm3'
  do i=1,nodes
   n=src%fit%nlo*(src%fit%nhi/src%fit%nlo)**(real(i-1,wp)/(nodes-1))
   if(i==1)n=src%fit%nlo
   if(i==nodes)n=src%fit%nhi
   call nv1973_uniform_beta_state(src,n,s,status)
   if(status/='ok')exit
   yp=s%np/n
   ! Compare energies and analytic charge curvature across the complete
   ! proton-fraction interval on a mixed fine-logarithmic/linear probe grid.
   emincheck=huge(1._wp)
   mntol=huge(1._wp)
   ymax=1-1e-10_wp
   do j=0,256
    if(j<=128)then
     ytest=10._wp**(-11+9._wp*j/128)
    else
     ytest=.01_wp+(ymax-.01_wp)*(j-128)/128
    endif
    call nv1973_uniform_composition(src,n,n*ytest,probe,status)
    if(status/='ok')exit
    emincheck=min(emincheck,probe%e_offset-s%e_offset)
    mntol=min(mntol,probe%beta_curvature)
    if(probe%beta_residual*(ytest-yp)<-1e-9_wp)then
     status='unexpected additional uniform composition stationary point'
     exit
    endif
   enddo
   if(status/='ok')exit
   if(emincheck< -1e-10_wp)then
    status='lower energy composition on independent sample'
    exit
   endif
   n1=max(src%fit%nlo,n*(1-delta))
   n2=min(src%fit%nhi,n*(1+delta))
   call nv1973_uniform_beta_state(src,n1,lo,status)
   if(status/='ok')exit
   call nv1973_uniform_beta_state(src,n2,hi,status)
   if(status/='ok')exit
   fd=n*n*(hi%e_offset-lo%e_offset)/(n2-n1)
   cs2=(hi%pressure-lo%pressure)/(hi%energy-lo%energy)
   euler=(s%energy+s%pressure)/n-s%mu_n
   if(abs(fd-s%pressure)>1e-4_wp*s%pressure.or.cs2<=0.or.cs2>=1.or.abs(euler)>1e-9_wp)then
    status='uniform beta thermodynamic derivative check failed'
    exit
   endif
   if(i>1)then
    if(s%energy<=energy_old.or.s%pressure<=pressure_old)then
     status='uniform beta decreases with density'
     exit
    endif
   endif
   energy_old=s%energy
   pressure_old=s%pressure
   eos%nb(i)=n
   eos%q(1,i)=s%pressure/n
   eos%q(3,i)=s%mu_n/eos%mn-1
   eos%q(4,i)=-s%mu_e/eos%mn
   eos%q(5,i)=s%beta_residual/eos%mn
   eos%q(6,i)=s%e_offset/eos%mn
   eos%q(7,i)=s%e_offset/eos%mn
   eos%yn(i)=1-yp
   eos%yp(i)=yp
   eos%ye(i)=yp
   write(u,'(*(es25.16e3,:,","))')n,s%energy,s%pressure,yp,s%mu_n,s%mu_p,s%mu_e,s%beta_residual, &
     s%beta_curvature,abs(fd-s%pressure)/s%pressure,euler,cs2,emincheck,mntol
  enddo
  close(u)
  if(status/='ok')return
  call validate_eos(eos,ok,status)
  if(.not.ok)return
  call write_compose(trim(out)//'/uniform_beta',eos,status)
  if(status/='ok')return
  call write_metadata(trim(out)//'/uniform_beta/metadata.nml',eos,status)
  if(status/='ok')return
  call read_compose(trim(out)//'/uniform_beta',back,status)
  if(status/='ok')return
  if(.not.back%beta_equilibrated.or.any(abs(back%q-eos%q)>1e-13_wp))status='uniform beta round-trip mismatch'
 end subroutine
end module
