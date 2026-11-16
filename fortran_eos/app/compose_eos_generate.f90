! Native physical generation entry point; Python transports inputs and hashes only.
program compose_eos_generate
 use compose_eos_types,only:wp,eos_data_t,allocate_eos
 use compose_eos_io,only:write_compose,read_compose
 use compose_eos_core_native
 use compose_eos_core_domain
 use compose_eos_gdfm_native
 use compose_eos_spg_native
 use compose_eos_tf_models,only:tf_model_t,load_tf_model
 use compose_eos_spatial_tf
 use compose_eos_lower_crust_native
 use compose_eos_crust_attachment
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 character(1024)::request_path,out,mass_catalogue=''
 character(32)::model='',family='',transition_policy='fixed_density_source'
 integer::n_points=161,n_density=0,n_mesh=0,meshes(64)=0,outer_points=65,lower_points=65
 real(wp)::n_min=0,n_max=0,sigma0=0,surface_b=0,densities(4096)=0
 logical::include_lower_crust=.false.
 namelist/generation/model,family,n_points,n_min,n_max,sigma0,surface_b,transition_policy, &
  n_density,densities,n_mesh,meshes,include_lower_crust,mass_catalogue,outer_points,lower_points
 type(eos_data_t)::eos,outer_eos,lower_eos
 type(core_state_t),allocatable::core(:)
 type(core_state_t)::audit,plus1,plus2,minus1,minus2
 type(spg_native_result_t),allocatable::spg(:)
 type(gdfm_native_model_t)::gm
 type(gdfm_native_beta_t)::uniform
 type(spg_native_result_t)::drip,endpoint
 type(tf_state_t),allocatable::attempts(:),all_attempts(:),chosen(:)
 type(tf_state_t)::selected
 type(tf_model_t)::tf_model
 type(mass_catalogue_t)::catalogue
 type(lower_crust_state_t),allocatable::outer_states(:),lower_states(:)
 type(lower_source_seam_t)::seam
 real(wp),allocatable::grid(:)
 real(wp)::lo,hi,mn,mp,cross,n,h(2),err(2),relative(2)=[2e-4_wp,1e-4_wp]
 integer::u,ios,i,j,k,a,idx,unresolved,lower_count,au,pu,su,phase
 integer,parameter::dims(6)=[3,2,1,2,3,3]
 character(8),parameter::seeds(6)=[character(8)::'droplet','droplet','droplet','hole','hole','uniform']
 character(256)::status,family_status
 character(1024)::directory
 if(command_argument_count()/=2)then
  write(*,'(a)')'usage: compose_eos_generate REQUEST_NML OUTPUT_DIR';stop 2
 endif
 call get_command_argument(1,request_path);call get_command_argument(2,out)
 open(newunit=u,file=trim(request_path),status='old',iostat=ios)
 if(ios/=0)then;write(*,'(a)')'cannot read request';stop 2;endif
 read(u,nml=generation,iostat=ios);close(u)
 if(ios/=0)then;write(*,'(a)')'invalid generation namelist';stop 2;endif
 call mkdir(trim(out))
 if(include_lower_crust.and.family=='core')call fail('core cannot include lower crust')
 select case(trim(family))
 case('core')
  call core_model_bounds(trim(model),lo,hi,mn,mp,status);call require_ok(status)
  if(n_points<3.or.n_points>100000.or.n_min<lo.or.n_max>hi.or.n_max<=n_min)call fail('invalid core support request')
  allocate(grid(n_points))
  do i=1,n_points
   grid(i)=n_min*(n_max/n_min)**(real(i-1,wp)/real(n_points-1,wp))
  enddo
  grid(1)=n_min;grid(n_points)=n_max
  call generate_core_native(trim(model),grid,eos,core,status);call require_ok(status)
  open(newunit=u,file=trim(out)//'/states.csv',status='replace')
  write(u,'(a)')'n_b,epsilon,pressure,mu_n,mu_p,mu_e,Yn,Yp,Ye,Ymu,chi,gibbs_residual,charge_residual,euler_residual'
  do i=1,size(core)
   write(u,'(*(es25.16e3,:,","))')core(i)%n_b,core(i)%epsilon,core(i)%pressure,core(i)%mu_n,core(i)%mu_p,core(i)%mu_e, &
    core(i)%yn,core(i)%yp,core(i)%ye,core(i)%ymu,core(i)%chi,core(i)%gibbs_residual,core(i)%charge_residual,core(i)%euler_residual
  enddo
  close(u)
  open(newunit=u,file=trim(out)//'/first_law.csv',status='replace')
  write(u,'(a)')'n_b,h1,h2,error1,error2'
  do i=1,7
   n=exp(log(n_min)+(.1_wp+.8_wp*real(i-1,wp)/6)*log(n_max/n_min))
   call solve_core_state(trim(model),n,audit,status);call require_ok(status)
   do j=1,2
    h(j)=min(n*relative(j),(n-n_min)/4,(n_max-n)/4)
    if(h(j)<32*epsilon(n)*n)call fail('density interval too narrow for derivative audit')
    call solve_core_state(trim(model),n+h(j),plus1,status);call require_ok(status)
    call solve_core_state(trim(model),n+2*h(j),plus2,status);call require_ok(status)
    call solve_core_state(trim(model),n-h(j),minus1,status);call require_ok(status)
    call solve_core_state(trim(model),n-2*h(j),minus2,status);call require_ok(status)
    err(j)=(-plus2%epsilon+8*plus1%epsilon-8*minus1%epsilon+minus2%epsilon)/(12*h(j))-audit%mu_n
   enddo
   write(u,'(*(es25.16e3,:,","))')n,h,err
  enddo
  close(u)
  call write_core_domain_diagnostics()
 case('spg')
  if(n_points<3.or.n_points>100000.or.sigma0<=0.or.surface_b<=-16)call fail('invalid SPG request')
  if(transition_policy/='fixed_density_source')call fail('unsupported SPG transition policy')
  call load_gdfm_native(trim(model),gm,status);call require_ok(status)
  call spg_native_drip(gm,sigma0,surface_b,lo,drip,status);call require_ok(status)
  call spg_native_transition(gm,sigma0,surface_b,cross,endpoint,uniform,status);call require_ok(status)
  hi=cross
  allocate(spg(n_points));call allocate_eos(eos,n_points)
  eos%mn=gdfm_mass;eos%mp=gdfm_mass
  eos%composition_file_present=.true.;eos%composition_available=.true.;eos%beta_equilibrated=.true.
  open(newunit=u,file=trim(out)//'/states.csv',status='replace')
  write(u,'(a)')'n_b,epsilon,pressure,mu_b,mu_e,Yn,Yp,Ye,Ymu,Yh,Aav,Zav,KKT,curvature,gas_proton_margin,muon_margin,radius'
  do i=1,n_points
   if(i==1)then
    spg(i)=drip
   elseif(i==n_points)then
    spg(i)=endpoint
   else
    call spg_native_solve(gm,sigma0,surface_b,lo*(hi/lo)**(real(i-1,wp)/real(n_points-1,wp)),spg(i),status,start=spg(i-1)%z)
    call require_ok(status)
   endif
   if(.not.spg(i)%stationary)call fail('SPG requested row not stationary')
   call fill_row(eos,i,spg(i)%nb,spg(i)%epsilon,spg(i)%pressure,spg(i)%mu_b,spg(i)%mu_e,1)
   eos%yn(i)=spg(i)%yn;eos%yp(i)=spg(i)%yp;eos%ye(i)=spg(i)%ye;eos%ymu(i)=spg(i)%ymu
   eos%yh(i)=spg(i)%yh;eos%abar(i)=spg(i)%aav;eos%zbar(i)=spg(i)%zav
   call complete_masks(eos,i)
   write(u,'(*(es25.16e3,:,","))')spg(i)%nb,spg(i)%epsilon,spg(i)%pressure,spg(i)%mu_b,spg(i)%mu_e, &
    spg(i)%yn,spg(i)%yp,spg(i)%ye,spg(i)%ymu,spg(i)%yh,spg(i)%aav,spg(i)%zav, &
    spg(i)%max_projected_gradient,spg(i)%min_scaled_curvature,spg(i)%gas_proton_margin,spg(i)%muon_margin,spg(i)%radius
   flush(u)
  enddo
  close(u)
  open(newunit=u,file=trim(out)//'/transition.json',status='replace')
  write(u,'(a)')'{"policy":"fixed_density_source","source_surface_calibration":"candidate_unresolved_original",'
  write(u,'(a,es25.16e3,a)')'"n_b":',cross,','
  write(u,'(a,es25.16e3,a)')'"delta_epsilon":',endpoint%epsilon-uniform%epsilon,','
  write(u,'(a,es25.16e3,a)')'"delta_pressure":',endpoint%pressure-uniform%pressure,','
  write(u,'(a,es25.16e3,a)')'"delta_mu_b":',endpoint%mu_b-uniform%mu_n,'}'
  close(u)
 case('tf')
  call load_tf_model(trim(model),tf_model,status);call require_ok(status)
  if(n_density<3.or.n_density>4096.or.n_mesh<1.or.n_mesh>64)call fail('invalid TF request dimensions')
  if(any(meshes(:n_mesh)<8).or.any(densities(:n_density)<=0))call fail('invalid TF mesh or density')
  do i=2,n_density
   if(densities(i)<=densities(i-1))call fail('TF densities must increase')
  enddo
  do j=2,n_mesh
   if(meshes(j)<=meshes(j-1))call fail('TF meshes must increase')
  enddo
  allocate(chosen(n_density))
  open(newunit=au,file=trim(out)//'/attempts.csv',status='replace')
  write(au,'(a)')'node,mesh,dimension,seed,attempt,status,geometry,epsilon,pressure,energy_per_baryon,KKT,radius_stationarity,scalar_residual,vector_residual,electron_residual,baryon_error,charge_error,radius,evaluation_status,optimizer_mode,iterations,evaluations'
  open(newunit=pu,file=trim(out)//'/profiles.csv',status='replace')
  write(pu,'(a)')'node,mesh,attempt,cell,nn,np,ne,sigma,omega,rho,V,mu_n,mu_p,mu_e'
  open(newunit=su,file=trim(out)//'/selected.csv',status='replace')
  write(su,'(a)')'node,n_b,mesh,attempt,geometry,epsilon,pressure,mu_b,mu_e,KKT,radius_stationarity'
  do j=1,n_mesh
   do i=1,n_density
    allocate(all_attempts(0))
    do k=1,6
     call solve_tf_family_native(trim(model),densities(i),meshes(j),dims(k),trim(seeds(k)),selected,family_status,attempts)
     all_attempts=[all_attempts,attempts]
    enddo
    do a=1,size(all_attempts)
     call write_attempt(au,pu,i,meshes(j),a,all_attempts(a))
    enddo
    flush(au);flush(pu)
    call select_tf_candidates_native(all_attempts,dims,seeds,densities(i),idx,status,unresolved,lower_count)
    if(idx<1)call fail('TF '//trim(status)//' at density '//real_string(densities(i)))
    chosen(i)=all_attempts(idx)
    write(su,'(i0,",",es25.16e3,",",i0,",",i0,",",a,6(",",es25.16e3))')i,densities(i),meshes(j),idx,trim(chosen(i)%geometry), &
     chosen(i)%epsilon,chosen(i)%pressure,chosen(i)%mu_B,chosen(i)%mu_Q_lepton,chosen(i)%max_chemical_KKT_residual_MeV,chosen(i)%radius_stationarity_MeV
    flush(su)
    write(*,'(a,a,a,es13.5,a,i0,a,a)')trim(model),' selected ','n=',densities(i),' N=',meshes(j),' ',trim(chosen(i)%geometry)
    deallocate(all_attempts)
   enddo
   call allocate_eos(eos,n_density);eos%mn=tf_model%mn;eos%mp=tf_model%mp
   eos%composition_file_present=.true.;eos%beta_equilibrated=.true.
   do i=1,n_density
    select case(trim(chosen(i)%geometry))
    case('droplet');phase=1
    case('rod');phase=2
    case('slab');phase=3
    case('tube');phase=4
    case('bubble');phase=5
    case('uniform');phase=6
    case default;phase=9
    end select
    call fill_row(eos,i,densities(i),chosen(i)%epsilon,chosen(i)%pressure,chosen(i)%mu_B,chosen(i)%mu_Q_lepton,phase)
    eos%ye(i)=chosen(i)%Y_e;eos%ymu(i)=chosen(i)%Y_mu
    eos%ye_present(i)=.true.;eos%ymu_present(i)=.true.;eos%row_status(i)=1
    if(phase==6)then
     eos%yn(i)=chosen(i)%Y_n;eos%yp(i)=chosen(i)%Y_p;eos%abar(i)=1;eos%zbar(i)=0;eos%yh(i)=0
     call complete_masks(eos,i)
    endif
   enddo
   directory=trim(out)//'/mesh_'//int_string(meshes(j));call save_checked(trim(directory),eos)
  enddo
  close(au);close(pu);close(su)
 case default
  call fail('unknown generation family')
 end select
 call save_checked(trim(out),eos)
 if(include_lower_crust)then
  call load_mass_catalogue(trim(mass_catalogue),catalogue,status);call require_ok(status)
  call prepare_lower_source_native(catalogue,outer_points,lower_points,outer_eos,lower_eos,outer_states,lower_states,seam,status)
  call require_ok(status)
  call save_checked(trim(out)//'/outer',outer_eos)
  call save_checked(trim(out)//'/lower_source_fit',lower_eos)
  open(newunit=u,file=trim(out)//'/crust_segments.json',status='replace')
  write(u,'(a)')'{"construction":"documented_boundary_separate_native_segments","nominal_lower_inner_boundary_fm3":0.002,'
  write(u,'(a)')'"cross_segment_interpolation":"UNAPPROVED","source_equivalence":"representation_and_calibration_dependent",'
  write(u,'(a)')'"outer_lower_seam":{'
  write(u,'(a,es25.16e3,a)')'"delta_n":',seam%delta_nb,','
  write(u,'(a,es25.16e3,a)')'"delta_P":',seam%delta_pressure,','
  write(u,'(a,es25.16e3,a)')'"delta_mu":',seam%delta_mu_b,','
  write(u,'(a,es25.16e3,a)')'"delta_epsilon":',seam%delta_epsilon,'},'
  write(u,'(a,es25.16e3,a)')'"lower_last_density":',lower_eos%nb(lower_eos%n),','
  write(u,'(a,es25.16e3,a)')'"inner_first_density":',eos%nb(1),'}'
  close(u)
  call write_attachment_diagnostics()
 endif
 open(newunit=u,file=trim(out)//'/SUPPORT.json',status='replace')
 write(u,'(a)')'{"scope":"computed_native_samples","global_minimum_proven":false,"continuous_mesh_limit_proven":false,'
 write(u,'(a)')'"cross_phase_interpolation":"UNAPPROVED","cross_segment_interpolation":"UNAPPROVED",'
 write(u,'(a)')'"physical_acceptance":"independent_of_CompOSE_numeric_agreement"}'
 close(u)
 open(newunit=u,file=trim(out)//'/RESULT.json',status='replace')
 write(u,'(a)')'{"backend":"fortran","status":"REQUESTED_ROWS_NATIVE_VERIFIED","native_readback":true,'
 write(u,'(a)')'"model":"'//trim(model)//'","family":"'//trim(family)//'",'
 write(u,'(a,i0,a)')'"rows":',eos%n,',"held_nodes":0,'
 write(u,'(a)',advance='no')'"meshes":['
 do j=1,n_mesh
  if(j>1)write(u,'(a)',advance='no')','
  write(u,'(i0)',advance='no')meshes(j)
 enddo
 write(u,'(a)')'],"lower_crust_segments":'//bool_string(include_lower_crust)//'}'
 close(u)
 write(*,'(a,i0)')'REQUESTED_ROWS_NATIVE_VERIFIED rows=',eos%n
contains
 subroutine write_core_domain_diagnostics()
 integer::v,l
 real(wp)::cs,cm,fl,ca,frozen,hyy,beta,hstep
 character(256)::fd_status,analytic_status
 character(:),allocatable::fd_values,analytic_value
 open(newunit=v,file=trim(out)//'/sound_speed.csv',status='replace')
 write(v,'(a)')'n_b,cs2_over_c2_fd,cs2_over_c2_mu,firstlaw_relative,fd_status,cs2_over_c2_analytic,analytic_status'
 do l=1,size(core)
  hstep=2e-5_wp*core(l)%n_b
  fd_status='stencil outside declared model support'
  if(core(l)%n_b-2*hstep>=lo.and.core(l)%n_b+2*hstep<=hi) &
   call core_sound_speed_fd(trim(model),core(l)%n_b,2e-5_wp,cs,cm,fl,fd_status)
  fd_values=',,'
  if(fd_status=='ok')fd_values=real_string(cs)//','//real_string(cm)//','//real_string(fl)
  call core_sound_speed_analytic(trim(model),core(l)%n_b,ca,frozen,hyy,beta,analytic_status)
  if(model=='APR_APR'.and.index(analytic_status,'mixed-phase')==1) &
   call core_sound_speed_gibbs(core(l)%n_b,ca,analytic_status)
  analytic_value=''
  if(analytic_status=='ok')analytic_value=real_string(ca)
  write(v,'(a)')real_string(core(l)%n_b)//','//fd_values//','//csvquote(fd_status)//','// &
   analytic_value//','//csvquote(analytic_status)
 enddo
 close(v)
 end subroutine
 subroutine write_attachment_diagnostics()
 integer::v
 type(crust_attachment_t)::d
 real(wp)::left(4),right(4)
 character(256)::s
 open(newunit=v,file=trim(out)//'/attachment_diagnostics.json',status='replace')
 write(v,'(a)')'{"interpretation":"source-prescription endpoints; no equilibrium or interpolation certification",'
 left=eos_endpoint(outer_eos,outer_eos%n);right=eos_endpoint(lower_eos,1)
 call source_attachment_diagnostics(left,right,BSK_DRIP,d,s)
 write(v,'(a)')'"outer_to_lower":'
 call write_attachment_record(v,left,right,d,s)
 write(v,'(a)')', "lower_to_inner":'
 left=eos_endpoint(lower_eos,lower_eos%n);right=eos_endpoint(eos,1)
 call source_attachment_diagnostics(left,right,BSK_LOWER_MAX,d,s)
 call write_attachment_record(v,left,right,d,s)
 write(v,'(a)')'}'
 close(v)
 end subroutine
 function eos_endpoint(e,k)result(x)
 type(eos_data_t),intent(in)::e
 integer,intent(in)::k
 real(wp)::x(4)
 x=[e%nb(k),e%nb(k)*e%mn*(1+e%q(7,k)),e%nb(k)*e%q(1,k),e%mn*(1+e%q(3,k))]
 end function
 subroutine write_attachment_record(v,left,right,d,s)
 integer,intent(in)::v
 real(wp),intent(in)::left(4),right(4)
 type(crust_attachment_t),intent(in)::d
 character(*),intent(in)::s
 ! Raw endpoint differences remain reportable even for a malformed attachment.
 write(v,'(a)')'{"status":'//csvquote(s)//','
 write(v,'(a,es25.16e3,a)')'"delta_n":',right(1)-left(1),','
 write(v,'(a,es25.16e3,a)')'"delta_epsilon":',right(2)-left(2),','
 write(v,'(a,es25.16e3,a)')'"delta_P":',right(3)-left(3),','
 write(v,'(a,es25.16e3,a)')'"delta_mu":',right(4)-left(4),','
 write(v,'(a)')'"source_prescription_valid":'//bool_string(d%source_prescription_valid)//','
 write(v,'(a)')'"density_overlap":'//bool_string(right(1)<left(1))//','
 write(v,'(a)')'"density_gap":'//bool_string(right(1)>left(1))//','
 write(v,'(a)')'"equilibrium_established":false,"cross_interval_supported":false}'
 end subroutine
 subroutine require_ok(s)
 character(*),intent(in)::s
 if(trim(s)/='ok')call fail(trim(s))
 end subroutine
 subroutine fail(message)
 character(*),intent(in)::message
 integer::v
 write(*,'(a)')'FAILED: '//trim(message)
 open(newunit=v,file=trim(out)//'/FAILURE.txt',status='replace')
 write(v,'(a)')trim(message);close(v)
 stop 1
 end subroutine
 function int_string(i)result(s)
 integer,intent(in)::i
 character(:),allocatable::s
 character(32)::b
 write(b,'(i0)')i;s=trim(b)
 end function
 function real_string(x)result(s)
 real(wp),intent(in)::x
 character(:),allocatable::s
 character(32)::b
 write(b,'(es24.16e3)')x;s=trim(adjustl(b))
 end function
 function bool_string(b)result(s)
 logical,intent(in)::b
 character(:),allocatable::s
 s='false';if(b)s='true'
 end function
 function shellquote(s)result(q)
 character(*),intent(in)::s
 character(:),allocatable::q
 integer::l
 q="'"
 do l=1,len_trim(s)
  if(s(l:l)=="'")then
   q=q//"'"//achar(34)//"'"//achar(34)//"'"
  else
   q=q//s(l:l)
  endif
 enddo
 q=q//"'"
 end function
 subroutine mkdir(path)
 character(*),intent(in)::path
 integer::code
 call execute_command_line('mkdir -p -- '//shellquote(path),exitstat=code)
 if(code/=0)then;write(*,'(a)')'cannot create output';stop 2;endif
 end subroutine
 subroutine fill_row(e,k,n,energy,p,mub,mue,phase)
 type(eos_data_t),intent(inout)::e
 integer,intent(in)::k,phase
 real(wp),intent(in)::n,energy,p,mub,mue
 if(.not.all(ieee_is_finite([n,energy,p,mub,mue])).or.n<=0)call fail('nonfinite physical row')
 e%nb(k)=n;e%q(1,k)=p/n;e%q(3,k)=mub/e%mn-1;e%q(4,k)=-mue/e%mn
 e%q(6,k)=energy/(n*e%mn)-1;e%q(7,k)=e%q(6,k);e%iphase(k)=phase
 end subroutine
 subroutine complete_masks(e,k)
 type(eos_data_t),intent(inout)::e
 integer,intent(in)::k
 e%yn_present(k)=.true.;e%yp_present(k)=.true.;e%ye_present(k)=.true.;e%ymu_present(k)=.true.
 e%abar_present(k)=.true.;e%zbar_present(k)=.true.;e%yh_present(k)=.true.
 e%c3_applicable(k)=.true.;e%c4_applicable(k)=.true.;e%row_status(k)=2
 end subroutine
 subroutine save_checked(path,e)
 character(*),intent(in)::path
 type(eos_data_t),intent(in)::e
 type(eos_data_t)::back
 character(256)::s
 call mkdir(path)
 call write_compose(path,e,s);call require_ok(s)
 call read_compose(path,back,s);call require_ok(s)
 if(back%n/=e%n)call fail('native row count readback mismatch')
 if(maxval(abs(back%nb-e%nb))>1e-13_wp*maxval(abs(e%nb)))call fail('native density readback mismatch')
 if(maxval(abs(back%q-e%q))>1e-12_wp*max(1._wp,maxval(abs(e%q))))call fail('native thermo readback mismatch')
 if(any(back%iphase/=e%iphase))call fail('native phase readback mismatch')
 if(any(back%yn_present.neqv.e%yn_present).or.any(back%yp_present.neqv.e%yp_present).or. &
 any(back%ye_present.neqv.e%ye_present).or.any(back%ymu_present.neqv.e%ymu_present).or. &
 any(back%abar_present.neqv.e%abar_present).or.any(back%zbar_present.neqv.e%zbar_present).or.any(back%yh_present.neqv.e%yh_present)) &
 call fail('native composition mask readback mismatch')
 if(maxval(abs(back%yn-e%yn))+maxval(abs(back%yp-e%yp))+maxval(abs(back%ye-e%ye))+maxval(abs(back%ymu-e%ymu))+ &
 maxval(abs(back%yh-e%yh))+maxval(abs(back%abar-e%abar))+maxval(abs(back%zbar-e%zbar))>1e-9_wp)call fail('native composition readback mismatch')
 end subroutine
 function csvquote(s)result(q)
 character(*),intent(in)::s
 character(:),allocatable::q
 integer::l
 q='"'
 do l=1,len_trim(s)
  if(s(l:l)=='"')q=q//'"'
  q=q//s(l:l)
 enddo
 q=q//'"'
 end function
 subroutine write_attempt(v,pv,node,mesh,id,s)
 integer,intent(in)::v,pv,node,mesh,id
 type(tf_state_t),intent(in)::s
 integer::l
 write(v,'(i0,",",i0,",",i0,",",a,",",i0,",",a,",",a,11(",",es25.16e3),",",a,3(",",i0))')node,mesh,s%dimension,trim(s%seed),id, &
 trim(s%numerical_status),trim(s%geometry),s%epsilon,s%pressure,s%energy_per_baryon,s%max_chemical_KKT_residual_MeV, &
 s%radius_stationarity_MeV,s%scalar_pde_residual,s%vector_pde_residual,s%electron_pde_residual, &
 s%baryon_constraint_error,s%charge_constraint_error,s%radius_fm,csvquote(s%evaluation_status),s%optimizer_mode,s%iterations,s%evaluations
 if(allocated(s%nn).and.allocated(s%mu_n))then
  do l=1,size(s%nn)
   write(pv,'(i0,",",i0,",",i0,",",i0,10(",",es25.16e3))')node,mesh,id,l,s%nn(l),s%np(l),s%ne(l), &
    s%phi(l),s%omega(l),s%rho(l),s%coulomb(l),s%mu_n(l),s%mu_p(l),s%mu_e(l)
  enddo
 endif
 end subroutine
end program
