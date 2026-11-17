! Native field re-evaluation of retained chord profiles, with no optimization.
program tf_finish_profile_readback
 use compose_eos_types,only:wp
 use compose_eos_spatial_tf
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 type(tf_grid_t)::g
 type(tf_state_t)::s
 character(512)::directory,path,arg
 character(160)::model,seed,geometry,status
 integer::n,count,i,u,d,j
 real(wp)::nb,r,raw_mu,raw_p,mu_kkt,denom,threshold
 real(wp),allocatable::nn(:),np(:)
 call get_command_argument(1,model);call get_command_argument(2,arg);read(arg,*)n
 call get_command_argument(3,directory);call get_command_argument(4,arg);read(arg,*)count
 allocate(nn(n),np(n))
 do i=1,count
  write(path,'(a,a,i3.3,a)')trim(directory),'/state_',i,'.txt'
  open(newunit=u,file=trim(path),status='old',action='read')
  read(u,*)d,seed,nb,r,geometry;read(u,*)nn;read(u,*)np;close(u)
  call make_tf_grid(n,d,g,status)
  if(status/='ok')stop 1
  s=tf_state_t()
  call evaluate_screened_tf_profile(trim(model),nb,g,nn,np,r,s,status)
  print '(a,1x,i0,1x,a)','EVALUATION',i,trim(status)
  if(status/='ok')stop 2
  ! Compare the same active-density multiplier as solve_spatial_tf, while
  ! retaining the evaluator's distinct neutron-volume mean separately.
  raw_mu=s%mu_B;raw_p=s%pressure;mu_kkt=0;denom=0
  threshold=max(1e-10_wp,nb*1e-7_wp)
  do j=1,n
   if(nn(j)>threshold)then
    denom=denom+g%weights(j);mu_kkt=mu_kkt+g%weights(j)*s%mu_n(j)
   endif
   if(np(j)>threshold)then
    denom=denom+g%weights(j);mu_kkt=mu_kkt+g%weights(j)*(s%mu_p(j)+s%mu_Q_lepton)
   endif
  enddo
  if(denom<=0)stop 6
  s%mu_B=mu_kkt/denom;s%pressure=nb*s%mu_B-s%epsilon
  print '(a,1x,i0,11(1x,es25.16e3))','PROFILE_READBACK',i,nb,s%epsilon,s%pressure,s%mu_B, &
   s%scalar_pde_residual,s%vector_pde_residual,s%electron_pde_residual, &
   s%baryon_constraint_error,s%charge_constraint_error,raw_p,raw_mu
  if(.not.all(ieee_is_finite([s%epsilon,s%pressure,s%mu_B])))stop 3
  if(s%scalar_pde_residual>1e-8_wp.or.s%vector_pde_residual>1e-8_wp.or.s%electron_pde_residual>2e-10_wp)stop 4
  if(s%baryon_constraint_error>max(1e-12_wp,1e-8_wp*nb).or.s%charge_constraint_error>max(1e-12_wp,1e-8_wp*nb))stop 5
 enddo
 print '(a)','PASS_NATIVE_RETAINED_PROFILE_FIELD_READBACK'
end program
