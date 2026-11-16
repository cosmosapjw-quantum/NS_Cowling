program native_spg_probe
 use compose_eos_types,only:wp
 use compose_eos_gdfm_native
 use compose_eos_spg_native
 implicit none
 type(gdfm_native_model_t)::model
 type(spg_native_result_t)::state,drip_state,cross_state
 type(spg_native_result_t),allocatable::attempts(:)
 type(gdfm_native_beta_t)::uniform
 character(160)::cmd,id,ss,arg
 real(wp)::v(7),e,g(4),eps,mun,mup,p,nb,nc,nu,cross,nd,zlast(4)
 integer::i,ios,j,npoints
 call get_command_argument(1,cmd);call get_command_argument(2,id)
 call load_gdfm_native(trim(id),model,ss)
 if(ss/='ok')then
  print *,trim(ss);stop 2
 endif
 v=0
 do i=3,command_argument_count()
  call get_command_argument(i,arg);read(arg,*,iostat=ios)v(i-2)
  if(ios/=0)stop 3
 enddo
 select case(trim(cmd))
 case('bulk')
  call gdfm_native_baryons(model,v(1),v(2),eps,mun,mup,p,ss)
  print '(4es25.16e3)',eps,mun,mup,p
 case('cell')
  call spg_native_cell(model,v(1),v(2),v(3),v(4:7),e,g,state,ss)
  print '(*(es25.16e3,1x))',e,g,state%yn,state%yp,state%ye,state%ymu,state%yh
 case('solve')
  call spg_native_solve(model,v(1),v(2),v(3),state,ss)
  print '(*(es25.16e3,1x))',state%epsilon,state%pressure,state%mu_b,state%mu_e,state%max_projected_gradient, &
    state%yn,state%yp,state%ye,state%ymu,state%yh,state%z,state%min_scaled_curvature
 case('attempts')
  call spg_native_solve(model,v(1),v(2),v(3),state,ss,attempts=attempts)
  do i=1,size(attempts)
   print '(*(es25.16e3,1x))',attempts(i)%energy_per_baryon,attempts(i)%z, &
     attempts(i)%gradient,attempts(i)%max_projected_gradient,attempts(i)%gas_proton_margin
   print '(a)',trim(attempts(i)%status)
  enddo
 case('beta')
  call gdfm_native_beta(model,v(1),uniform,ss)
  print '(*(es25.16e3,1x))',uniform%epsilon,uniform%pressure,uniform%mu_n,uniform%mu_e,uniform%yp,uniform%ymu
 case('drip')
  call spg_native_drip(model,v(1),v(2),nb,state,ss)
  print '(*(es25.16e3,1x))',nb,state%pressure,state%mu_b,state%max_projected_gradient,state%z
 case('crossing')
  call spg_native_crossing(model,v(1),v(2),v(3:4),cross,state,uniform,ss)
  print '(*(es25.16e3,1x))',cross,state%epsilon-uniform%epsilon,state%pressure,state%max_projected_gradient
 case('transition')
  call spg_native_transition(model,v(1),v(2),cross,state,uniform,ss)
  print '(*(es25.16e3,1x))',cross,state%epsilon-uniform%epsilon,state%pressure,state%max_projected_gradient
 case('grid')
  npoints=nint(v(3));if(npoints<3)stop 3
  call spg_native_drip(model,v(1),v(2),nd,drip_state,ss)
  if(ss/='ok')then
   print '(a)',trim(ss);stop 2
  endif
  call spg_native_transition(model,v(1),v(2),cross,cross_state,uniform,ss)
  if(ss/='ok')then
   print '(a)',trim(ss);stop 2
  endif
  zlast=drip_state%z
  do j=1,npoints
   nb=exp(log(nd)+real(j-1,wp)/real(npoints-1,wp)*log(cross/nd))
   if(j==1)then
    state=drip_state
   elseif(j==npoints)then
    state=cross_state
   else
    call spg_native_solve(model,v(1),v(2),nb,state,ss,start=zlast)
   endif
   if(ss/='ok')then
    print *,nb
    print '(a)',trim(ss);stop 2
   endif
   zlast=state%z
   print '(*(es25.16e3,1x))',state%nb,state%epsilon,state%pressure,state%mu_b,state%mu_e, &
    state%yn,state%yp,state%ye,state%ymu,state%yh,state%aav,state%zav, &
    state%max_projected_gradient,state%min_scaled_curvature,state%z
  enddo
 case('endpoint')
  call spg_native_endpoint(model,v(1),v(2),v(3:4),nc,nu,cross,ss)
  print '(*(es25.16e3,1x))',nc,nu,cross
 case default
  ss='unknown command'
 end select
 print '(a)',trim(ss)
 if(ss/='ok')stop 2
end program
