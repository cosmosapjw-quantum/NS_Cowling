program test_firstlaw_interval
  use compose_eos_types
  use compose_eos_matching
  implicit none
  real(wp)::u(401),nb(401),p(401),f(401),mu(401),q(401),cs(401),ux(1),nx(1),px(1),fx(1),mx(1),qx(1),cx(1)
  real(wp)::kn(4),vv(4),ss(4),pp,dp,anti,x,dx,fc,grid_max,dense_max,umax,f1,f2,mu1,mu2,c1,c2
  integer::i
  logical::ok,point_ok,valid
  character(len=512)::msg,vmsg
  dx=log(2._wp)
  call make_pressure_spline(5._wp,.1_wp,10._wp,.1_wp,.1_wp,kn,vv,ss)
  grid_max=-huge(1._wp);dense_max=grid_max;umax=0._wp
  do i=0,10000
    x=real(i,wp)/10000._wp
    call spline_eval(kn,vv,ss,x,pp,dp,anti)
    fc=(dp/dx-pp)/(.1_wp*exp(x*dx))-dx*pressure_work_integral(5._wp,x,.1_wp,10._wp,.1_wp,.1_wp,log(.1_wp),dx)
    if(mod(i,25)==0)grid_max=max(grid_max,fc)
    if(fc>dense_max)then;dense_max=fc;umax=x;endif
  enddo
  print '(a,3es24.16)','MARGIN_MAXIMA_GRID_DENSE_LOCATION ',grid_max,dense_max,umax
  if(dense_max<=grid_max+1.e-8_wp.or.grid_max<=0._wp)stop 2
  f1=.5_wp*(grid_max+dense_max)
  f2=f1+dx*pressure_work_integral(5._wp,1._wp,.1_wp,10._wp,.1_wp,.1_wp,log(.1_wp),dx)
  mu1=f1+1._wp;mu2=f2+50._wp
  c1=.1_wp/(.1_wp*mu1*dx);c2=.1_wp/(.2_wp*mu2*dx)
  do i=1,401;u(i)=real(i-1,wp)/400._wp;enddo
  call construct_bridge('first_law_constrained_hermite_pressure_bridge',u,.1_wp,.2_wp,.1_wp,10._wp, &
    f1,f2,mu1,mu2,1._wp,50._wp,c1,c2,nb,p,f,mu,q,cs,ok,msg)
  if(.not.ok)then
    print '(a,a)','INTERVAL_REJECTED ',trim(msg)
    if(index(msg,'interval')==0)stop 3
  else
    call validate_bridge(nb,p,f,mu,cs,valid,vmsg)
    ux(1)=umax
    call construct_bridge('first_law_constrained_hermite_pressure_bridge',ux,.1_wp,.2_wp,.1_wp,10._wp, &
      f1,f2,mu1,mu2,1._wp,50._wp,c1,c2,nx,px,fx,mx,qx,cx,point_ok,msg)
    print '(a,l1,2es24.16)','401_ACCEPTED_MAX_CS2_INTERIOR_CS2 ',valid,maxval(cs),cx(1)
    if(valid.and.point_ok.and.cx(1)>1._wp)then
      print '(a)','FAIL: 401-sample gate accepted a superluminal first-law bridge'
      stop 1
    endif
    stop 4
  endif
  f1=f1+100._wp;f2=f2+100._wp
  mu1=f1+1._wp;mu2=f2+50._wp
  c1=.1_wp/(.1_wp*mu1*dx);c2=.1_wp/(.2_wp*mu2*dx)
  call construct_bridge('first_law_constrained_hermite_pressure_bridge',u,.1_wp,.2_wp,.1_wp,10._wp, &
    f1,f2,mu1,mu2,1._wp,50._wp,c1,c2,nb,p,f,mu,q,cs,ok,msg)
  if(.not.ok)then
    print '(a,a)','FAIL_VALID_CONTROL ',trim(msg)
    stop 5
  endif
  call validate_bridge(nb,p,f,mu,cs,valid,vmsg)
  if(.not.valid.or.maxval(cs)>=1._wp)stop 6
  print '(a,es24.16)','VALID_FIRSTLAW_CS2_MAX ',maxval(cs)
  print '(a)','first-law interval counterexample and valid control: PASS'
end program
