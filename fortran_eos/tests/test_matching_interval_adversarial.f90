program test_matching_interval_adversarial
  use compose_eos_types
  use compose_eos_matching
  implicit none
  real(wp)::u(401),nb(401),p(401),f(401),mu(401),q(401),cs(401)
  real(wp)::ux(1),nx(1),px(1),fx(1),mx(1),qx(1),cx(1)
  real(wp),parameter::p2=2.1041881294006204_wp,d0=2.1043168390302784_wp,d1=1.5780789687536598_wp
  integer::j
  logical::ok,valid,point_ok
  character(len=512)::msg,vmsg
  do j=1,401;u(j)=real(j-1,wp)/400._wp;enddo
  call construct_bridge('quintic_P',u,1._wp,3._wp,.01_wp,p2,.99_wp,2._wp-p2/3._wp, &
      1._wp,2._wp,.01_wp,p2/3._wp,1._wp/d0,1.5_wp/d1,nb,p,f,mu,q,cs,ok,msg)
  if(ok)then
    call validate_bridge(nb,p,f,mu,cs,valid,vmsg)
    ux(1)=.8919032228381553_wp
    call construct_bridge('quintic_P',ux,1._wp,3._wp,.01_wp,p2,.99_wp,2._wp-p2/3._wp, &
        1._wp,2._wp,.01_wp,p2/3._wp,1._wp/d0,1.5_wp/d1,nx,px,fx,mx,qx,cx,point_ok,msg)
    print '(a,l1,a,es24.16,a,es24.16)','401_NODE_ACCEPTED=',valid,' CS2_NODE_MAX=',maxval(cs),' CS2_INTERIOR=',cx(1)
    if(valid.and.cx(1)>1._wp)then
      print '(a)','FAIL: acausal interval accepted by matching constructor'
      stop 1
    endif
  else
    print '(a,a)','CERTIFICATE_REJECTION=',trim(msg)
    if(index(msg,'interval')==0)stop 2
  endif
  print '(a)','matching interval counterexample: PASS'
end program
