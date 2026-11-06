! Regular-centre series for the stated cold barotropic TOV/Cowling equations.
! These coefficients change initialization accuracy, not the differential law.
module ns_center_series
 use ns_precision, only: wp
 use ns_physical_constants, only: g_si,c_si
 use ns_status
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public :: tov_center_second, cowling_center_second
contains
 subroutine tov_center_second(ec,pc,chi,r,m,p,nu_offset,status)
  real(wp),intent(in)::ec,pc,chi,r
  real(wp),intent(out)::m,p,nu_offset
  type(status_t),intent(out)::status
  real(wp)::g1,g3,m3,m5,p2,p4,e2,pi
  m=0._wp;p=0._wp;nu_offset=0._wp
  call fail_status(status,NS_INVALID_INPUT,'physics','invalid central-series input')
  if(.not.all(ieee_is_finite([ec,pc,chi,r])))return
  if(min(ec,pc,chi,r)<=0._wp.or.chi>1._wp)return
  pi=acos(-1._wp)
  g1=4._wp*pi*g_si/(c_si**4)*(ec/3._wp+pc)
  m3=4._wp*pi*ec/(3._wp*c_si**2)
  p2=-(ec+pc)*g1/2._wp;e2=p2/chi
  m5=4._wp*pi*e2/(5._wp*c_si**2)
  g3=g_si/c_si**2*(m5+4._wp*pi*p2/c_si**2+ &
      (m3+4._wp*pi*pc/c_si**2)*2._wp*g_si*m3/c_si**2)
  p4=-((ec+pc)*g3+(1._wp+1._wp/chi)*p2*g1)/4._wp
  if(abs(p2*r*r/pc)>.05_wp.or.g1*r*r>.05_wp)then
   call fail_status(status,NS_OUT_OF_DOMAIN,'numerics','start radius outside small central-series domain');return
  endif
  m=r**3*(m3+m5*r*r);p=pc+r*r*(p2+p4*r*r)
  nu_offset=r*r*(g1/2._wp+g3*r*r/4._wp)
  if(.not.all(ieee_is_finite([m,p,nu_offset])).or.min(m,p)<=0._wp)then
   call fail_status(status,NS_NONFINITE_RESULT,'numerics','central-series evaluation failed');return
  endif
  call clear_status(status)
 end subroutine
 subroutine cowling_center_second(ec,pc,chi,nu_c,omega,r,w_factor,v_factor,status)
  real(wp),intent(in)::ec,pc,chi,nu_c,omega,r
  real(wp),intent(out)::w_factor,v_factor
  type(status_t),intent(out)::status
  real(wp)::g1,lambda2,kc,w2,v2,pi
  w_factor=0._wp;v_factor=0._wp
  call fail_status(status,NS_INVALID_INPUT,'physics','invalid ell=2 Cowling series input')
  if(.not.all(ieee_is_finite([ec,pc,chi,nu_c,omega,r])))return
  if(min(ec,pc,chi,r)<=0._wp.or.chi>1._wp.or.omega<0._wp)return
  pi=acos(-1._wp)
  g1=4._wp*pi*g_si/(c_si**4)*(ec/3._wp+pc)
  lambda2=4._wp*pi*g_si*ec/(3._wp*c_si**4)
  kc=(omega/c_si)**2*exp(-2._wp*nu_c)
  w2=(4._wp*(g1-kc/2._wp)/chi+18._wp*lambda2+6._wp*g1)/14._wp
  v2=(2._wp*g1+2._wp*lambda2+2._wp*w2)/4._wp
  if(.not.all(ieee_is_finite([w2,v2])))return
  if(max(abs(w2*r*r),abs(v2*r*r))>.05_wp)then
   call fail_status(status,NS_OUT_OF_DOMAIN,'numerics','oscillatory series outside small-start domain');return
  endif
  w_factor=1._wp+w2*r*r;v_factor=1._wp+v2*r*r
  call clear_status(status)
 end subroutine
end module
