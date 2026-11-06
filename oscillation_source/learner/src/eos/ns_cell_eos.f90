! Local-coordinate evaluation of an existing admitted cubic P(mu).
! No new fit, no altered native knots, and no interpolation of chi.
module ns_cell_eos
 use ns_precision, only: wp
 use ns_status
 use ns_eos_types
 use ns_matching_types, only: mu_bridge_t
 use ns_table_eos, only: table_eos_t
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public :: cell_state,cell_pressure_coordinate,table_pressure_coordinate
contains
 subroutine cell_state(cell,t,s,dn_dmu,status)
  type(mu_bridge_t),intent(in)::cell
  real(wp),intent(in)::t
  type(eos_state_t),intent(out)::s
  real(wp),intent(out)::dn_dmu
  type(status_t),intent(out)::status
  real(wp)::h,mu,p,n,dn,e,chi
  call clear_eos_state(s);dn_dmu=0._wp
  call fail_status(status,NS_OUT_OF_DOMAIN,'physics','local EOS coordinate must lie in [0,1]')
  if(.not.cell%ready.or..not.ieee_is_finite(t))return
  if(t<0._wp.or.t>1._wp)return
  h=cell%mu_high_j-cell%mu_low_j
  mu=cell%mu_low_j+h*t
  p=cell%a0+t*(cell%a1+t*(cell%a2+t*cell%a3))
  n=(cell%a1+t*(2._wp*cell%a2+3._wp*t*cell%a3))/h
  dn=(2._wp*cell%a2+6._wp*t*cell%a3)/(h*h)
  e=mu*n-p;chi=n/(mu*dn)
  call fail_status(status,NS_MATCHING_INCONSISTENT,'physics','nonpositive or acausal local thermodynamics')
  if(.not.all(ieee_is_finite([p,n,dn,e,mu,chi])))return
  if(min(p,n,dn,e,mu,chi)<=0._wp.or.chi>1._wp+32._wp*epsilon(1._wp))return
  ! chi is not clipped: a roundoff-sized excess in an admitted analytic cell is returned.
  s=eos_state_t(p,e,n,mu,chi,SIDE_REGULAR,.true.);dn_dmu=dn
  call clear_status(status)
 end subroutine
 subroutine cell_pressure_coordinate(cell,p,t,status)
  type(mu_bridge_t),intent(in)::cell
  real(wp),intent(in)::p
  real(wp),intent(out)::t
  type(status_t),intent(out)::status
  real(wp)::lo,hi,f,df,tn,dp,tol
  integer::j
  t=0._wp
  call fail_status(status,NS_OUT_OF_DOMAIN,'physics','pressure outside the selected native cell')
  if(.not.cell%ready.or..not.ieee_is_finite(p))return
  if(p<cell%p_low_pa.or.p>cell%p_high_pa)return
  if(p==cell%p_low_pa)then;call clear_status(status);return;endif
  if(p==cell%p_high_pa)then;t=1._wp;call clear_status(status);return;endif
  lo=0._wp;hi=1._wp;dp=p-cell%a0
  t=dp/(cell%p_high_pa-cell%p_low_pa)
  tol=16._wp*epsilon(1._wp)*max(p,cell%p_low_pa)
  do j=1,100
   ! Subtract the pressure anchor before solving, not two almost equal chemical potentials.
   f=t*(cell%a1+t*(cell%a2+t*cell%a3))-dp
   if(abs(f)<=tol)then;call clear_status(status);return;endif
   if(f>0._wp)then;hi=t;else;lo=t;endif
   df=cell%a1+t*(2._wp*cell%a2+3._wp*t*cell%a3)
   tn=t-f/df
   if(tn<=lo.or.tn>=hi.or..not.ieee_is_finite(tn))tn=lo+0.5_wp*(hi-lo)
   if(tn==t)exit
   t=tn
  enddo
  call fail_status(status,NS_STAGNATION,'numerics','local EOS inverse did not close pressure residual')
 end subroutine
 subroutine table_pressure_coordinate(model,p,i,t,s,status)
  type(table_eos_t),intent(in)::model
  real(wp),intent(in)::p
  integer,intent(out)::i
  real(wp),intent(out)::t
  type(eos_state_t),intent(out)::s
  type(status_t),intent(out)::status
  integer::lo,hi,mid
  real(wp)::dn
  i=0;t=0._wp;call clear_eos_state(s)
  call fail_status(status,NS_OUT_OF_DOMAIN,'physics','pressure outside admitted table')
  if(.not.model%ready.or..not.ieee_is_finite(p))return
  if(p<model%raw%pressure_pa(1).or.p>model%raw%pressure_pa(model%raw%n))return
  lo=1;hi=model%raw%n-1
  do while(lo<hi)
   mid=(lo+hi)/2
   if(p<=model%raw%pressure_pa(mid+1))then;hi=mid;else;lo=mid+1;endif
  enddo
  i=lo
  call cell_pressure_coordinate(model%cells(i),p,t,status)
  if(.not.is_ok(status))return
  call cell_state(model%cells(i),t,s,dn,status)
 end subroutine
end module
