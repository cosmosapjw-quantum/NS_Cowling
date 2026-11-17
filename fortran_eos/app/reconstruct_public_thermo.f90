! Minimal energy-primary public-table reconstruction helper.
! Units: n_B [fm^-3], F-m_n*c^2 [MeV].  No physical constant is changed.
! Input: number of rows, then n_B and source binding energy per baryon.
! Output: dF/dln(n_B), d2F/dln(n_B)^2, n_B*dF/dn_B from an independent
! cubic-in-density control.  Source energy knots are never changed.
!
! The two C2 splines use not-a-knot boundary conditions. The tridiagonal
! equations follow from equality of second derivatives at internal knots
! and equality of third derivatives in the first/last two intervals.
program reconstruct_public_thermo
  use, intrinsic :: iso_fortran_env, only: real64, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: wp=real64
  integer :: n,i,j,ios,u,v
  character(len=4096) :: input,output
  real(wp),allocatable :: nb(:),binding(:),x(:),sx(:),sn(:),second(:)
  real(wp) :: h,c,d
  if(command_argument_count()/=2) error stop 'usage: reconstruct_public_thermo INPUT OUTPUT'
  call get_command_argument(1,input)
  call get_command_argument(2,output)
  open(newunit=u,file=trim(input),status='old',action='read',iostat=ios)
  if(ios/=0) error stop 'cannot open reconstruction input'
  read(u,*,iostat=ios)n
  if(ios/=0.or.n<4) error stop 'energy-primary C2 reconstruction requires >=4 rows'
  allocate(nb(n),binding(n),x(n),sx(n),sn(n),second(n))
  do i=1,n
    read(u,*,iostat=ios)nb(i),binding(i)
    if(ios/=0) error stop 'invalid reconstruction input row'
  enddo
  close(u)
  if(.not.all(ieee_is_finite(nb)).or..not.all(ieee_is_finite(binding))) &
    error stop 'nonfinite source values'
  if(any(nb<=0._wp).or.any(nb(2:)<=nb(:n-1))) &
    error stop 'nonpositive or unordered source densities'
  x=log(nb/nb(1))
  ! Removing a constant before differentiation avoids cancellation of the
  ! rest-energy scale. It has no effect on the physical energy reference.
  binding=binding-binding(1)
  call not_a_knot_slopes(x,binding,sx)
  call not_a_knot_slopes(nb,binding,sn)
  do i=1,n
    j=min(i,n-1)
    h=x(j+1)-x(j)
    c=3._wp*(binding(j+1)-binding(j))-h*(2._wp*sx(j)+sx(j+1))
    d=2._wp*(binding(j)-binding(j+1))+h*(sx(j)+sx(j+1))
    second(i)=2._wp*c/(h*h)
    if(i==n)second(i)=second(i)+6._wp*d/(h*h)
  enddo
  if(.not.all(ieee_is_finite(sx)).or..not.all(ieee_is_finite(sn)).or. &
     .not.all(ieee_is_finite(second))) error stop 'nonfinite spline output'
  open(newunit=v,file=trim(output),status='replace',action='write',iostat=ios)
  if(ios/=0) error stop 'cannot open reconstruction output'
  do i=1,n
    write(v,'(3(es26.17e3,1x))')sx(i),second(i),nb(i)*sn(i)
  enddo
  close(v)
contains
  subroutine not_a_knot_slopes(grid,values,slope)
    real(wp),intent(in) :: grid(:),values(:)
    real(wp),intent(out) :: slope(:)
    real(wp),allocatable :: step(:),secant(:),lower(:),diagonal(:),upper(:),rhs(:)
    real(wp) :: span,factor
    integer :: k,m
    m=size(grid)
    allocate(step(m-1),secant(m-1),lower(m),diagonal(m),upper(m),rhs(m))
    step=grid(2:)-grid(:m-1)
    secant=(values(2:)-values(:m-1))/step
    lower=0._wp;upper=0._wp;diagonal=0._wp;rhs=0._wp
    do k=2,m-1
      lower(k)=step(k)
      diagonal(k)=2._wp*(step(k-1)+step(k))
      upper(k)=step(k-1)
      rhs(k)=3._wp*(step(k)*secant(k-1)+step(k-1)*secant(k))
    enddo
    span=step(1)+step(2)
    diagonal(1)=step(2);upper(1)=span
    rhs(1)=((step(1)+2._wp*span)*step(2)*secant(1)+step(1)**2*secant(2))/span
    span=step(m-2)+step(m-1)
    lower(m)=span;diagonal(m)=step(m-2)
    rhs(m)=(step(m-1)**2*secant(m-2)+(2._wp*span+step(m-1))*step(m-2)*secant(m-1))/span
    do k=2,m
      if(abs(diagonal(k-1))<=tiny(1._wp)) error stop 'singular spline system'
      factor=lower(k)/diagonal(k-1)
      diagonal(k)=diagonal(k)-factor*upper(k-1)
      rhs(k)=rhs(k)-factor*rhs(k-1)
    enddo
    if(abs(diagonal(m))<=tiny(1._wp)) error stop 'singular spline system'
    slope(m)=rhs(m)/diagonal(m)
    do k=m-1,1,-1
      slope(k)=(rhs(k)-upper(k)*slope(k+1))/diagonal(k)
    enddo
  end subroutine not_a_knot_slopes
end program reconstruct_public_thermo
