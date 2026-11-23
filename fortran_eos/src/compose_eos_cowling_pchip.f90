! Source-domain one-sided PCHIP slope construction and field evaluation.
! SciPy PchipInterpolator slope definitions (Fritsch-Butland weighted
! harmonic mean, one-sided three-point endpoint limiter) migrated to
! Fortran. Geometry is phase-local, never interpolated across Maxwell.
! Units preserved: independent s=sqrt(h_start-h), dimensionless.
!
! Designed for 5 sampled background fields per phase: radius, Phi,
! lambda, Phi_prime and a PCHIP-compatible spare c_s² row. The source
! exact Cowling RHS overrides that last field analytically.
module compose_eos_cowling_pchip
  use iso_c_binding, only: c_double,c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: cowling_native_pchip_slopes,cowling_native_pchip_values
  public :: cowling_native_pchip_slopes_fields,cowling_native_pchip_values_fields
contains
  pure real(c_double) function end_slope(h0,h1,del0,del1) result(m)
    real(c_double),intent(in)::h0,h1,del0,del1
    m=((2._c_double*h0+h1)*del0-h0*del1)/(h0+h1)
    if(sign(1._c_double,m)/=sign(1._c_double,del0))then
      m=0._c_double
    else if(sign(1._c_double,del0)/=sign(1._c_double,del1) &
            .and.abs(m)>3._c_double*abs(del0))then
      m=3._c_double*del0
    endif
  end function

  ! Fortran column-major values(5,n), corresponding to Python y(n,5).
  ! No cancellation of physical EOS sharp jumps: the caller supplies
  ! exactly one phase's samples and validates the original source.
  subroutine cowling_native_pchip_slopes(n,s,y,d,status) bind(c)
    integer(c_int),value::n
    real(c_double),intent(in)::s(n),y(5,n)
    real(c_double),intent(out)::d(5,n)
    integer(c_int),intent(out)::status
    call cowling_native_pchip_slopes_fields(5_c_int,n,s,y,d,status)
  end subroutine

  ! Compact field-count ABI: two fields for TOV, five for Cowling.
  ! No padding, saved state, or allocation is needed by the slope kernel.
  subroutine cowling_native_pchip_slopes_fields(nfields,n,s,y,d,status) bind(c)
    integer(c_int),value::nfields,n
    real(c_double),intent(in)::s(n),y(nfields,n)
    real(c_double),intent(out)::d(nfields,n)
    integer(c_int),intent(out)::status
    real(c_double)::hl,hr,dl,dr,w1,w2
    integer::i,j
    status=1
    d=0._c_double
    if(nfields<1.or.n<2)return
    if(.not.all(ieee_is_finite(s)).or..not.all(ieee_is_finite(y)))return
    if(any(s(2:n)<=s(1:n-1)))return
    if(n==2)then
      d(:,1)=(y(:,2)-y(:,1))/(s(2)-s(1))
      d(:,2)=d(:,1)
    else
      do j=1,nfields
        hl=s(2)-s(1)
        hr=s(3)-s(2)
        dl=(y(j,2)-y(j,1))/hl
        dr=(y(j,3)-y(j,2))/hr
        d(j,1)=end_slope(hl,hr,dl,dr)
        do i=2,n-1
          hl=s(i)-s(i-1)
          hr=s(i+1)-s(i)
          dl=(y(j,i)-y(j,i-1))/hl
          dr=(y(j,i+1)-y(j,i))/hr
          if(dl==0._c_double.or.dr==0._c_double)then
            d(j,i)=0._c_double
          else if(sign(1._c_double,dl)/=sign(1._c_double,dr))then
            d(j,i)=0._c_double
          else
            w1=2._c_double*hr+hl
            w2=hr+2._c_double*hl
            d(j,i)=(w1+w2)/(w1/dl+w2/dr)
          endif
        enddo
        hl=s(n)-s(n-1)
        hr=s(n-1)-s(n-2)
        dl=(y(j,n)-y(j,n-1))/hl
        dr=(y(j,n-1)-y(j,n-2))/hr
        d(j,n)=end_slope(hl,hr,dl,dr)
      enddo
    endif
    if(.not.all(ieee_is_finite(d)))then
      d=0._c_double
      status=2
      return
    endif
    status=0
  end subroutine

  ! Independent C ABI profile-field interpolator. O(nquery) for nearby
  ! ordered points, binary fallback otherwise; never extrapolates.
  subroutine cowling_native_pchip_values(n,s,y,d,nquery,x,v,status) bind(c)
    integer(c_int),value::n,nquery
    real(c_double),intent(in)::s(n),y(5,n),d(5,n),x(nquery)
    real(c_double),intent(out)::v(5,nquery)
    integer(c_int),intent(out)::status
    call cowling_native_pchip_values_fields(5_c_int,n,s,y,d,nquery,x,v,status)
  end subroutine

  ! Batch evaluation validates one complete phase table per call.
  ! Nearby queries reuse an interval; arbitrary order uses binary fallback.
  subroutine cowling_native_pchip_values_fields(nfields,n,s,y,d,nquery,x,v,status) bind(c)
    integer(c_int),value::nfields,n,nquery
    real(c_double),intent(in)::s(n),y(nfields,n),d(nfields,n),x(nquery)
    real(c_double),intent(out)::v(nfields,nquery)
    integer(c_int),intent(out)::status
    real(c_double)::dx,z,z2,z3,q
    integer::i,lo,hi,mid
    status=1
    v=0._c_double
    if(nfields<1.or.n<2.or.nquery<1)return
    if(.not.all(ieee_is_finite(s)).or..not.all(ieee_is_finite(y)) &
       .or..not.all(ieee_is_finite(d)).or..not.all(ieee_is_finite(x)))return
    if(any(s(2:n)<=s(1:n-1)))return
    if(any(x<s(1)-1.e-12_c_double).or.any(x>s(n)+1.e-12_c_double))return
    lo=1
    do i=1,nquery
      q=min(s(n),max(s(1),x(i)))
      if(q<=s(1))then
        v(:,i)=y(:,1)
      else if(q>=s(n))then
        v(:,i)=y(:,n)
      else
        if(.not.(s(lo)<=q.and.q<s(lo+1)))then
          if(lo<n-1)then
            if(s(lo+1)<=q.and.q<s(lo+2))lo=lo+1
          endif
          if(lo>1)then
            if(s(lo-1)<=q.and.q<s(lo))lo=lo-1
          endif
          if(.not.(s(lo)<=q.and.q<s(lo+1)))then
            lo=1;hi=n
            do while(hi-lo>1)
              mid=(lo+hi)/2
              if(s(mid)<=q)then
                lo=mid
              else
                hi=mid
              endif
            enddo
          endif
        endif
        dx=s(lo+1)-s(lo)
        z=(q-s(lo))/dx
        z2=z*z
        z3=z2*z
        v(:,i)=(2*z3-3*z2+1)*y(:,lo)+(z3-2*z2+z)*dx*d(:,lo) &
              +(-2*z3+3*z2)*y(:,lo+1)+(z3-z2)*dx*d(:,lo+1)
      endif
    enddo
    if(.not.all(ieee_is_finite(v)))then
      v=0._c_double;status=2;return
    endif
    status=0
  end subroutine
end module compose_eos_cowling_pchip
