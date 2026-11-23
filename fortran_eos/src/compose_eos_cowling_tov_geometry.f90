! Phase-local TOV dense output and batched Cowling geometry.
! Signature (-,+,+,+); r,m are geometrized km, P is km^-2, h is
! dimensionless. Original EOS sampling and phase identity are upstream.
! These kernels have no saved state and do not enter DVODE.
module compose_eos_cowling_tov_geometry
  use iso_c_binding, only: c_double,c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: tov_native_phase_values,tov_native_geometry,tov_native_pchip_slopes
  interface
    ! Python float ** 3 calls libm pow; keep its rounding, rather than
    ! replacing the cube with two rounded floating multiplications.
    function c_pow(x,power) bind(c,name="pow") result(y)
      import :: c_double
      real(c_double),value::x,power
      real(c_double)::y
    end function
    ! libm log1p preserves weak-field terms lost by LOG(1+x).
    function c_log1p(x) bind(c,name="log1p") result(y)
      import :: c_double
      real(c_double),value::x
      real(c_double)::y
    end function
  end interface
contains
  ! Preserve the SciPy 1.17 PCHIP/CubicHermiteSpline arithmetic sequence
  ! used by the source-locked TOV baseline. Generic Cowling PCHIP stays
  ! unchanged; mathematical equivalence alone does not fix adaptive histories.
  subroutine tov_native_pchip_slopes(n,s,y,d,at_knots,status) bind(c)
    integer(c_int),value::n,at_knots
    real(c_double),intent(in)::s(n),y(2,n)
    real(c_double),intent(out)::d(2,n)
    integer(c_int),intent(out)::status
    real(c_double)::hl,hr,dl,dr,w1,w2,whmean,t,c0,c1,res
    integer::i,j
    status=1
    d=0._c_double
    if(n<2.or.at_knots<0.or.at_knots>1)return
    if(.not.all(ieee_is_finite(s)).or..not.all(ieee_is_finite(y)))return
    if(any(s(2:n)<=s(1:n-1)))return
    if(n==2)then
      d(:,1)=(y(:,2)-y(:,1))/(s(2)-s(1))
      d(:,2)=d(:,1)
    else
      do j=1,2
        hl=s(2)-s(1);hr=s(3)-s(2)
        dl=(y(j,2)-y(j,1))/hl;dr=(y(j,3)-y(j,2))/hr
        d(j,1)=tov_end_slope(hl,hr,dl,dr)
        do i=2,n-1
          hl=s(i)-s(i-1);hr=s(i+1)-s(i)
          dl=(y(j,i)-y(j,i-1))/hl;dr=(y(j,i+1)-y(j,i))/hr
          if(dl==0._c_double.or.dr==0._c_double)then
            d(j,i)=0._c_double
          else if(sign(1._c_double,dl)/=sign(1._c_double,dr))then
            d(j,i)=0._c_double
          else
            w1=2._c_double*hr+hl;w2=hr+2._c_double*hl
            whmean=(w1/dl+w2/dr)/(w1+w2)
            d(j,i)=1._c_double/whmean
          endif
        enddo
        hl=s(n)-s(n-1);hr=s(n-1)-s(n-2)
        dl=(y(j,n)-y(j,n-1))/hl;dr=(y(j,n-1)-y(j,n-2))/hr
        d(j,n)=tov_end_slope(hl,hr,dl,dr)
      enddo
    endif
    if(.not.all(ieee_is_finite(d)))then
      d=0._c_double;status=2;return
    endif
    if(at_knots==1)then
      ! Interior knots select the interval on their right, giving d there.
      ! The final knot is evaluated at the previous interval's full width.
      hl=s(n)-s(n-1)
      do j=1,2
        dl=(y(j,n)-y(j,n-1))/hl
        t=(d(j,n-1)+d(j,n)-2._c_double*dl)/hl
        c0=t/hl;c1=(dl-d(j,n-1))/hl-t
        res=d(j,n-1)
        res=res+c1*hl*2._c_double
        res=res+c0*(hl*hl)*3._c_double
        d(j,n)=res
      enddo
    endif
    if(.not.all(ieee_is_finite(d)))then
      d=0._c_double;status=2;return
    endif
    status=0
  end subroutine

  pure real(c_double) function tov_end_slope(h0,h1,m0,m1) result(d)
    real(c_double),intent(in)::h0,h1,m0,m1
    d=((2._c_double*h0+h1)*m0-h0*m1)/(h0+h1)
    if(m0==0._c_double)then
      d=0._c_double
    else if(d==0._c_double)then
      d=0._c_double
    else if(sign(1._c_double,d)/=sign(1._c_double,m0))then
      d=0._c_double
    else if(sign(1._c_double,m0)/=sign(1._c_double,m1).and.abs(d)>3._c_double*abs(m0))then
      d=3._c_double*m0
    endif
  end function

  subroutine tov_native_phase_values(n,hstart,s,y,d,nquery,h,values,status) bind(c)
    integer(c_int),value::n,nquery
    real(c_double),value::hstart
    real(c_double),intent(in)::s(n),y(2,n),d(2,n),h(nquery)
    real(c_double),intent(out)::values(2,nquery)
    integer(c_int),intent(out)::status
    real(c_double)::q,width,dx,delta,t,c0,c1,z,res
    integer::i,j,lo,hi,mid
    status=1
    values=0._c_double
    if(n<2.or.nquery<1)return
    if(.not.ieee_is_finite(hstart).or..not.all(ieee_is_finite(h)))return
    if(.not.all(ieee_is_finite(s)).or..not.all(ieee_is_finite(y)).or..not.all(ieee_is_finite(d)))return
    if(hstart<0._c_double.or.any(h<0._c_double).or.any(h>hstart))return
    if(any(s(2:n)<=s(1:n-1)))return
    lo=1
    do i=1,nquery
      q=sqrt(max(hstart-h(i),0._c_double))
      if(q<s(1)-1.e-12_c_double.or.q>s(n)+1.e-12_c_double)then
        values=0._c_double
        return
      endif
      q=min(s(n),max(s(1),q))
      if(q>=s(n))then
        lo=n-1
      else if(.not.(s(lo)<=q.and.q<s(lo+1)))then
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
      width=s(lo+1)-s(lo)
      dx=q-s(lo)
      do j=1,2
        delta=(y(j,lo+1)-y(j,lo))/width
        t=(d(j,lo)+d(j,lo+1)-2._c_double*delta)/width
        c0=t/width;c1=(delta-d(j,lo))/width-t
        ! SciPy _ppoly.evaluate_poly1 sums ascending powers, not Horner.
        res=y(j,lo)
        res=res+d(j,lo)*dx
        z=dx*dx
        res=res+c1*z
        z=z*dx
        res=res+c0*z
        values(j,i)=res
      enddo
    enddo
    if(.not.all(ieee_is_finite(values)))then
      values=0._c_double;status=2;return
    endif
    status=0
  end subroutine

  subroutine tov_native_geometry(n,h,h_surface,m_total,r_surface,r_m, &
       pressure_geom,cs2,values,status) bind(c)
    integer(c_int),value::n
    real(c_double),value::h_surface,m_total,r_surface
    real(c_double),intent(in)::h(n),r_m(2,n),pressure_geom(n),cs2(n)
    real(c_double),intent(out)::values(5,n)
    integer(c_int),intent(out)::status
    real(c_double)::surface_phi,r,m,pi
    integer::i
    status=1
    values=0._c_double
    if(n<1)return
    if(.not.all(ieee_is_finite([h_surface,m_total,r_surface])))return
    if(.not.all(ieee_is_finite(h)).or..not.all(ieee_is_finite(r_m)))return
    if(.not.all(ieee_is_finite(pressure_geom)).or..not.all(ieee_is_finite(cs2)))return
    if(h_surface<0._c_double.or.r_surface<=0._c_double.or.m_total<0._c_double)return
    if(m_total>=.5_c_double*r_surface)return
    if(any(h<h_surface).or.any(r_m(1,:)<=0._c_double).or.any(r_m(2,:)<0._c_double))return
    if(any(r_m(2,:)>=.5_c_double*r_m(1,:)))return
    if(any(pressure_geom<0._c_double).or.any(cs2<=0._c_double))return
    surface_phi=.5_c_double*c_log1p(-2._c_double*m_total/r_surface)
    pi=acos(-1._c_double)
    do i=1,n
      r=r_m(1,i);m=r_m(2,i)
      values(1,i)=r
      values(2,i)=surface_phi-h(i)+h_surface
      values(3,i)=-.5_c_double*c_log1p(-2._c_double*m/r)
      values(4,i)=(m+4._c_double*pi*c_pow(r,3._c_double)*pressure_geom(i))/(r*(r-2._c_double*m))
      values(5,i)=cs2(i)
    enddo
    if(.not.all(ieee_is_finite(values)))then
      values=0._c_double
      status=2
      return
    endif
    status=0
  end subroutine
end module compose_eos_cowling_tov_geometry
