! C16 Fortran local mode-shape overlap diagnostic.
! The radial Cowling W displacement is compared phase-by-phase in a
! normalized radial coordinate u \\in [0,1] for each Maxwell phase.
! The phase boundary is never crossed by interpolation. The absolute
! normalized inner products are invariant to an overall eigenmode sign.
!
! WARNING: This is a descriptive cross-star shape diagnostic, NOT a
! covariant mode-inner-product / avoided-crossing theorem. Physical
! eigenmode admission still requires the original Maxwell junction,
! matched determinant and phasewise node counts.
module compose_eos_cowling_mode_overlap
  use iso_c_binding,only:c_int,c_double
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  integer,parameter::ngrid=64
  public::cowling_native_mode_shape_overlap
contains
  pure real(c_double) function interp_linear_phase(r,w,i0,i1,u) result(value)
    real(c_double),intent(in)::r(:),w(:),u
    integer,intent(in)::i0,i1
    real(c_double)::x,dx
    integer::lo,hi,mid
    if(u<=0._c_double)then
      value=w(i0);return
    endif
    if(u>=1._c_double)then
      value=w(i1);return
    endif
    x=r(i0)+u*(r(i1)-r(i0))
    lo=i0;hi=i1
    do while(hi-lo>1)
      mid=(hi+lo)/2
      if(r(mid)<=x)then
        lo=mid
      else
        hi=mid
      endif
    enddo
    dx=r(hi)-r(lo)
    if(dx<=0._c_double)then
      value=w(hi)
    else
      value=w(lo)+(w(hi)-w(lo))*(x-r(lo))/dx
    endif
  end function

  subroutine one_phase(r1,w1,a,b,r2,w2,c,d,norm1,norm2,cross)
    real(c_double),intent(in)::r1(:),w1(:),r2(:),w2(:)
    integer,intent(in)::a,b,c,d
    real(c_double),intent(out)::norm1,norm2,cross
    real(c_double)::u,t,v,weight
    integer::i
    norm1=0._c_double;norm2=0._c_double;cross=0._c_double
    do i=0,ngrid
      u=real(i,c_double)/real(ngrid,c_double)
      t=interp_linear_phase(r1,w1,a,b,u)
      v=interp_linear_phase(r2,w2,c,d,u)
      weight=1._c_double
      if(i==0.or.i==ngrid)weight=.5_c_double
      norm1=norm1+weight*t*t
      norm2=norm2+weight*v*v
      cross=cross+weight*t*v
    enddo
    norm1=norm1/real(ngrid,c_double)
    norm2=norm2/real(ngrid,c_double)
    cross=cross/real(ngrid,c_double)
  end subroutine

  subroutine cowling_native_mode_shape_overlap(n1,phase1,r1,w1, &
       n2,phase2,r2,w2,total,inner,outer,status) bind(c)
    integer(c_int),value::n1,n2
    integer(c_int),intent(in)::phase1(n1),phase2(n2)
    real(c_double),intent(in)::r1(n1),w1(n1),r2(n2),w2(n2)
    real(c_double),intent(out)::total,inner,outer
    integer(c_int),intent(out)::status
    integer::a,b,c,d,p
    real(c_double)::aa(2),bb(2),ab(2)
    total=0._c_double;inner=0._c_double;outer=0._c_double;status=1
    if(n1<6.or.n2<6.or.n1>100000.or.n2>100000)return
    if(phase1(1)/=0.or.phase2(1)/=0)return
    if(phase1(n1)/=1.or.phase2(n2)/=1)return
    if(any(phase1(2:n1)<phase1(1:n1-1)).or. &
       any(phase2(2:n2)<phase2(1:n2-1)))return
    if(any(phase1<0).or.any(phase1>1).or.any(phase2<0).or.any(phase2>1))return
    if(.not.all(ieee_is_finite(r1)).or..not.all(ieee_is_finite(r2)))return
    if(.not.all(ieee_is_finite(w1)).or..not.all(ieee_is_finite(w2)))return
    if(any(r1<0._c_double).or.any(r2<0._c_double))return
    if(any(r1(2:n1)<r1(1:n1-1)).or.any(r2(2:n2)<r2(1:n2-1)))return
    do p=0,1
      a=findloc(phase1,int(p,c_int),dim=1)
      c=findloc(phase2,int(p,c_int),dim=1)
      b=a+count(phase1==p)-1
      d=c+count(phase2==p)-1
      if(a<=0.or.c<=0.or.b-a<2.or.d-c<2)return
      if(r1(b)<=r1(a).or.r2(d)<=r2(c))return
      call one_phase(r1,w1,a,b,r2,w2,c,d,aa(p+1),bb(p+1),ab(p+1))
    enddo
    if(any(aa<=0._c_double).or.any(bb<=0._c_double))then
      status=2;return
    endif
    inner=min(1._c_double,abs(ab(1))/sqrt(aa(1)*bb(1)))
    outer=min(1._c_double,abs(ab(2))/sqrt(aa(2)*bb(2)))
    total=min(1._c_double,abs(sum(ab))/sqrt(sum(aa)*sum(bb)))
    if(.not.all(ieee_is_finite([inner,outer,total])))then
      inner=0._c_double;outer=0._c_double;total=0._c_double
      status=3;return
    endif
    status=0
  end subroutine
end module compose_eos_cowling_mode_overlap
