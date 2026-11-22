! C18: phasewise, finite-scale-stable Cowling W-shape overlap.
!
! W profiles are compared in dimensionless areal radius within each
! Maxwell phase. This is a descriptive diagnostic, NOT a covariant
! pulsation inner product or an avoided-crossing certification.
! Only the jump between phase 0 and phase 1 may repeat the areal radius.
! No DVODE state, COMMON blocks, saved variables, or Python numerics.
module compose_eos_cowling_mode_overlap
  use iso_c_binding, only:c_int,c_double
  use, intrinsic::ieee_arithmetic, only:ieee_is_finite
  implicit none
  private
  integer,parameter::ngrid=64
  public::cowling_native_mode_shape_overlap
contains
  ! u increases monotonically over the fixed integration grid.
  ! For short physical mode profiles, a moving interval cursor avoids
  ! 2*ngrid binary searches per phase. For large arrays, binary lookup
  ! retains logarithmic cost instead of scanning every input interval.
  real(c_double) function interp_scaled(r,w,i0,i1,u,scale,hint) result(y)
    real(c_double),intent(in)::r(:),w(:),u,scale
    integer,intent(in)::i0,i1
    integer,intent(inout)::hint
    real(c_double)::x,t
    integer::lo,hi,mid
    if(u<=0._c_double)then
      y=w(i0)/scale
      return
    endif
    if(u>=1._c_double)then
      y=w(i1)/scale
      return
    endif
    x=r(i0)+u*(r(i1)-r(i0))
    if(i1-i0<=512)then
      lo=hint
      do while(lo<i1-1)
        if(r(lo+1)>=x)exit
        lo=lo+1
      enddo
      hint=lo
      hi=lo+1
    else
      lo=i0;hi=i1
      do while(hi-lo>1)
        mid=(hi+lo)/2
        if(r(mid)<=x)then
          lo=mid
        else
          hi=mid
        endif
      enddo
    endif
    t=(x-r(lo))/(r(hi)-r(lo))
    ! Each endpoint is normalized BEFORE subtraction to avoid overflow.
    y=(w(lo)/scale)*(1._c_double-t)+(w(hi)/scale)*t
  end function

  subroutine one_phase(r1,w1,a,b,r2,w2,c,d,aa,bb,ab,s1,s2)
    real(c_double),intent(in)::r1(:),w1(:),r2(:),w2(:)
    integer,intent(in)::a,b,c,d
    real(c_double),intent(out)::aa,bb,ab,s1,s2
    real(c_double)::u,t,v,weight
    integer::i,k1,k2
    aa=0._c_double;bb=0._c_double;ab=0._c_double
    s1=maxval(abs(w1(a:b)));s2=maxval(abs(w2(c:d)))
    if(s1<=0._c_double.or.s2<=0._c_double)return
    k1=a;k2=c
    do i=0,ngrid
      u=real(i,c_double)/real(ngrid,c_double)
      t=interp_scaled(r1,w1,a,b,u,s1,k1)
      v=interp_scaled(r2,w2,c,d,u,s2,k2)
      weight=1._c_double
      if(i==0.or.i==ngrid)weight=.5_c_double
      aa=aa+weight*t*t
      bb=bb+weight*v*v
      ab=ab+weight*t*v
    enddo
    aa=aa/real(ngrid,c_double)
    bb=bb/real(ngrid,c_double)
    ab=ab/real(ngrid,c_double)
  end subroutine

  subroutine cowling_native_mode_shape_overlap(n1,phase1,r1,w1, &
       n2,phase2,r2,w2,total,inner,outer,status) bind(c)
    integer(c_int),value::n1,n2
    integer(c_int),intent(in)::phase1(n1),phase2(n2)
    real(c_double),intent(in)::r1(n1),w1(n1),r2(n2),w2(n2)
    real(c_double),intent(out)::total,inner,outer
    integer(c_int),intent(out)::status
    integer::m1,m2,p,a,b,c,d
    real(c_double)::aa(2),bb(2),ab(2),s1(2),s2(2),g1,g2
    real(c_double)::q1(2),q2(2),an,bn,cn
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
    ! The Maxwell junction may have r_high=r_low, but neither phase may
    ! have non-increasing internal radius values.
    if(any(r1(2:n1)<r1(1:n1-1)).or.any(r2(2:n2)<r2(1:n2-1)))return
    m1=count(phase1==0_c_int);m2=count(phase2==0_c_int)
    if(m1<3.or.n1-m1<3.or.m2<3.or.n2-m2<3)return
    if(any(r1(2:m1)<=r1(1:m1-1)))return
    if(any(r1(m1+2:n1)<=r1(m1+1:n1-1)))return
    if(any(r2(2:m2)<=r2(1:m2-1)))return
    if(any(r2(m2+2:n2)<=r2(m2+1:n2-1)))return
    do p=1,2
      if(p==1)then
        a=1;b=m1;c=1;d=m2
      else
        a=m1+1;b=n1;c=m2+1;d=n2
      endif
      call one_phase(r1,w1,a,b,r2,w2,c,d,aa(p),bb(p),ab(p),s1(p),s2(p))
    enddo
    if(any(aa<=0._c_double).or.any(bb<=0._c_double))then
      status=2;return
    endif
    inner=min(1._c_double,abs(ab(1))/sqrt(aa(1)*bb(1)))
    outer=min(1._c_double,abs(ab(2))/sqrt(aa(2)*bb(2)))
    ! Phase-local scales retain the true global relative phase weighting.
    ! Dividing by global per-star maxima avoids overflow/underflow even
    ! when the two physical W profiles have vastly different amplitudes.
    g1=maxval(s1);g2=maxval(s2)
    q1=s1/g1;q2=s2/g2
    an=sum(aa*q1*q1);bn=sum(bb*q2*q2)
    cn=sum(ab*q1*q2)
    if(an<=0._c_double.or.bn<=0._c_double)then
      status=2;inner=0._c_double;outer=0._c_double;return
    endif
    total=min(1._c_double,abs(cn)/sqrt(an*bn))
    if(.not.all(ieee_is_finite([inner,outer,total])))then
      inner=0._c_double;outer=0._c_double;total=0._c_double
      status=3;return
    endif
    status=0
  end subroutine
end module compose_eos_cowling_mode_overlap
