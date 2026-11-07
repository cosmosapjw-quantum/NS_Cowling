module compose_eos_verify
  use compose_eos_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  type,public :: verify_report_t
    logical :: c1=.false.,c2=.false.,c2_applicable=.true.,c3=.false.,c4=.false.,c3_available=.false.
    logical :: c5=.false.,c5_applicable=.true.,c6=.false.,c7=.false.,c8=.false.,c9=.false.
    logical :: c10=.false.,c11=.false.,overall=.false.
    integer :: derivative_points=0,derivative_masked=0,acausal_count=0,acausal_firstlaw_count=0
    real(wp) :: c1_max_q2=0._wp,c1_max_q6q7=0._wp,c1_max_mul=0._wp
    real(wp) :: c2_max_rel=0._wp,c3_max_abs=0._wp,c4_max_abs=0._wp
    real(wp) :: min_dpdnb=0._wp,min_delta_p=0._wp
    real(wp) :: c7_max_rel=0._wp,c8_max_rel=0._wp,c9_max_rel=0._wp
    real(wp) :: min_cs2=0._wp,max_cs2=0._wp,min_cs2_fl=0._wp,max_cs2_fl=0._wp
    real(wp) :: min_dmubdnb=0._wp,max_cs2_difference=0._wp
    character(len=256) :: failed_checks=''
  end type verify_report_t
  public :: check_consistency,write_verify_json
contains
  pure function relerr(a,b) result(r)
    real(wp),intent(in)::a(:),b(:)
    real(wp)::r(size(a))
    integer::i
    do i=1,size(a)
      r(i)=abs(a(i)-b(i))/max(abs(a(i)),abs(b(i)),tiny(1._wp))
    enddo
  end function relerr

  pure function derivative_nonuniform(x,y) result(d)
    real(wp),intent(in)::x(:),y(:)
    real(wp)::d(size(x)),h1,h2
    integer::n,i
    n=size(x);d=0._wp
    if(n<3)return
    h1=x(2)-x(1);h2=x(3)-x(2)
    d(1)=-(2._wp*h1+h2)/(h1*(h1+h2))*y(1)+(h1+h2)/(h1*h2)*y(2)-h1/(h2*(h1+h2))*y(3)
    do i=2,n-1
      h1=x(i)-x(i-1);h2=x(i+1)-x(i)
      d(i)=-h2/(h1*(h1+h2))*y(i-1)+(h2-h1)/(h1*h2)*y(i)+h1/(h2*(h1+h2))*y(i+1)
    enddo
    h1=x(n)-x(n-1);h2=x(n-1)-x(n-2)
    d(n)=(2._wp*h1+h2)/(h1*(h1+h2))*y(n)-(h1+h2)/(h1*h2)*y(n-1)+h1/(h2*(h1+h2))*y(n-2)
  end function derivative_nonuniform

  subroutine check_consistency(eos,report,algebraic_tolerance,max_log_step,cs2_tolerance,composition_tolerance)
    type(eos_data_t),intent(in)::eos
    type(verify_report_t),intent(out)::report
    real(wp),intent(in),optional::algebraic_tolerance,max_log_step,cs2_tolerance,composition_tolerance
    real(wp)::atol,dxmax,cstol,comptol,muqv,mulmax
    real(wp),allocatable::yl(:),dp_m(:),dmub_m(:),p_m(:),nb_m(:),df_m(:),deps_m(:),mub_m(:),cs2_m(:),csfl_m(:)
    real(wp),allocatable::p(:),f(:),fb(:),eps(:),mub(:),muq(:),mul(:),x(:)
    real(wp),allocatable::dp(:),dmub(:),df(:),deps(:),r(:),cs2(:),csfl(:)
    logical,allocatable::mask(:)
    logical::c2app,finite_primary
    integer::n,i,nbad
    character(len=8)::names(11)
    atol=2.e-6_wp;dxmax=0.25_wp;cstol=1.e-6_wp;comptol=2.e-5_wp
    if(present(algebraic_tolerance))atol=algebraic_tolerance
    if(present(max_log_step))dxmax=max_log_step
    if(present(cs2_tolerance))cstol=cs2_tolerance
    if(present(composition_tolerance))comptol=composition_tolerance
    report=verify_report_t();n=eos%n
    if(n<3.or..not.allocated(eos%nb))then
      report%failed_checks='invalid_input';return
    endif
    allocate(p(n),f(n),fb(n),eps(n),mub(n),muq(n),mul(n),x(n))
    p=eos%q(1,:)*eos%nb
    f=eos%nb*(eos%q(6,:)+1._wp)*eos%mn
    fb=(eos%q(6,:)+1._wp)*eos%mn
    eps=eos%nb*(eos%q(7,:)+1._wp)*eos%mn
    mub=(eos%q(3,:)+1._wp)*eos%mn
    muq=eos%q(4,:)*eos%mn
    mul=eos%q(5,:)*eos%mn
    x=log(eos%nb)
    finite_primary=all(ieee_is_finite(p)).and.all(ieee_is_finite(f)).and.all(ieee_is_finite(eps)).and. &
                   all(ieee_is_finite(mub)).and.all(ieee_is_finite(muq)).and.all(ieee_is_finite(mul))

    report%c1_max_q2=maxval(abs(eos%q(2,:)))
    report%c1_max_q6q7=maxval(abs(eos%q(6,:)-eos%q(7,:)))
    report%c1_max_mul=maxval(abs(mul))
    report%c1=abs(eos%temperature_mev)<=1.e-12_wp.and.report%c1_max_q2<=atol.and. &
      report%c1_max_q6q7<=atol.and.(.not.eos%beta_equilibrated.or.report%c1_max_mul<=atol*eos%mn)

    c2app=.true.
    if(eos%il==0)then
      r=relerr(mub,(p+f)/eos%nb-eos%charge_fraction_grid*muq)
    else if(maxval(abs(mul))<=1.e-8_wp*eos%mn)then
      r=relerr(mub,(p+f)/eos%nb)
    else if(eos%composition_available)then
      allocate(yl(n));yl=eos%ye+eos%ymu
      r=relerr(mub,(p+f)/eos%nb-yl*mul)
    else
      c2app=.false.;r=0._wp
    endif
    report%c2_applicable=c2app
    report%c2_max_rel=0._wp
    if(c2app)report%c2_max_rel=maxval(r)
    report%c2=c2app.and.report%c2_max_rel<=atol

    report%c3_available=eos%composition_available
    if(eos%composition_available)then
      do i=1,n
        if(eos%abar(i)>0._wp)then
          muqv=eos%zbar(i)/eos%abar(i)
        else
          muqv=0._wp
        endif
        report%c3_max_abs=max(report%c3_max_abs,abs(eos%yp(i)+muqv*eos%yh(i)-eos%ye(i)-eos%ymu(i)))
        report%c4_max_abs=max(report%c4_max_abs,abs(eos%yn(i)+eos%yp(i)+eos%yh(i)-1._wp))
      enddo
      report%c3=report%c3_max_abs<=comptol
      report%c4=report%c4_max_abs<=comptol
    else
      report%c3=.true.;report%c4=.true.
    endif
    report%c5_applicable=eos%beta_equilibrated
    report%c5=(.not.eos%beta_equilibrated).or.report%c1_max_mul<=atol*eos%mn

    report%min_delta_p=minval(p(2:)-p(:n-1))
    report%min_dpdnb=minval((p(2:)-p(:n-1))/(eos%nb(2:)-eos%nb(:n-1)))
    report%c6=report%min_delta_p>=-1.e-12_wp*maxval(abs(p))

    dp=derivative_nonuniform(eos%nb,p)
    dmub=derivative_nonuniform(eos%nb,mub)
    df=derivative_nonuniform(eos%nb,fb)
    deps=derivative_nonuniform(eos%nb,eps)
    allocate(mask(n),cs2(n),csfl(n));mask=.true.;cs2=0._wp;csfl=0._wp
    do i=1,n-1
      if(log(eos%nb(i+1)/eos%nb(i))>dxmax)mask(i:i+1)=.false.
      if(eos%composition_available)then
        if(eos%iphase(i+1)/=eos%iphase(i))mask(i:i+1)=.false.
      endif
    enddo
    do i=1,n
      if(deps(i)<=0._wp.or.mub(i)<=0._wp)mask(i)=.false.
      if(mask(i))then
        cs2(i)=dp(i)/deps(i);csfl(i)=dp(i)/mub(i)
      endif
    enddo
    report%derivative_points=count(mask);report%derivative_masked=n-report%derivative_points
    if(report%derivative_points>0)then
      dp_m=pack(dp,mask);dmub_m=pack(dmub,mask);p_m=pack(p,mask);nb_m=pack(eos%nb,mask)
      df_m=pack(df,mask);deps_m=pack(deps,mask);mub_m=pack(mub,mask)
      cs2_m=pack(cs2,mask);csfl_m=pack(csfl,mask)
      report%c7_max_rel=maxval(relerr(dp_m,nb_m*dmub_m))
      report%c8_max_rel=maxval(relerr(p_m,nb_m**2*pack(df,mask)))
      report%c9_max_rel=maxval(relerr(deps_m,mub_m))
      report%min_cs2=minval(cs2_m);report%max_cs2=maxval(cs2_m)
      report%min_cs2_fl=minval(csfl_m);report%max_cs2_fl=maxval(csfl_m)
      report%max_cs2_difference=maxval(abs(cs2_m-csfl_m))
      report%min_dmubdnb=minval(dmub_m)
      report%acausal_count=count(cs2_m>1._wp+cstol)
      report%acausal_firstlaw_count=count(csfl_m>1._wp+cstol)
      report%c7=report%c7_max_rel<=atol*100._wp
      report%c8=report%c8_max_rel<=atol*100._wp
      report%c9=report%c9_max_rel<=atol*100._wp
      report%c10=report%min_cs2>=-cstol.and.report%max_cs2<1._wp+cstol.and. &
        report%min_cs2_fl>=-cstol.and.report%max_cs2_fl<1._wp+cstol
      report%c11=report%min_dmubdnb>=-1.e-10_wp
    endif
    report%overall=finite_primary.and.report%c1.and.report%c2.and.report%c3.and.report%c4.and. &
      report%c5.and.report%c6.and.report%c10.and.report%c11
    report%failed_checks=''
    if(.not.finite_primary)call add_failure(report,'FINITE')
    if(.not.report%c1)call add_failure(report,'C1')
    if(.not.report%c2)call add_failure(report,'C2')
    if(.not.report%c3)call add_failure(report,'C3')
    if(.not.report%c4)call add_failure(report,'C4')
    if(.not.report%c5)call add_failure(report,'C5')
    if(.not.report%c6)call add_failure(report,'C6')
    if(.not.report%c10)call add_failure(report,'C10')
    if(.not.report%c11)call add_failure(report,'C11')
  end subroutine check_consistency

  subroutine add_failure(report,name)
    type(verify_report_t),intent(inout)::report
    character(len=*),intent(in)::name
    integer::k
    k=len_trim(report%failed_checks)
    if(k==0)then
      report%failed_checks=trim(name)
    else
      report%failed_checks=trim(report%failed_checks)//','//trim(name)
    endif
  end subroutine add_failure

  pure function logical_text(x) result(s)
    logical,intent(in)::x
    character(len=5)::s
    if(x)then;s='true ';else;s='false';endif
  end function logical_text

  subroutine write_verify_json(path,eos,report)
    character(len=*),intent(in)::path
    type(eos_data_t),intent(in)::eos
    type(verify_report_t),intent(in)::report
    integer::u,ios
    open(newunit=u,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0)return
    write(u,'(a)')'{'
    write(u,'(a)')'"schema":"compose-eos-fortran-verification-v1",'
    write(u,'(a,i0,a)')'"n":',eos%n,','
    write(u,'(a,es25.16e3,a)')'"temperature_mev":',eos%temperature_mev,','
    write(u,'(a)')'"C1":{"pass":'//trim(logical_text(report%c1))//'},'
    write(u,'(a)')'"C2":{"applicable":'//trim(logical_text(report%c2_applicable))//',"pass":'//trim(logical_text(report%c2))//'},'
    if(report%c3_available)then
      write(u,'(a)')'"C3":{"status":"CHECKED","pass":'//trim(logical_text(report%c3))//'},'
      write(u,'(a)')'"C4":{"status":"CHECKED","pass":'//trim(logical_text(report%c4))//'},'
    else
      write(u,'(a)')'"C3":{"status":"NOT_AVAILABLE","pass":null},'
      write(u,'(a)')'"C4":{"status":"NOT_AVAILABLE","pass":null},'
    endif
    if(report%c5_applicable)then
      write(u,'(a)')'"C5":{"applicable":true,"pass":'//trim(logical_text(report%c5))//'},'
    else
      write(u,'(a)')'"C5":{"applicable":false,"pass":null},'
    endif
    write(u,'(a)')'"C6":{"pass":'//trim(logical_text(report%c6))//'},'
    write(u,'(a,es25.16e3,a,i0,a)')'"C7":{"status":"FINITE_DIFFERENCE_DIAGNOSTIC","max_rel":',report%c7_max_rel,',"n_points":',report%derivative_points,'},'
    write(u,'(a,es25.16e3,a,i0,a)')'"C8":{"status":"FINITE_DIFFERENCE_DIAGNOSTIC","max_rel":',report%c8_max_rel,',"n_points":',report%derivative_points,'},'
    write(u,'(a,es25.16e3,a,i0,a)')'"C9":{"status":"FINITE_DIFFERENCE_DIAGNOSTIC","max_rel":',report%c9_max_rel,',"n_points":',report%derivative_points,'},'
    write(u,'(a)')'"C10":{"pass":'//trim(logical_text(report%c10))//'},'
    write(u,'(a)')'"C11":{"pass":'//trim(logical_text(report%c11))//'},'
    write(u,'(a,i0,a)')'"derivative_points":',report%derivative_points,','
    write(u,'(a,i0,a)')'"derivative_masked_points":',report%derivative_masked,','
    write(u,'(a,es25.16e3,a)')'"C1_max_Q2":',report%c1_max_q2,','
    write(u,'(a,es25.16e3,a)')'"C1_max_Q6_minus_Q7":',report%c1_max_q6q7,','
    write(u,'(a,es25.16e3,a)')'"C1_max_mu_l":',report%c1_max_mul,','
    write(u,'(a,es25.16e3,a)')'"C2_max_relative_residual":',report%c2_max_rel,','
    write(u,'(a,es25.16e3,a)')'"C3_max_absolute_residual":',report%c3_max_abs,','
    write(u,'(a,es25.16e3,a)')'"C4_max_absolute_residual":',report%c4_max_abs,','
    write(u,'(a,es25.16e3,a)')'"C7_max_relative_residual":',report%c7_max_rel,','
    write(u,'(a,es25.16e3,a)')'"C8_max_relative_residual":',report%c8_max_rel,','
    write(u,'(a,es25.16e3,a)')'"C9_max_relative_residual":',report%c9_max_rel,','
    write(u,'(a,es25.16e3,a)')'"min_cs2":',report%min_cs2,','
    write(u,'(a,es25.16e3,a)')'"max_cs2":',report%max_cs2,','
    write(u,'(a,es25.16e3,a)')'"min_cs2_first_law":',report%min_cs2_fl,','
    write(u,'(a,es25.16e3,a)')'"max_cs2_first_law":',report%max_cs2_fl,','
    write(u,'(a,es25.16e3,a)')'"min_dmu_b_dnb":',report%min_dmubdnb,','
    write(u,'(a,i0,a)')'"acausal_point_count":',report%acausal_count,','
    write(u,'(a,i0,a)')'"acausal_point_count_first_law":',report%acausal_firstlaw_count,','
    write(u,'(a,es25.16e3,a)')'"min_dP_dnb":',report%min_dpdnb,','
    write(u,'(a,es25.16e3,a)')'"min_delta_P":',report%min_delta_p,','
    write(u,'(a)')'"failed_checks":"'//trim(report%failed_checks)//'",'
    write(u,'(a)')'"overall":{"status":"'//merge('PASS','FAIL',report%overall)//'","pass":'//trim(logical_text(report%overall))//'}'
    write(u,'(a)')'}'
    close(u)
  end subroutine write_verify_json
end module compose_eos_verify
