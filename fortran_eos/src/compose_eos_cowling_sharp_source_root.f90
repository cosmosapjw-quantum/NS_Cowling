! Source-exact, one-Maxwell-interface Cowling eigenfrequency polishing.
! Every bracket trial, EOS c_s^2 evaluation, strong-knot DVODE restart and
! slow/rapid junction transport executes inside the Fortran library.
! Geometry uses independent one-sided PCHIP, source c_s^2 uses ORIGINAL
! contract laws. Neither frozen-composition nor full-GR perturbations.
! The legacy DVODE COMMON is non-reentrant: use one process per worker.
module compose_eos_cowling_sharp_source_root
  use iso_c_binding, only: c_double,c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_dvode, only: cowling_native_setup, &
        cowling_native_segment, cowling_native_segment_restarts, &
        cowling_native_interface_transfer, &
        cowling_native_source_callback_is_active
  use compose_eos_cowling_source_exact, only: &
        cowling_native_source_select_restarts, &
        cowling_native_source_phase_setup, &
        cowling_native_source_register,cowling_native_source_unregister
  implicit none
  private
  real(c_double),parameter :: pi=3.1415926535897932384626433832795_c_double
  real(c_double),parameter :: c_km_s=299792.458_c_double
  public :: cowling_native_sharp_source_brent, exact_onejump_det
contains
  ! Integrate exactly one physical phase leg. Build and schedule all chosen
  ! one-sided derivative-discontinuity restarts in Fortran, preserving DVODE
  ! initial state continuity and creating a FRESH DVODE history per segment.
  subroutine exact_one_leg(omega,ell,rtol,atol,max_ds,mf, &
      n,s,v,d,hstart,side,cuts,start_s,finish_s,y0,y1,nfev,status)
    real(c_double),intent(in) :: omega,rtol,atol,max_ds,hstart,start_s,finish_s
    integer(c_int),intent(in) :: ell,mf,n,side
    real(c_double),intent(in) :: s(n),v(5,n),d(5,n),cuts(:),y0(2)
    real(c_double),intent(out)::y1(2)
    integer(c_int),intent(out)::nfev,status
    real(c_double),allocatable::targets(:),outputs(:,:)
    integer(c_int),allocatable::cut_indices(:)
    integer(c_int)::st,register_status,integrator_status,cleanup_status
    integer :: j,npts,ncuts,k
    nfev=0;y1=0._c_double;status=1
    if(.not.(start_s/=finish_s))return
    ncuts=count(cuts>min(start_s,finish_s).and. &
                cuts<max(start_s,finish_s))
    npts=ncuts+2
    allocate(targets(npts),outputs(2,npts),cut_indices(ncuts))
    targets(1)=start_s;targets(npts)=finish_s
    k=0
    if(start_s<finish_s)then
      do j=1,size(cuts)
        if(cuts(j)<=start_s.or.cuts(j)>=finish_s)cycle
        k=k+1;targets(k+1)=cuts(j)
      enddo
    else
      do j=size(cuts),1,-1
        if(cuts(j)<=finish_s.or.cuts(j)>=start_s)cycle
        k=k+1;targets(k+1)=cuts(j)
      enddo
    endif
    if(k/=ncuts)return
    do j=1,ncuts
      cut_indices(j)=int(j+1,c_int)
    enddo
    call cowling_native_setup(n,s,v,d,st)
    if(st/=0)then
      status=10+st;return
    endif
    call cowling_native_source_phase_setup(n,s,v,d,hstart,side,st)
    if(st/=0)then
      status=20+st;return
    endif
    call cowling_native_source_register(register_status)
    if(register_status/=0)then
      status=30+register_status;return
    endif
    if(ncuts>0)then
      call cowling_native_segment_restarts(omega,ell,rtol,atol,max_ds,mf, &
        int(npts,c_int),targets,int(ncuts,c_int),cut_indices,y0,outputs, &
        nfev,integrator_status)
    else
      call cowling_native_segment(omega,ell,rtol,atol,max_ds,mf, &
        int(npts,c_int),targets,y0,outputs,nfev,integrator_status)
    endif
    ! Release the Fortran-to-Fortran source pointer even on integration failure.
    call cowling_native_source_unregister(cleanup_status)
    if(cleanup_status/=0)then
      status=41;nfev=0;return
    endif
    if(integrator_status/=0)then
      status=50+integrator_status;nfev=0;return
    endif
    if(.not.all(ieee_is_finite(outputs(:,npts))))then
      status=69;nfev=0;return
    endif
    y1=outputs(:,npts)
    status=0
  end subroutine

  ! Returns the original signed matching determinant for the three legs:
  ! center -> match (high); surface -> interface (low);
  ! outer->inner Maxwell transfer; interface -> match (high).
  subroutine exact_onejump_det(f_hz,ell,rtol,atol,max_ds,mf, &
      nh,sh,vh,dh,nl,sl,vl,dl,hstart,smatch,geom,conversion, &
      cuts_h,cuts_l,residual,nfev,status)
    real(c_double),intent(in)::f_hz,rtol,atol,max_ds,hstart,smatch
    integer(c_int),intent(in)::ell,mf,nh,nl,conversion
    real(c_double),intent(in)::sh(nh),vh(5,nh),dh(5,nh)
    real(c_double),intent(in)::sl(nl),vl(5,nl),dl(5,nl),geom(7)
    real(c_double),intent(in)::cuts_h(:),cuts_l(:)
    real(c_double),intent(out)::residual
    integer(c_int),intent(out)::nfev,status
    real(c_double)::omega,asurf,bsurf,ui(2),ulo(2),uo(2),after(2)
    real(c_double)::y0(2),traction,kinematic,den
    integer(c_int)::n1,n2,n3,s1,s2,s3,transfer_status
    residual=0._c_double;nfev=0;status=1
    omega=2._c_double*pi*f_hz/c_km_s
    call exact_one_leg(omega,ell,rtol,atol,max_ds,mf,nh,sh,vh,dh, &
       hstart,2_c_int,cuts_h,sh(1),smatch, &
       [1._c_double,-1._c_double/real(ell,c_double)],ui,n1,s1)
    if(s1/=0)then
      status=10+s1;return
    endif
    if(vl(3,nl)-2._c_double*vl(2,nl)>600._c_double)then
      status=21;return
    endif
    asurf=omega*omega*exp(vl(3,nl)-2._c_double*vl(2,nl))*vl(1,nl)
    bsurf=vl(4,nl)
    if(.not.ieee_is_finite(asurf).or.asurf<=tiny(asurf) &
       .or..not.ieee_is_finite(bsurf).or.bsurf<=0._c_double)then
      status=22;return
    endif
    y0=[1._c_double,-bsurf/asurf]
    if(.not.all(ieee_is_finite(y0)))then
      status=23;return
    endif
    call exact_one_leg(omega,ell,rtol,atol,max_ds,mf,nl,sl,vl,dl, &
       hstart,1_c_int,cuts_l,sl(nl),sl(1),y0,ulo,n2,s2)
    if(s2/=0)then
      status=30+s2;return
    endif
    ! Both sides have common geometry and pressure, different epsilon.
    call cowling_native_interface_transfer(omega,geom(1),geom(2),geom(3), &
       geom(4),geom(5),geom(6),geom(7),conversion,ulo,after, &
       traction,kinematic,transfer_status)
    if(transfer_status/=0)then
      status=41;return
    endif
    if(max(traction,kinematic)>1.e-9_c_double)then
      status=42;return
    endif
    call exact_one_leg(omega,ell,rtol,atol,max_ds,mf,nh,sh,vh,dh, &
       hstart,2_c_int,cuts_h,sh(nh),smatch,after,uo,n3,s3)
    if(s3/=0)then
      status=50+s3;return
    endif
    den=sqrt(sum(ui*ui))*sqrt(sum(uo*uo))
    if(.not.ieee_is_finite(den).or.den<=tiny(den))then
      status=61;return
    endif
    residual=(ui(1)*uo(2)-ui(2)*uo(1))/den
    if(.not.ieee_is_finite(residual))then
      residual=0._c_double;status=62;return
    endif
    nfev=n1+n2+n3
    status=0
  end subroutine

  ! Native exact-contract strong-knot scheduling and bracketed Brent-Dekker.
  ! source EOS contract MUST have been supplied through the existing
  ! cowling_native_source_eos_setup API in this serialized process first.
  ! Both strong-knot lists and every root iteration stay in Fortran.
  subroutine cowling_native_sharp_source_brent( &
      f_lo,f_hi,ell,rtol,atol,max_ds,mf, &
      nh,sh,vh,dh,nl,sl,vl,dl,hstart,smatch,geom,conversion, &
      xtol_hz,relative_jump,root_hz,root_residual,iterations,total_nfev, &
      restart_count,status) bind(c)
    real(c_double),value::f_lo,f_hi,rtol,atol,max_ds,hstart,smatch,xtol_hz,relative_jump
    integer(c_int),value::ell,mf,nh,nl,conversion
    real(c_double),intent(in)::sh(nh),vh(5,nh),dh(5,nh)
    real(c_double),intent(in)::sl(nl),vl(5,nl),dl(5,nl),geom(7)
    real(c_double),intent(out)::root_hz,root_residual
    integer(c_int),intent(out)::iterations,total_nfev,restart_count,status
    real(c_double),allocatable::candidate_h(:),candidate_l(:)
    real(c_double)::a,b,c,fa,fb,fc,d,e,s,p,q,r,tol,xmid
    integer(c_int)::nc_h,nc_l,sel_h,sel_l,nfev,shot_status
    integer::iteration
    root_hz=0._c_double;root_residual=0._c_double
    iterations=0;total_nfev=0;restart_count=0;status=1
    if(cowling_native_source_callback_is_active())then
      status=2;return
    endif
    if(nh<4.or.nl<4.or.ell<2)return
    if(mf/=10.and.mf/=22)return
    if(conversion/=1.and.conversion/=2)return
    if(.not.all(ieee_is_finite([f_lo,f_hi,rtol,atol,max_ds, &
                hstart,smatch,xtol_hz,relative_jump])))return
    if(.not.all(ieee_is_finite(sh)).or..not.all(ieee_is_finite(sl)))return
    if(.not.all(ieee_is_finite(vh)).or..not.all(ieee_is_finite(vl)))return
    if(.not.all(ieee_is_finite(dh)).or..not.all(ieee_is_finite(dl)))return
    if(.not.all(ieee_is_finite(geom)))return
    if(f_lo<=0._c_double.or.f_hi<=f_lo)return
    if(xtol_hz<=0._c_double.or.xtol_hz>=f_hi-f_lo)return
    if(rtol<1.e-13_c_double.or.rtol>1.e-3_c_double)return
    if(atol<=0._c_double.or.max_ds<=0._c_double)return
    if(relative_jump<=0._c_double.or.relative_jump>=1._c_double)return
    if(hstart<=0._c_double.or.abs(sh(1))>1.e-12_c_double)return
    if(smatch<=sh(1).or.smatch>=sh(nh))return
    if(abs(sl(1)-sh(nh))>1.e-12_c_double)return
    if(sl(nl)<=sl(1))return
    if(hstart-sl(nl)**2<0._c_double)return
    if(any(sh(2:nh)<=sh(1:nh-1)))return
    if(any(sl(2:nl)<=sl(1:nl-1)))return
    if(any(abs(vh(1:4,nh)-vl(1:4,1))> &
      1.e-8_c_double*max(1._c_double,abs(geom(1:4)))))return
    if(any(abs(vh(1:4,nh)-geom(1:4))> &
      1.e-8_c_double*max(1._c_double,abs(geom(1:4)))))return
    if(geom(1)<=0._c_double.or.geom(4)<=0._c_double.or.geom(5)<0._c_double)return
    if(geom(6)<=0._c_double.or.geom(7)<=geom(6))return
    allocate(candidate_h(nh),candidate_l(nl))
    call cowling_native_source_select_restarts(nh,hstart, &
       hstart-sh(1)**2,hstart-sh(nh)**2,relative_jump,2_c_int, &
       candidate_h,nc_h,sel_h)
    if(sel_h/=0)then
      status=71+sel_h
      return
    endif
    call cowling_native_source_select_restarts(nl,hstart, &
       hstart-sl(1)**2,hstart-sl(nl)**2,relative_jump,1_c_int, &
       candidate_l,nc_l,sel_l)
    if(sel_l/=0)then
      status=81+sel_l
      return
    endif
    restart_count=nc_h+nc_l
    a=f_lo;b=f_hi
    call exact_onejump_det(a,ell,rtol,atol,max_ds,mf, &
       nh,sh,vh,dh,nl,sl,vl,dl,hstart,smatch,geom,conversion, &
       candidate_h(1:nc_h),candidate_l(1:nc_l),fa,nfev,shot_status)
    if(shot_status/=0.or..not.ieee_is_finite(fa))then
      status=100+shot_status;return
    endif
    total_nfev=total_nfev+nfev
    call exact_onejump_det(b,ell,rtol,atol,max_ds,mf, &
       nh,sh,vh,dh,nl,sl,vl,dl,hstart,smatch,geom,conversion, &
       candidate_h(1:nc_h),candidate_l(1:nc_l),fb,nfev,shot_status)
    if(shot_status/=0.or..not.ieee_is_finite(fb))then
      status=100+shot_status;total_nfev=0;return
    endif
    total_nfev=total_nfev+nfev
    if(fa==0._c_double.or.fb==0._c_double)then
      if(fa==0._c_double)then
        root_hz=a;root_residual=fa
      else
        root_hz=b;root_residual=fb
      endif
      status=0
      return
    endif
    if(sign(1._c_double,fa)==sign(1._c_double,fb))then
      status=3;total_nfev=0;return
    endif
    c=a;fc=fa;d=b-a;e=d
    do iteration=1,128
      if(abs(fc)<abs(fb))then
        a=b;b=c;c=a
        fa=fb;fb=fc;fc=fa
      endif
      tol=2._c_double*epsilon(1._c_double)*abs(b)+0.5_c_double*xtol_hz
      xmid=0.5_c_double*(c-b)
      if(abs(xmid)<=tol.or.fb==0._c_double)then
        if(abs(fb)>1.e-5_c_double)then
          total_nfev=0;status=33;return
        endif
        root_hz=b;root_residual=fb;iterations=iteration-1
        status=0;return
      endif
      if(abs(e)>=tol.and.abs(fa)>abs(fb))then
        s=fb/fa
        if(a==c)then
          p=2._c_double*xmid*s
          q=1._c_double-s
        else
          q=fa/fc
          r=fb/fc
          p=s*(2._c_double*xmid*q*(q-r)-(b-a)*(r-1._c_double))
          q=(q-1._c_double)*(r-1._c_double)*(s-1._c_double)
        endif
        if(p>0._c_double)q=-q
        p=abs(p)
        if(q/=0._c_double.and. &
             2._c_double*p<3._c_double*xmid*q-abs(tol*q).and. &
             p<abs(.5_c_double*e*q))then
          e=d;d=p/q
        else
          d=xmid;e=d
        endif
      else
        d=xmid;e=d
      endif
      a=b;fa=fb
      if(abs(d)>tol)then
        b=b+d
      else
        b=b+sign(tol,xmid)
      endif
      call exact_onejump_det(b,ell,rtol,atol,max_ds,mf, &
        nh,sh,vh,dh,nl,sl,vl,dl,hstart,smatch,geom,conversion, &
        candidate_h(1:nc_h),candidate_l(1:nc_l),fb,nfev,shot_status)
      if(shot_status/=0.or..not.ieee_is_finite(fb))then
        total_nfev=0;status=100+shot_status;return
      endif
      total_nfev=total_nfev+nfev
      if(sign(1._c_double,fb)==sign(1._c_double,fc))then
        c=a;fc=fa;d=b-a;e=d
      endif
    enddo
    total_nfev=0;status=32
  end subroutine
end module compose_eos_cowling_sharp_source_root
