! Original source-locked cold EOS sound-speed evaluator.
! All thermodynamic evaluation is native Fortran (h=log(mu/mu_surface)).
! Source law sectors: surface polytrope, pressure-primary power laws,
! causal exponential n/mu bridge and upper cubic Hermite P(mu) intervals.
! ONE Maxwell interface: one-sided choice at the shared chemical potential.
! No interpolation across upper EOS sound-speed derivative knots.
! The geometry is separate, phase-local PCHIP; the c_s^2 field is NEVER
! taken from interpolated geometry tables. Runtime tables are process-local.
module compose_eos_cowling_source_exact
  use iso_c_binding, only: c_double,c_int,c_funloc,c_null_funptr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_dvode, only: cowling_native_set_source_callback, &
      cowling_native_source_callback_is_active
  implicit none
  private
  integer,save :: n_outer=0,n_upper=0,jump_row=0,n_geom=0,bin_hint=1,phase_side=0
  real(c_double),save :: mu_rest=0._c_double,mu_surface_anchor=0._c_double
  real(c_double),save :: mu_bridge_left=0._c_double,mu_bridge_right=0._c_double
  real(c_double),save :: bridge_width=0._c_double,bridge_yl=0._c_double
  real(c_double),save :: bridge_dy=0._c_double,bridge_kappa=0._c_double
  real(c_double),save :: surface_gamma=0._c_double
  real(c_double),save :: hstart_phase=0._c_double
  real(c_double),allocatable,save :: outer(:,:),outer_mu(:),upper(:,:),upper_mu(:)
  real(c_double),allocatable,save :: geom_s(:),geom_v(:,:),geom_d(:,:)
  public :: cowling_native_source_eos_setup,cowling_native_source_cs2_eval
  public :: cowling_native_source_phase_setup,cowling_native_source_register
  public :: cowling_native_source_unregister
  public :: cowling_native_source_select_restarts
contains
  pure real(c_double) function expm1_accurate(z) result(y)
    real(c_double),intent(in) :: z
    if(abs(z)<1.e-5_c_double)then
      y=z*(1._c_double+z*(.5_c_double+z*(1._c_double/6._c_double+ &
        z*(1._c_double/24._c_double+z*(1._c_double/120._c_double+ &
        z/720._c_double)))))
    else
      y=exp(z)-1._c_double
    endif
  end function
  ! Inputs are physical contract rows, *not* reconstructed PCHIP data.
  ! surface = [mu_surface(=m),n0,p0,gamma].
  ! outer (4,N)=[n,p,F=epsilon/n,gamma_interval] (last gamma not used).
  ! bridge=[n_L,P_L,epsilon_L,n_R,P_R,epsilon_R,kappa].
  ! upper (3,N)=[n,P,epsilon] including the two Maxwell plateau rows.
  ! index_jump is the original zero-based low-density upper-table row.
  subroutine cowling_native_source_eos_setup(no,outer_in,nu,upper_in, &
      surface,bridge,index_jump,status) bind(c)
    integer(c_int),value :: no,nu,index_jump
    real(c_double),intent(in)::outer_in(4,no),upper_in(3,nu),surface(4),bridge(7)
    integer(c_int),intent(out)::status
    real(c_double)::mu_l,mu_r,scale
    integer::i,err
    status=1
    n_outer=0;n_upper=0
    if(no<2.or.nu<4.or.index_jump<1.or.index_jump>=nu-2)return
    if(.not.all(ieee_is_finite(surface)).or..not.all(ieee_is_finite(bridge)))return
    if(.not.all(ieee_is_finite(outer_in)).or..not.all(ieee_is_finite(upper_in)))return
    if(any(surface<=0._c_double).or.surface(4)<=1._c_double)return
    if(any(outer_in(1:3,:)<=0._c_double).or. &
       any(outer_in(4,:no-1)<=0._c_double))return
    if(any(upper_in<=0._c_double).or.any(bridge(1:6)<=0._c_double))return
    if(any(outer_in(1,2:no)<=outer_in(1,1:no-1)).or. &
       any(outer_in(2,2:no)<=outer_in(2,1:no-1)))return
    if(any(upper_in(1,2:nu)<=upper_in(1,1:nu-1)))return
    ! Exactly one physically declared Maxwell doublet can differ by a
    ! final-bit pressure encoding. All other source intervals must strictly
    ! increase in P, and the declared doublet is checked independently.
    do i=1,nu-1
      if(i==index_jump+1)cycle
      if(upper_in(2,i+1)<=upper_in(2,i))return
    enddo
    if(allocated(outer))deallocate(outer)
    if(allocated(outer_mu))deallocate(outer_mu)
    if(allocated(upper))deallocate(upper)
    if(allocated(upper_mu))deallocate(upper_mu)
    allocate(outer(4,no),outer_mu(no),upper(3,nu),upper_mu(nu),stat=err)
    if(err/=0)return
    outer=outer_in;upper=upper_in
    outer_mu=outer(3,:)+outer(2,:)/outer(1,:)
    upper_mu=(upper(3,:)+upper(2,:))/upper(1,:)
    if(any(outer_mu(2:no)<=outer_mu(1:no-1)))return
    jump_row=index_jump+1
    if(any(upper_mu(2:jump_row)<=upper_mu(1:jump_row-1)))return
    if(any(upper_mu(jump_row+2:nu)<=upper_mu(jump_row+1:nu-1)))return
    scale=max(1._c_double,abs(upper_mu(jump_row)),abs(upper_mu(jump_row+1)))
    if(abs(upper_mu(jump_row)-upper_mu(jump_row+1))>2.e-12_c_double*scale)return
    scale=max(1._c_double,abs(upper(2,jump_row)),abs(upper(2,jump_row+1)))
    if(abs(upper(2,jump_row)-upper(2,jump_row+1))>2.e-12_c_double*scale)return
    mu_rest=surface(1)
    surface_gamma=surface(4)
    mu_surface_anchor=mu_rest+surface_gamma/(surface_gamma-1._c_double)* &
       surface(3)/surface(2)
    mu_l=(bridge(3)+bridge(2))/bridge(1)
    mu_r=(bridge(6)+bridge(5))/bridge(4)
    bridge_width=mu_r-mu_l
    if(bridge_width<=0._c_double.or.bridge(7)<-200._c_double.or. &
       bridge(7)>200._c_double)return
    mu_bridge_left=mu_l;mu_bridge_right=mu_r
    bridge_yl=bridge(1)/mu_l
    bridge_dy=bridge(4)/mu_r-bridge_yl
    bridge_kappa=bridge(7)
    if(bridge_dy<0._c_double)return
    if(abs(outer_mu(1)-mu_surface_anchor)>1.e-8_c_double*mu_surface_anchor)return
    if(abs(outer_mu(no)-mu_l)>1.e-8_c_double*mu_l)return
    if(abs(mu_r-upper_mu(1))>1.e-8_c_double*mu_r)return
    if(mu_rest<=0._c_double.or.mu_surface_anchor<=mu_rest) return
    if(any(upper_mu(1:jump_row)<mu_bridge_right))return
    n_outer=no;n_upper=nu
    status=0
  end subroutine

  ! Binary interval selection equivalent to source searchsorted(...,right).
  integer function bracket_right(x,array) result(i)
    real(c_double),intent(in)::x,array(:)
    integer::lo,hi,mid,n
    n=size(array)
    lo=1;hi=n
    if(x<=array(1))then
      i=1;return
    endif
    if(x>=array(n))then
      i=n-1;return
    endif
    do while(hi-lo>1)
      mid=(lo+hi)/2
      if(array(mid)<=x)then
        lo=mid
      else
        hi=mid
      endif
    enddo
    i=lo
  end function

  subroutine cowling_native_source_cs2_eval(h,side,cs2,sector,status) bind(c)
    real(c_double),value :: h
    integer(c_int),value :: side
    real(c_double),intent(out)::cs2
    integer(c_int),intent(out)::sector,status
    real(c_double)::mu,g,q0,tl,tr,ht,band,ratio
    real(c_double)::t,phi,dphi,y,a,b,dp,d0,d1,a2,a3,n,dndmu,span
    integer::i
    cs2=0._c_double;sector=0;status=1
    if(n_outer<2.or.n_upper<4.or.(side/=1.and.side/=2))return
    if(.not.ieee_is_finite(h).or.h<0._c_double)return
    if(h>log(upper_mu(n_upper)/mu_rest)+2.e-13_c_double)return
    mu=mu_rest*exp(h)
    if(.not.ieee_is_finite(mu))return
    if(mu<=mu_surface_anchor)then
      sector=1
      cs2=(surface_gamma-1._c_double)*(mu-mu_rest)/mu
    else if(mu<=outer_mu(n_outer))then
      sector=2
      i=bracket_right(mu,outer_mu)
      g=outer(4,i)
      q0=outer(2,i)/outer(1,i)
      cs2=((g-1._c_double)*(mu-outer(3,i))+q0)/mu
    else if(mu<=mu_bridge_right)then
      sector=3
      t=max(0._c_double,min(1._c_double,(mu-mu_bridge_left)/bridge_width))
      if(bridge_kappa==0._c_double)then
        phi=t;dphi=1._c_double
      else
        phi=expm1_accurate(bridge_kappa*t)/expm1_accurate(bridge_kappa)
        dphi=bridge_kappa*exp(bridge_kappa*t)/expm1_accurate(bridge_kappa)
      endif
      y=bridge_yl+bridge_dy*phi
      cs2=y/(y+mu*bridge_dy*dphi/bridge_width)
    else
      sector=4
      tl=log(upper_mu(jump_row)/mu_rest)
      tr=log(upper_mu(jump_row+1)/mu_rest)
      band=8._c_double*epsilon(1._c_double)*max(1._c_double,abs(tl))
      if(h>=min(tl,tr)-band.and.h<=max(tl,tr)+band)then
        if(side==1)then
          i=jump_row-1
        else
          i=jump_row+1
        endif
      else
        i=bracket_right(mu,upper_mu)
      endif
      if(i==jump_row.or.i<1.or.i>=n_upper)return
      span=upper_mu(i+1)-upper_mu(i)
      if(span<=0._c_double)return
      t=max(0._c_double,min(1._c_double,(mu-upper_mu(i))/span))
      dp=upper(2,i+1)-upper(2,i)
      d0=span*upper(1,i)
      d1=span*upper(1,i+1)
      a2=3._c_double*dp-2._c_double*d0-d1
      a3=-2._c_double*dp+d0+d1
      n=(d0+t*(2._c_double*a2+3._c_double*a3*t))/span
      dndmu=(2._c_double*a2+6._c_double*a3*t)/(span*span)
      if(n<=0._c_double.or.dndmu<=0._c_double)return
      cs2=n/(mu*dndmu)
    endif
    if(.not.ieee_is_finite(cs2))then
      cs2=0._c_double
      sector=0
      return
    endif
    if(cs2<0._c_double.or.cs2>1._c_double+1.e-10_c_double)then
      cs2=0._c_double
      sector=0
      return
    endif
    status=0
  end subroutine

  ! Physical one-sided derivative jumps which require fresh DVODE history.
  ! This reproduces the previous source-contract criterion (not a PCHIP
  ! numerical cutoff): compare equilibrium c_s² at h +/- dh strictly
  ! inside one physical phase. Outputs s are sorted ascending.
  subroutine cowling_native_source_select_restarts(capacity,hstart,hbegin,hend, &
           relative_jump,side,cuts,ncuts,status) bind(c)
    integer(c_int),value :: capacity,side
    real(c_double),value :: hstart,hbegin,hend,relative_jump
    real(c_double),intent(out)::cuts(capacity)
    integer(c_int),intent(out)::ncuts,status
    real(c_double)::h,dh,csl,csh,scale,sx,mu
    integer(c_int)::rlo,rhi,istatlo,isthi
    integer::index,j,k,kind
    ncuts=0;status=1;cuts=0._c_double
    if(n_outer<2.or.n_upper<4.or.capacity<2)return
    if(side/=1.and.side/=2)return
    if(.not.all(ieee_is_finite([hstart,hbegin,hend,relative_jump])))return
    if(.not.(hstart>=hbegin.and.hbegin>hend.and.hend>=0._c_double))return
    if(.not.(relative_jump>0._c_double.and.relative_jump<1._c_double))return
    do index=1,n_outer+n_upper+3
      if(index==1)then
        mu=mu_rest
      else if(index==2)then
        mu=mu_surface_anchor
      else if(index<=n_outer+2)then
        mu=outer_mu(index-2)
      else if(index==n_outer+3)then
        mu=mu_bridge_right
      else
        mu=upper_mu(index-(n_outer+3))
      endif
      if(mu<=0._c_double)cycle
      h=log(mu/mu_rest)
      if(h<=hend.or.h>=hbegin)cycle
      dh=max(2.e-11_c_double,1.e-9_c_double*max(h,1.e-5_c_double))
      if(h-dh<=hend.or.h+dh>=hbegin)cycle
      call cowling_native_source_cs2_eval(h-dh,side,csl,rlo,istatlo)
      call cowling_native_source_cs2_eval(h+dh,side,csh,rhi,isthi)
      if(istatlo/=0.or.isthi/=0)then
        status=4
        ncuts=0;cuts=0._c_double
        return
      endif
      if(min(csl,csh)<=0._c_double)then
        status=5
        ncuts=0;cuts=0._c_double
        return
      endif
      if(abs(csh-csl)/max(abs(csh),abs(csl))<relative_jump)cycle
      sx=sqrt(max(0._c_double,hstart-h))
      ! Floating-point duplicates at a phase boundary are not extra
      ! physical transitions. Source knots never merge across the jump.
      if(ncuts>0)then
        if(any(abs(cuts(1:ncuts)-sx)<5.e-13_c_double))cycle
      endif
      if(ncuts>=capacity)then
        status=3
        ncuts=0;cuts=0._c_double
        return
      endif
      ncuts=ncuts+1
      cuts(ncuts)=sx
    enddo
    ! Low/source mu order produces decreasing s: keep caller-independent
    ! increasing s; reverse direction only when building ODE targets.
    do j=2,ncuts
      sx=cuts(j)
      k=j-1
      do while(k>=1)
        if(cuts(k)<=sx)exit
        cuts(k+1)=cuts(k)
        k=k-1
      enddo
      cuts(k+1)=sx
    enddo
    status=0
  end subroutine

  ! Geometry PCHIP stays separate from exact source c_s². No source law
  ! is interpolated across an EOS knot or a Maxwell interface.
  subroutine cowling_native_source_phase_setup(n,s,v,d,hstart,side,status) bind(c)
    integer(c_int),value :: n,side
    real(c_double),intent(in) :: s(n),v(5,n),d(5,n)
    real(c_double),value :: hstart
    integer(c_int),intent(out) :: status
    integer::err
    status=1;n_geom=0
    if(n<4.or.n_outer<2.or.side<1.or.side>2)return
    if(.not.ieee_is_finite(hstart).or.hstart<=0._c_double)return
    if(.not.all(ieee_is_finite(s)).or..not.all(ieee_is_finite(v)) &
       .or..not.all(ieee_is_finite(d)))return
    if(any(s(2:n)<=s(1:n-1)).or.s(1)<0._c_double)return
    if(allocated(geom_s))deallocate(geom_s)
    if(allocated(geom_v))deallocate(geom_v)
    if(allocated(geom_d))deallocate(geom_d)
    allocate(geom_s(n),geom_v(5,n),geom_d(5,n),stat=err)
    if(err/=0)return
    geom_s=s;geom_v=v;geom_d=d
    hstart_phase=hstart;phase_side=side;n_geom=n;bin_hint=1
    status=0
  end subroutine

  ! Standard Hermite table with cached last interval for nearby DVODE
  ! evaluations. No Python callback nor Python numerical EOS evaluation.
  subroutine native_source_full_callback(s,v,status) bind(c)
    real(c_double),value :: s
    real(c_double),intent(out)::v(5)
    integer(c_int),intent(out)::status
    real(c_double)::dx,z,z2,z3,cs2
    integer(c_int)::region,physics_status
    integer::lo,hi,mid
    status=1;v=0._c_double
    if(n_geom<4.or..not.ieee_is_finite(s))return
    if(s<geom_s(1)-1.e-11_c_double.or.s>geom_s(n_geom)+1.e-11_c_double)return
    if(s<=geom_s(1))then
      v=geom_v(:,1)
    else if(s>=geom_s(n_geom))then
      v=geom_v(:,n_geom)
    else
      lo=max(1,min(bin_hint,n_geom-1))
      if(.not.(geom_s(lo)<=s.and.s<geom_s(lo+1)))then
        if(lo>1)then
          if(geom_s(lo-1)<=s.and.s<geom_s(lo))lo=lo-1
        endif
        if(lo<n_geom-1)then
          if(geom_s(lo+1)<=s.and.s<geom_s(lo+2))lo=lo+1
        endif
        if(.not.(geom_s(lo)<=s.and.s<geom_s(lo+1)))then
          lo=1;hi=n_geom
          do while(hi-lo>1)
            mid=(lo+hi)/2
            if(geom_s(mid)<=s)then
              lo=mid
            else
              hi=mid
            endif
          enddo
        endif
      endif
      bin_hint=lo
      dx=geom_s(lo+1)-geom_s(lo)
      z=(s-geom_s(lo))/dx;z2=z*z;z3=z2*z
      v=(2*z3-3*z2+1)*geom_v(:,lo)+(z3-2*z2+z)*dx*geom_d(:,lo) &
       +(-2*z3+3*z2)*geom_v(:,lo+1)+(z3-z2)*dx*geom_d(:,lo+1)
    endif
    call cowling_native_source_cs2_eval(hstart_phase-s*s, &
         int(phase_side,c_int),cs2,region,physics_status)
    if(physics_status/=0)return
    v(5)=cs2
    if(.not.all(ieee_is_finite(v)).or.v(1)<=0._c_double &
        .or.v(4)<=0._c_double.or.v(5)<=0._c_double)return
    status=0
  end subroutine

  subroutine cowling_native_source_register(status) bind(c)
    integer(c_int),intent(out)::status
    status=1
    if(n_outer<2.or.n_geom<4)return
    if(cowling_native_source_callback_is_active())then
      status=2;return
    endif
    call cowling_native_set_source_callback(c_funloc(native_source_full_callback),status)
  end subroutine

  subroutine cowling_native_source_unregister(status) bind(c)
    integer(c_int),intent(out)::status
    call cowling_native_set_source_callback(c_null_funptr,status)
  end subroutine
end module compose_eos_cowling_source_exact
