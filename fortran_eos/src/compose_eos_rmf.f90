module compose_eos_rmf
  use compose_eos_types, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  real(wp), parameter :: hbarc = 197.3269804_wp
  real(wp), parameter :: pi = acos(-1._wp)
  real(wp), parameter :: m_e = 0.51099895_wp
  real(wp), parameter :: m_mu = 105.6583755_wp

  type, public :: ddrmf_model_t
    character(len=16) :: id = ''
    real(wp) :: mn = 0._wp, mp = 0._wp
    real(wp) :: rho_sat = 0._wp
    real(wp) :: ms = 0._wp, mw = 0._wp, mr = 0._wp
    real(wp) :: gs0 = 0._wp, gw0 = 0._wp, gr0 = 0._wp
    real(wp) :: as = 0._wp, bs = 0._wp, cs = 0._wp, ds = 0._wp
    real(wp) :: aw = 0._wp, bw = 0._wp, cw = 0._wp, dw = 0._wp
    real(wp) :: ar = 0._wp
  end type ddrmf_model_t

  type, public :: ddrmf_state_t
    real(wp) :: rho_b = 0._wp, yp = 0._wp
    real(wp) :: rho_n = 0._wp, rho_p = 0._wp, rho_3 = 0._wp
    real(wp) :: mstar_n = 0._wp, mstar_p = 0._wp
    real(wp) :: sigma_scalar = 0._wp, rho_scalar = 0._wp
    real(wp) :: sigma_r = 0._wp
    real(wp) :: mu_n = 0._wp, mu_p = 0._wp
    real(wp) :: eps_b = 0._wp, p_b = 0._wp, p_thermo = 0._wp
    real(wp) :: pressure_identity_residual = 0._wp
    real(wp) :: scalar_field_residual = 0._wp
  end type ddrmf_state_t

  type, public :: ddrmf_beta_result_t
    type(ddrmf_state_t) :: baryons
    real(wp) :: yp = 0._wp, ye = 0._wp, ymu = 0._wp
    real(wp) :: rho_e = 0._wp, rho_mu = 0._wp
    real(wp) :: mu_e = 0._wp, mu_mu = -1._wp
    real(wp) :: eps_leptons = 0._wp, p_leptons = 0._wp
    real(wp) :: eps_total = 0._wp, p_total = 0._wp
    real(wp) :: beta_residual = 0._wp, charge_residual = 0._wp
    real(wp) :: muon_onset_rho_p = 0._wp, muon_onset_gap = 0._wp
    logical :: muon_active = .false.
  end type ddrmf_beta_result_t

  type, public :: ddrmf_diagnostic_t
    type(ddrmf_beta_result_t) :: beta
    real(wp) :: symmetry_second_derivative = 0._wp
    real(wp) :: symmetry_pnm_minus_snm = 0._wp
    real(wp) :: symmetry_step = 0.05_wp
  end type ddrmf_diagnostic_t

  public :: load_ddrmf_model, solve_ddrmf_homogeneous
  public :: solve_ddrmf_beta_equilibrium, evaluate_ddrmf

contains

  subroutine load_ddrmf_model(model_id, model, status)
    character(len=*), intent(in) :: model_id
    type(ddrmf_model_t), intent(out) :: model
    character(len=*), intent(out) :: status

    model = ddrmf_model_t()
    status = 'unknown DD-RMF model: '//trim(model_id)
    select case(trim(model_id))
    case('GPPVA_DD2')
      model%id='GPPVA_DD2'
      model%mn=939.56536_wp; model%mp=938.27203_wp
      model%rho_sat=0.149065_wp
      model%ms=546.212459_wp; model%mw=783._wp; model%mr=763._wp
      model%gs0=10.686681_wp; model%gw0=13.342362_wp; model%gr0=3.626940_wp
      model%as=1.357630_wp; model%bs=0.634442_wp; model%cs=1.005358_wp; model%ds=0.575810_wp
      model%aw=1.369718_wp; model%bw=0.496475_wp; model%cw=0.817753_wp; model%dw=0.638452_wp
      model%ar=0.518903_wp
    case('GPPVA_DDME2')
      model%id='GPPVA_DDME2'
      model%mn=939._wp; model%mp=939._wp
      model%rho_sat=0.152_wp
      model%ms=550.1238_wp; model%mw=783._wp; model%mr=763._wp
      model%gs0=10.5396_wp; model%gw0=13.0189_wp; model%gr0=3.6836_wp
      model%as=1.3881_wp; model%bs=1.0943_wp; model%cs=1.7057_wp; model%ds=0.4421_wp
      model%aw=1.3892_wp; model%bw=0.9240_wp; model%cw=1.4620_wp; model%dw=0.4775_wp
      model%ar=0.5647_wp
    case('GPPVA_TW')
      model%id='GPPVA_TW'
      model%mn=939._wp; model%mp=939._wp
      model%rho_sat=0.153_wp
      model%ms=550._wp; model%mw=783._wp; model%mr=763._wp
      model%gs0=10.72854_wp; model%gw0=13.29015_wp; model%gr0=3.66098_wp
      model%as=1.365469_wp; model%bs=0.226061_wp; model%cs=0.409704_wp; model%ds=0.901995_wp
      model%aw=1.402488_wp; model%bw=0.172577_wp; model%cw=0.344293_wp; model%dw=0.983955_wp
      model%ar=0.515_wp
    case default
      return
    end select
    status='ok'
  end subroutine load_ddrmf_model

  subroutine coupling_values(model, rho_b, gs, gw, gr, dgs, dgw, dgr)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in) :: rho_b
    real(wp), intent(out) :: gs, gw, gr, dgs, dgw, dgr
    real(wp) :: x, us, uw, fs, fw, dfsdx, dfwdx
    x=rho_b/model%rho_sat
    us=x+model%ds
    uw=x+model%dw
    fs=model%as*(1._wp+model%bs*us*us)/(1._wp+model%cs*us*us)
    fw=model%aw*(1._wp+model%bw*uw*uw)/(1._wp+model%cw*uw*uw)
    dfsdx=2._wp*model%as*(model%bs-model%cs)*us/(1._wp+model%cs*us*us)**2
    dfwdx=2._wp*model%aw*(model%bw-model%cw)*uw/(1._wp+model%cw*uw*uw)**2
    gs=model%gs0*fs
    gw=model%gw0*fw
    gr=model%gr0*exp(-model%ar*(x-1._wp))
    dgs=model%gs0*dfsdx/model%rho_sat
    dgw=model%gw0*dfwdx/model%rho_sat
    dgr=-model%ar*gr/model%rho_sat
  end subroutine coupling_values

  subroutine kinetic_thermo(rho, mass, eps, press, scalar, ef)
    real(wp), intent(in) :: rho, mass
    real(wp), intent(out) :: eps, press, scalar, ef
    real(wp) :: kf, x, root, asinhx, pref
    eps=0._wp; press=0._wp; scalar=0._wp; ef=mass
    if(rho<=0._wp) return
    kf=hbarc*(3._wp*pi*pi*rho)**(1._wp/3._wp)
    ef=sqrt(kf*kf+mass*mass)
    if(mass<=0._wp) return
    x=kf/mass
    root=sqrt(1._wp+x*x)
    asinhx=log(x+root)
    pref=mass**4/(8._wp*pi*pi*hbarc**3)
    eps=pref*(x*(2._wp*x*x+1._wp)*root-asinhx)
    press=(pref/3._wp)*(x*(2._wp*x*x-3._wp)*root+3._wp*asinhx)
    scalar=mass**3/(2._wp*pi*pi*hbarc**3)*(x*root-asinhx)
  end subroutine kinetic_thermo

  subroutine solve_ddrmf_homogeneous(model, rho_b, yp, state, status)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in) :: rho_b, yp
    type(ddrmf_state_t), intent(out) :: state
    character(len=*), intent(out) :: status
    real(wp) :: rn, rp, r3, gs, gw, gr, dgs, dgw, dgr
    real(wp) :: cs2, s_lo, s_hi, s_mid, f_mid, nsn, nsp, en, ep, pn, pp
    real(wp) :: sn_tmp, sp_tmp, ef_n, ef_p, vn, vp, sigma_r
    integer :: it

    state=ddrmf_state_t()
    status='invalid density or proton fraction'
    if(.not.ieee_is_finite(rho_b+yp).or.rho_b<=0._wp) return
    if(yp<0._wp.or.yp>1._wp) return
    rn=rho_b*(1._wp-yp); rp=rho_b*yp; r3=rp-rn
    call coupling_values(model,rho_b,gs,gw,gr,dgs,dgw,dgr)
    if(min(gs,gw,gr,model%ms,model%mw,model%mr)<=0._wp) then
      status='nonpositive DD-RMF coupling or meson mass'; return
    endif

    cs2=(gs/model%ms)**2*hbarc**3
    s_lo=0._wp; s_hi=min(model%mn,model%mp)-1.e-9_wp
    do it=1,160
      s_mid=0.5_wp*(s_lo+s_hi)
      call kinetic_thermo(rn,model%mn-s_mid,en,pn,sn_tmp,ef_n)
      call kinetic_thermo(rp,model%mp-s_mid,ep,pp,sp_tmp,ef_p)
      f_mid=s_mid-cs2*(sn_tmp+sp_tmp)
      if(f_mid>0._wp) then
        s_hi=s_mid
      else
        s_lo=s_mid
      endif
      if(abs(s_hi-s_lo)<=2.e-13_wp*max(1._wp,s_mid)) exit
    enddo
    s_mid=0.5_wp*(s_lo+s_hi)
    call kinetic_thermo(rn,model%mn-s_mid,en,pn,nsn,ef_n)
    call kinetic_thermo(rp,model%mp-s_mid,ep,pp,nsp,ef_p)
    state%scalar_field_residual=s_mid-cs2*(nsn+nsp)
    vn=(gw/model%mw)**2*hbarc**3*rho_b
    vp=(gr/model%mr)**2*hbarc**3*r3
    sigma_r=(dgw/gw)*vn*rho_b+(dgr/gr)*vp*r3-(dgs/gs)*s_mid*(nsn+nsp)
    state%mu_n=ef_n+vn-vp+sigma_r
    state%mu_p=ef_p+vn+vp+sigma_r
    state%eps_b=en+ep+0.5_wp*s_mid*(nsn+nsp)+0.5_wp*vn*rho_b+0.5_wp*vp*r3
    state%p_b=pn+pp-0.5_wp*s_mid*(nsn+nsp)+0.5_wp*vn*rho_b+0.5_wp*vp*r3+rho_b*sigma_r
    state%p_thermo=state%mu_n*rn+state%mu_p*rp-state%eps_b
    state%pressure_identity_residual=state%p_thermo-state%p_b
    state%rho_b=rho_b; state%yp=yp; state%rho_n=rn; state%rho_p=rp; state%rho_3=r3
    state%mstar_n=model%mn-s_mid; state%mstar_p=model%mp-s_mid
    state%sigma_scalar=s_mid; state%rho_scalar=nsn+nsp; state%sigma_r=sigma_r
    if(.not.all(ieee_is_finite([state%mu_n,state%mu_p,state%eps_b,state%p_b,state%sigma_r]))) then
      status='nonfinite DD-RMF homogeneous state'; return
    endif
    if(abs(state%scalar_field_residual)>1.e-8_wp.or.abs(state%pressure_identity_residual)>2.e-7_wp) then
      status='DD-RMF field or pressure identity residual exceeds tolerance'; return
    endif
    status='ok'
  end subroutine solve_ddrmf_homogeneous

  subroutine lepton_state(rho_p, rho_e, rho_mu, mu_l, eps_l, p_l, muon_active, rho_p_onset)
    real(wp), intent(in) :: rho_p
    real(wp), intent(out) :: rho_e, rho_mu, mu_l, eps_l, p_l, rho_p_onset
    logical, intent(out) :: muon_active
    real(wp) :: mu_lo, mu_hi, mu_mid, ne, nmu, dummy, ee, pe, se, efe, emu, pmu, smu, efmu
    integer :: it
    rho_e=0._wp; rho_mu=0._wp; eps_l=0._wp; p_l=0._wp; muon_active=.false.
    mu_l=m_e
    rho_p_onset=(sqrt(m_mu*m_mu-m_e*m_e)/hbarc)**3/(3._wp*pi*pi)
    if(rho_p<=0._wp) return
    mu_hi=sqrt((hbarc*(3._wp*pi*pi*rho_p)**(1._wp/3._wp))**2+m_e*m_e)
    if(mu_hi<=m_mu) then
      mu_l=mu_hi; rho_e=rho_p; rho_mu=0._wp
    else
      mu_lo=m_mu
      do it=1,140
        mu_mid=0.5_wp*(mu_lo+mu_hi)
        ne=(sqrt(max(0._wp,mu_mid*mu_mid-m_e*m_e))/hbarc)**3/(3._wp*pi*pi)
        nmu=(sqrt(max(0._wp,mu_mid*mu_mid-m_mu*m_mu))/hbarc)**3/(3._wp*pi*pi)
        if(ne+nmu<rho_p) then
          mu_lo=mu_mid
        else
          mu_hi=mu_mid
        endif
        if(abs(mu_hi-mu_lo)<2.e-13_wp*max(1._wp,mu_mid)) exit
      enddo
      mu_l=0.5_wp*(mu_lo+mu_hi)
      rho_e=(sqrt(max(0._wp,mu_l*mu_l-m_e*m_e))/hbarc)**3/(3._wp*pi*pi)
      rho_mu=(sqrt(max(0._wp,mu_l*mu_l-m_mu*m_mu))/hbarc)**3/(3._wp*pi*pi)
      muon_active=rho_mu>1.e-16_wp
    endif
    call kinetic_thermo(rho_e,m_e,ee,pe,se,efe)
    call kinetic_thermo(rho_mu,m_mu,emu,pmu,smu,efmu)
    eps_l=ee+emu; p_l=pe+pmu
  end subroutine lepton_state

  subroutine beta_residual(model,rho_b,yp,residual,state,rho_e,rho_mu,mu_l,eps_l,p_l,active,rho_p_onset,status)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in) :: rho_b, yp
    real(wp), intent(out) :: residual, rho_e, rho_mu, mu_l, eps_l, p_l, rho_p_onset
    type(ddrmf_state_t), intent(out) :: state
    logical, intent(out) :: active
    character(len=*), intent(out) :: status
    call solve_ddrmf_homogeneous(model,rho_b,yp,state,status)
    if(status/='ok') return
    call lepton_state(rho_b*yp,rho_e,rho_mu,mu_l,eps_l,p_l,active,rho_p_onset)
    residual=state%mu_n-state%mu_p-mu_l
    if(.not.ieee_is_finite(residual)) status='nonfinite beta-equilibrium residual'
  end subroutine beta_residual

  subroutine solve_ddrmf_beta_equilibrium(model,rho_b,result,status)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in) :: rho_b
    type(ddrmf_beta_result_t), intent(out) :: result
    character(len=*), intent(out) :: status
    real(wp) :: lo,hi,mid, flo,fhi,fmid,re,rm,ml,el,pl,rponset
    type(ddrmf_state_t) :: st
    logical :: active
    integer :: it
    result=ddrmf_beta_result_t()
    status='invalid beta-equilibrium density'
    if(.not.ieee_is_finite(rho_b).or.rho_b<=0._wp) return
    lo=1.e-10_wp; hi=0.5_wp
    call beta_residual(model,rho_b,lo,flo,st,re,rm,ml,el,pl,active,rponset,status)
    if(status/='ok') return
    call beta_residual(model,rho_b,hi,fhi,st,re,rm,ml,el,pl,active,rponset,status)
    if(status/='ok') return
    if(flo*fhi>0._wp) then
      status='beta-equilibrium root is not bracketed on 1e-10 <= Yp <= 0.5'; return
    endif
    do it=1,140
      mid=0.5_wp*(lo+hi)
      call beta_residual(model,rho_b,mid,fmid,st,re,rm,ml,el,pl,active,rponset,status)
      if(status/='ok') return
      if(flo*fmid>0._wp) then
        lo=mid; flo=fmid
      else
        hi=mid; fhi=fmid
      endif
      if(abs(hi-lo)<2.e-13_wp) exit
    enddo
    mid=0.5_wp*(lo+hi)
    call beta_residual(model,rho_b,mid,fmid,st,re,rm,ml,el,pl,active,rponset,status)
    if(status/='ok') return
    result%baryons=st
    result%yp=mid; result%rho_e=re; result%rho_mu=rm; result%mu_e=ml
    result%muon_active=active
    if(active) then
      result%mu_mu=ml
    else
      result%mu_mu=-1._wp
    endif
    result%ye=re/rho_b; result%ymu=rm/rho_b
    result%eps_leptons=el; result%p_leptons=pl
    result%eps_total=st%eps_b+el; result%p_total=st%p_b+pl
    result%beta_residual=fmid
    result%charge_residual=rho_b*mid-re-rm
    result%muon_onset_rho_p=rponset
    result%muon_onset_gap=max(0._wp,m_mu-ml)
    if(abs(result%beta_residual)>2.e-8_wp.or.abs(result%charge_residual)>2.e-12_wp) then
      status='beta-equilibrium or charge-neutrality residual exceeds tolerance'; return
    endif
    status='ok'
  end subroutine solve_ddrmf_beta_equilibrium

  subroutine nuclear_binding_per_baryon(model,rho_b,yp,value,status)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in) :: rho_b,yp
    real(wp), intent(out) :: value
    character(len=*), intent(out) :: status
    type(ddrmf_state_t) :: st
    real(wp) :: rest
    call solve_ddrmf_homogeneous(model,rho_b,yp,st,status)
    value=0._wp
    if(status/='ok') return
    rest=(1._wp-yp)*model%mn+yp*model%mp
    value=st%eps_b/rho_b-rest
  end subroutine nuclear_binding_per_baryon

  subroutine symmetry_energy_ddrmf(model,rho_b,s2,sfd,status)
    type(ddrmf_model_t), intent(in) :: model
    real(wp), intent(in) :: rho_b
    real(wp), intent(out) :: s2,sfd
    character(len=*), intent(out) :: status
    real(wp), parameter :: h=0.05_wp
    real(wp) :: e0,ep1,em1,ep2,em2,epnm
    s2=0._wp; sfd=0._wp
    call nuclear_binding_per_baryon(model,rho_b,0.5_wp,e0,status); if(status/='ok') return
    call nuclear_binding_per_baryon(model,rho_b,(1._wp-h)/2._wp,ep1,status); if(status/='ok') return
    call nuclear_binding_per_baryon(model,rho_b,(1._wp+h)/2._wp,em1,status); if(status/='ok') return
    call nuclear_binding_per_baryon(model,rho_b,(1._wp-2._wp*h)/2._wp,ep2,status); if(status/='ok') return
    call nuclear_binding_per_baryon(model,rho_b,(1._wp+2._wp*h)/2._wp,em2,status); if(status/='ok') return
    call nuclear_binding_per_baryon(model,rho_b,0._wp,epnm,status); if(status/='ok') return
    s2=(-ep2+16._wp*ep1-30._wp*e0+16._wp*em1-em2)/(24._wp*h*h)
    sfd=epnm-e0
    status='ok'
  end subroutine symmetry_energy_ddrmf

  subroutine evaluate_ddrmf(model_id,rho_b,diagnostic,status)
    character(len=*), intent(in) :: model_id
    real(wp), intent(in) :: rho_b
    type(ddrmf_diagnostic_t), intent(out) :: diagnostic
    character(len=*), intent(out) :: status
    type(ddrmf_model_t) :: model
    diagnostic=ddrmf_diagnostic_t()
    call load_ddrmf_model(model_id,model,status); if(status/='ok') return
    call solve_ddrmf_beta_equilibrium(model,rho_b,diagnostic%beta,status); if(status/='ok') return
    call symmetry_energy_ddrmf(model,rho_b,diagnostic%symmetry_second_derivative, &
      diagnostic%symmetry_pnm_minus_snm,status)
  end subroutine evaluate_ddrmf

end module compose_eos_rmf
