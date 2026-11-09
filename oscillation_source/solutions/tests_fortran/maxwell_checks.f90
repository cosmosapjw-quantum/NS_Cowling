program maxwell_checks
    use ns_precision, only: wp
    use ns_status
    use ns_physical_constants, only: two_pi
    use ns_eos_types, only: eos_state_t,SIDE_REGULAR
    use ns_background_records, only: star_options_t,star_profile_t
    use ns_maxwell_tov, only: maxwell_star_t,integrate_maxwell_tov
    use ns_maxwell_cowling, only: maxwell_cowling_shot_t,shoot_maxwell_cowling
    use ns_cowling, only: cowling_options_t
    use ns_mode_spectrum, only: mode_options_t,mode_spectrum_t
    use ns_maxwell_mode_spectrum, only: maxwell_mode_spectrum_t
    use ns_production_solvers, only: solve_tov_star,solve_cowling_spectrum
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(wp),parameter :: mev_fm3_to_j_m3=1.602176634e32_wp
    real(wp),parameter :: mev_to_j=1.602176634e-13_wp
    real(wp),parameter :: fm3_to_m3=1e45_wp
    real(wp),parameter :: pt_mev_fm3=1._wp,pc_mev_fm3=1.4_wp,pcut_mev_fm3=0.2_wp
    real(wp) :: jump_eta,pt,pc,pcut,rel_m,rel_r,dp_rel,freq_grid(3)
    integer :: checks
    type(star_options_t) :: star_options
    type(maxwell_star_t) :: star
    type(star_profile_t) :: smooth_profile
    type(cowling_options_t) :: cowling_options
    type(maxwell_cowling_shot_t) :: slow_shot,rapid_shot
    type(mode_options_t) :: mode_options
    type(mode_spectrum_t) :: smooth_spectrum
    type(maxwell_mode_spectrum_t) :: slow_spectrum,rapid_spectrum
    type(status_t) :: status

    checks=0
    pt=pt_mev_fm3*mev_fm3_to_j_m3
    pc=pc_mev_fm3*mev_fm3_to_j_m3
    pcut=pcut_mev_fm3*mev_fm3_to_j_m3
    star_options%r0_m=1._wp
    star_options%r_limit_m=2e6_wp
    star_options%p_cut_pa=pcut
    star_options%nu_start=0._wp
    star_options%surface_policy='finite_cutoff_schwarzschild_time'
    star_options%eos_id='analytic_quadratic_cold_maxwell_fixture'
    star_options%ode%atol=[1e-5_wp,1e13_wp,1e-12_wp]
    star_options%ode%rtol=2e-9_wp
    star_options%ode%h_initial=1e-5_wp
    star_options%ode%h_min=1e-15_wp
    star_options%ode%h_max=0.05_wp
    star_options%ode%max_steps=100000

    ! A finite but small first-order density jump, below the Seidov rapid-onset
    ! threshold. This fixture tests a true nonzero discontinuity.
    jump_eta=5e-5_wp
    call solve_tov_star(eos_high,eos_low,pc,pt,star_options,star)
    if(.not.star%valid.or..not.is_ok(star%status)) then
        write(*,'(A,A,A,A,A,A)') 'MAXWELL_TOV_DIAGNOSTIC status=',trim(status_name(star%status%code)), &
            ' category=',trim(star%status%category),' detail=',trim(star%status%detail)
    end if
    call require(star%valid.and.is_ok(star%status),'native phase-segmented TOV')
    call require(star%core%p_pa(star%core%n)==pt,'core terminates at exact Pt')
    call require(star%envelope%p_pa(1)==pt,'envelope restarts at exact Pt')
    call require(star%core%r_m(star%core%n)==star%envelope%r_m(1),'radius continuity')
    call require(star%core%m_kg(star%core%n)==star%envelope%m_kg(1),'mass continuity')
    call require(star%jump%delta_density_m3>0._wp,'positive Maxwell density jump')
    call require(abs(star%jump%chemical_potential_mismatch_j/star%jump%chemical_potential_j)<1e-10_wp, &
        'chemical-potential continuity')
    call require(abs(star%jump%thermodynamic_residual_j_m3)/ &
        max(abs(star%jump%delta_energy_j_m3),tiny(1._wp))<1e-9_wp,'cold first-law jump identity')
    call require(star%jump%rapid_onset_stable_diagnostic,'small-jump Seidov onset diagnostic')

    cowling_options%ell=2
    cowling_options%amplitude=1._wp
    cowling_options%nu_shift=0._wp
    cowling_options%ode%atol=[1e-12_wp,1e-12_wp]
    cowling_options%ode%rtol=2e-8_wp
    cowling_options%ode%h_initial=1e-5_wp
    cowling_options%ode%h_min=1e-15_wp
    cowling_options%ode%h_max=0.025_wp
    cowling_options%ode%max_steps=100000
    call shoot_maxwell_cowling(star,two_pi*1800._wp,'slow',cowling_options,slow_shot)
    call require(slow_shot%valid.and.is_ok(slow_shot%status),'native slow-conversion Cowling shot')
    call shoot_maxwell_cowling(star,two_pi*1800._wp,'rapid',cowling_options,rapid_shot)
    call require(rapid_shot%valid.and.is_ok(rapid_shot%status),'native rapid-conversion Cowling shot')
    dp_rel=abs(slow_shot%delta_p_inner_pa-slow_shot%delta_p_outer_pa)/ &
        max(abs(slow_shot%delta_p_inner_pa),abs(slow_shot%delta_p_outer_pa),tiny(1._wp))
    call require(dp_rel<1e-9_wp,'slow interface conserves Delta P')
    dp_rel=abs(rapid_shot%delta_p_inner_pa-rapid_shot%delta_p_outer_pa)/ &
        max(abs(rapid_shot%delta_p_inner_pa),abs(rapid_shot%delta_p_outer_pa),tiny(1._wp))
    call require(dp_rel<1e-9_wp,'rapid interface conserves Delta P')
    call require(size(slow_shot%r_m)>=4.and.size(rapid_shot%r_m)>=4,'piecewise Cowling profiles retained')

    ! The production spectrum facade dispatches smooth profiles and a true
    ! phase-segmented Maxwell star without flattening the interface.
    mode_options%shot=cowling_options
    mode_options%root%xtol=1e-9_wp
    mode_options%root%max_iter=100
    mode_options%residual_tol=1e-5_wp
    mode_options%zero_tol=1e-10_wp
    freq_grid=two_pi*[1799._wp,1800._wp,1801._wp]
    call solve_cowling_spectrum(star,freq_grid,mode_options,'slow',slow_spectrum)
    call require(slow_spectrum%scan%valid.and.slow_spectrum%valid, &
        'production Maxwell spectrum dispatch completes slow-conversion scan')
    call require(size(slow_spectrum%scan%trials)==3.and.all(slow_spectrum%scan%trials%valid), &
        'production Maxwell spectrum retains all successful slow trials')
    call solve_cowling_spectrum(star,freq_grid,mode_options,'rapid',rapid_spectrum)
    call require(rapid_spectrum%scan%valid.and.rapid_spectrum%valid, &
        'production Maxwell spectrum dispatch completes rapid-conversion scan')

    ! Small-jump smooth-limit cross-check against the pre-existing single-EOS
    ! Fortran TOV implementation; the two phase callbacks tend continuously
    ! to the same cold barotrope as eta -> 0.
    jump_eta=1e-7_wp
    call integrate_maxwell_tov(eos_high,eos_low,pc,pt,star_options,star)
    call require(star%valid.and.is_ok(star%status),'small-jump native Maxwell TOV')
    call solve_tov_star(eos_smooth,pc,star_options,smooth_profile)
    call require(smooth_profile%valid.and.is_ok(smooth_profile%status),'single-EOS smooth TOV reference')

    call solve_cowling_spectrum(smooth_profile,freq_grid,mode_options,smooth_spectrum)
    call require(smooth_spectrum%scan%valid.and.smooth_spectrum%scan%all_trials_succeeded, &
        'production Cowling spectrum accepts smooth TOV profile')
    rel_m=abs(star%mass_kg-smooth_profile%mass_kg)/smooth_profile%mass_kg
    rel_r=abs(star%radius_m-smooth_profile%radius_m)/smooth_profile%radius_m
    call require(rel_m<2e-5_wp,'Maxwell TOV smooth-limit mass convergence')
    call require(rel_r<2e-5_wp,'Maxwell TOV smooth-limit radius convergence')

    print '(A,I0,A,ES12.4,A,ES12.4)', 'maxwell_checks: ',checks,' checks passed; smooth-limit dM/M=',rel_m,', dR/R=',rel_r
contains
    subroutine require(ok,label)
        logical,intent(in) :: ok
        character(len=*),intent(in) :: label
        if(.not.ok) then
            write(*,'(A)') 'FAIL: '//trim(label)
            stop 1
        end if
        checks=checks+1
    end subroutine

    subroutine eos_low(p,side,state,st)
        real(wp),intent(in) :: p
        integer,intent(in) :: side
        type(eos_state_t),intent(out) :: state
        type(status_t),intent(out) :: st
        call quadratic_eos(p,1._wp,939._wp,side,state,st)
    end subroutine

    subroutine eos_high(p,side,state,st)
        real(wp),intent(in) :: p
        integer,intent(in) :: side
        type(eos_state_t),intent(out) :: state
        type(status_t),intent(out) :: st
        real(wp) :: nt,stiffness,rest_mev
        nt=1._wp+jump_eta
        stiffness=1._wp/(nt*nt)
        rest_mev=941._wp-2._wp/nt
        call quadratic_eos(p,stiffness,rest_mev,side,state,st)
    end subroutine

    subroutine eos_smooth(p,side,state,st)
        real(wp),intent(in) :: p
        integer,intent(in) :: side
        type(eos_state_t),intent(out) :: state
        type(status_t),intent(out) :: st
        call quadratic_eos(p,1._wp,939._wp,side,state,st)
    end subroutine

    subroutine quadratic_eos(p,k_mev_fm3,rest_mev,side,state,st)
        real(wp),intent(in) :: p,k_mev_fm3,rest_mev
        integer,intent(in) :: side
        type(eos_state_t),intent(out) :: state
        type(status_t),intent(out) :: st
        real(wp) :: p_mev_fm3,n_fm3,mu_mev,eps_mev
        state=eos_state_t()
        call fail_status(st,NS_OUT_OF_DOMAIN,'physics','analytic Maxwell EOS requires positive pressure')
        if(p<=0._wp.or.k_mev_fm3<=0._wp) return
        p_mev_fm3=p/mev_fm3_to_j_m3
        n_fm3=sqrt(p_mev_fm3/k_mev_fm3)
        mu_mev=rest_mev+2._wp*k_mev_fm3*n_fm3
        eps_mev=mu_mev*n_fm3-p_mev_fm3
        if(.not.all(ieee_is_finite([p_mev_fm3,n_fm3,mu_mev,eps_mev])).or. &
           min(n_fm3,mu_mev,eps_mev)<=0._wp) then
            call fail_status(st,NS_NONFINITE_RESULT,'physics','analytic Maxwell EOS returned invalid state')
            return
        end if
        state%pressure_pa=p
        state%energy_j_m3=eps_mev*mev_fm3_to_j_m3
        state%n_b_m3=n_fm3*fm3_to_m3
        state%mu_j=mu_mev*mev_to_j
        state%chi=2._wp*k_mev_fm3*n_fm3/mu_mev
        state%side=side
        state%valid=.true.
        call clear_status(st)
    end subroutine
end program maxwell_checks
