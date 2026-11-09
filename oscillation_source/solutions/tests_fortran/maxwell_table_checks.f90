program maxwell_table_checks
    use ns_precision, only: wp
    use ns_status
    use ns_physical_constants, only: two_pi
    use ns_eos_types, only: eos_state_t,SIDE_INNER,SIDE_OUTER
    use ns_table_io, only: NORMALIZED_TABLE_HEADER,normalized_table_t,load_normalized_table
    use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos,eos_table_at_pressure
    use ns_background_records, only: star_options_t
    use ns_maxwell_tov, only: maxwell_star_t,integrate_maxwell_tov
    use ns_maxwell_cowling, only: maxwell_cowling_shot_t,shoot_maxwell_cowling
    use ns_cowling, only: cowling_options_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,parameter :: nrows=24
    real(wp),parameter :: mev_fm3_to_j_m3=1.602176634e32_wp
    real(wp),parameter :: mev_to_j=1.602176634e-13_wp
    real(wp),parameter :: fm3_to_m3=1e45_wp
    real(wp),parameter :: pt_mev_fm3=1._wp,pc_mev_fm3=1.4_wp,pcut_mev_fm3=0.2_wp
    real(wp),parameter :: eta=5e-5_wp
    character(len=*),parameter :: core_file='build/maxwell_core_fixture.csv'
    character(len=*),parameter :: envelope_file='build/maxwell_envelope_fixture.csv'
    real(wp) :: pt,pc,pcut,k_high,rest_high,mu_mismatch,thermo_res,thermo_scale,dp_rel
    integer :: checks
    type(normalized_table_t) :: core_raw,envelope_raw
    type(table_eos_t) :: core_model,envelope_model
    type(table_audit_t) :: core_audit,envelope_audit
    type(eos_state_t) :: core_t,env_t
    type(star_options_t) :: star_options
    type(maxwell_star_t) :: star
    type(cowling_options_t) :: cowling_options
    type(maxwell_cowling_shot_t) :: slow_shot,rapid_shot
    type(status_t) :: status

    checks=0
    pt=pt_mev_fm3*mev_fm3_to_j_m3
    pc=pc_mev_fm3*mev_fm3_to_j_m3
    pcut=pcut_mev_fm3*mev_fm3_to_j_m3
    k_high=1._wp/(1._wp+eta)**2
    rest_high=941._wp-2._wp/(1._wp+eta)

    call write_phase_table(envelope_file,1._wp,939._wp,sqrt(0.05_wp),1._wp)
    call write_phase_table(core_file,k_high,rest_high,1._wp+eta,(1._wp+eta)*sqrt(2._wp))
    call load_normalized_table(envelope_file,envelope_raw,status)
    call require(is_ok(status).and.envelope_raw%ready,'normalized envelope table loads')
    call load_normalized_table(core_file,core_raw,status)
    call require(is_ok(status).and.core_raw%ready,'normalized core table loads')
    call init_table_eos(envelope_raw,envelope_model,envelope_audit,status)
    call require(is_ok(status).and.envelope_model%ready.and.envelope_audit%fully_admitted, &
        'low-density table provider admitted')
    call init_table_eos(core_raw,core_model,core_audit,status)
    call require(is_ok(status).and.core_model%ready.and.core_audit%fully_admitted, &
        'high-density table provider admitted')

    call eos_table_at_pressure(core_model,pt,SIDE_INNER,core_t,status)
    call require(is_ok(status).and.core_t%valid,'high-density one-sided Pt query')
    call eos_table_at_pressure(envelope_model,pt,SIDE_OUTER,env_t,status)
    call require(is_ok(status).and.env_t%valid,'low-density one-sided Pt query')
    call require(abs(core_t%pressure_pa/pt-1._wp)<1e-14_wp,'core table pressure at Pt')
    call require(abs(env_t%pressure_pa/pt-1._wp)<1e-14_wp,'envelope table pressure at Pt')
    mu_mismatch=core_t%mu_j-env_t%mu_j
    call require(abs(mu_mismatch)/max(abs(core_t%mu_j),abs(env_t%mu_j))<1e-10_wp, &
        'normalized tables satisfy chemical-potential coexistence')
    call require(core_t%energy_j_m3>env_t%energy_j_m3.and.core_t%n_b_m3>env_t%n_b_m3, &
        'normalized tables retain positive sharp jumps')
    thermo_res=(core_t%energy_j_m3-env_t%energy_j_m3)- &
        0.5_wp*(core_t%mu_j+env_t%mu_j)*(core_t%n_b_m3-env_t%n_b_m3)
    thermo_scale=max(abs(core_t%energy_j_m3-env_t%energy_j_m3), &
        abs(0.5_wp*(core_t%mu_j+env_t%mu_j)*(core_t%n_b_m3-env_t%n_b_m3)),tiny(1._wp))
    call require(abs(thermo_res)/thermo_scale<1e-9_wp,'normalized-table cold first-law jump closure')

    star_options%r0_m=1._wp
    star_options%r_limit_m=2e6_wp
    star_options%p_cut_pa=pcut
    star_options%nu_start=0._wp
    star_options%surface_policy='finite_cutoff_schwarzschild_time'
    star_options%eos_id='normalized_table_analytic_maxwell_fixture'
    star_options%ode%atol=[1e-5_wp,1e13_wp,1e-12_wp]
    star_options%ode%rtol=2e-9_wp
    star_options%ode%h_initial=1e-5_wp
    star_options%ode%h_min=1e-15_wp
    star_options%ode%h_max=0.05_wp
    star_options%ode%max_steps=100000
    call integrate_maxwell_tov(high_eos,low_eos,pc,pt,star_options,star)
    call require(star%valid.and.is_ok(star%status),'table-backed phase-segmented TOV')
    call require(star%core%p_pa(star%core%n)==pt,'table-backed core ends at exact Pt')
    call require(star%envelope%p_pa(1)==pt,'table-backed envelope starts at exact Pt')
    call require(star%core%r_m(star%core%n)==star%envelope%r_m(1), &
        'table-backed interface radius continuity')
    call require(star%core%m_kg(star%core%n)==star%envelope%m_kg(1), &
        'table-backed interface mass continuity')

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
    call require(slow_shot%valid.and.is_ok(slow_shot%status),'table-backed slow Cowling shot')
    call shoot_maxwell_cowling(star,two_pi*1800._wp,'rapid',cowling_options,rapid_shot)
    call require(rapid_shot%valid.and.is_ok(rapid_shot%status),'table-backed rapid Cowling shot')
    dp_rel=abs(slow_shot%delta_p_inner_pa-slow_shot%delta_p_outer_pa)/ &
        max(abs(slow_shot%delta_p_inner_pa),abs(slow_shot%delta_p_outer_pa),tiny(1._wp))
    call require(dp_rel<1e-9_wp,'table-backed slow Delta P continuity')
    dp_rel=abs(rapid_shot%delta_p_inner_pa-rapid_shot%delta_p_outer_pa)/ &
        max(abs(rapid_shot%delta_p_inner_pa),abs(rapid_shot%delta_p_outer_pa),tiny(1._wp))
    call require(dp_rel<1e-9_wp,'table-backed rapid Delta P continuity')

    print '(A,I0,A)', 'maxwell_table_checks: ',checks,' checks passed'
    call delete_fixture(core_file)
    call delete_fixture(envelope_file)
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

    subroutine write_phase_table(filename,k_mev_fm3,rest_mev,n_min,n_max)
        character(len=*),intent(in) :: filename
        real(wp),intent(in) :: k_mev_fm3,rest_mev,n_min,n_max
        integer :: u,ios,i
        real(wp) :: n_fm3,p_mev_fm3,eps_mev,mu_mev,v(15),log_lo,log_hi
        open(newunit=u,file=filename,status='replace',action='write',iostat=ios)
        call require(ios==0,'open normalized fixture '//trim(filename))
        write(u,'(A)',iostat=ios) trim(NORMALIZED_TABLE_HEADER)
        call require(ios==0,'write normalized fixture header')
        log_lo=log(n_min);log_hi=log(n_max)
        do i=1,nrows
            n_fm3=exp(log_lo+(log_hi-log_lo)*real(i-1,wp)/real(nrows-1,wp))
            p_mev_fm3=k_mev_fm3*n_fm3**2
            eps_mev=rest_mev*n_fm3+p_mev_fm3
            mu_mev=rest_mev+2._wp*k_mev_fm3*n_fm3
            v=[n_fm3,n_fm3*fm3_to_m3,p_mev_fm3,p_mev_fm3*mev_fm3_to_j_m3, &
                eps_mev,eps_mev*mev_fm3_to_j_m3,mu_mev,mu_mev*mev_to_j, &
                p_mev_fm3/n_fm3,0._wp,mu_mev/939._wp-1._wp,0._wp,0._wp,0._wp,0._wp]
            write(u,'(I0,3(",",I0),15(",",ES25.16E3))',iostat=ios) i,0,i,1,v
            call require(ios==0,'write normalized fixture row')
        end do
        close(u,iostat=ios)
        call require(ios==0,'close normalized fixture')
    end subroutine

    subroutine delete_fixture(filename)
        character(len=*),intent(in) :: filename
        integer :: u,ios
        open(newunit=u,file=filename,status='old',iostat=ios)
        if(ios==0) close(u,status='delete')
    end subroutine

    subroutine high_eos(p,side,state,st)
        real(wp),intent(in) :: p
        integer,intent(in) :: side
        type(eos_state_t),intent(out) :: state
        type(status_t),intent(out) :: st
        call eos_table_at_pressure(core_model,p,side,state,st)
    end subroutine

    subroutine low_eos(p,side,state,st)
        real(wp),intent(in) :: p
        integer,intent(in) :: side
        type(eos_state_t),intent(out) :: state
        type(status_t),intent(out) :: st
        call eos_table_at_pressure(envelope_model,p,side,state,st)
    end subroutine
end program maxwell_table_checks