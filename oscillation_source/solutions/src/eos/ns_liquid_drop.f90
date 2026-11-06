! T12 pedagogical liquid-drop crust toy at fixed baryon density.
! Scope: one nuclear species + neutral zero-temperature electron background below neutron drip.
! Not included: pairing, lattice energy, neutron drip/free neutrons, shell corrections, pasta, full-star EOS dispatch.
module ns_liquid_drop
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use ns_precision, only: wp
    use ns_physical_constants, only: inverse_fm3_in_si, mev_in_j
    use ns_electron_gas, only: electron_state_t, electron_state
    use ns_status
    implicit none
    private

    type, public :: liquid_drop_model_t
        logical :: ready = .false.
        real(wp) :: a_volume_mev = 15.75_wp
        real(wp) :: a_surface_mev = 17.8_wp
        real(wp) :: a_coulomb_mev = 0.711_wp
        real(wp) :: a_asymmetry_mev = 23.7_wp
        real(wp) :: n_b_min_fm3 = 1.0e-10_wp
        real(wp) :: n_b_max_fm3 = 2.5e-4_wp
    end type

    type, public :: liquid_drop_terms_t
        logical :: valid = .false.
        real(wp) :: volume_mev = 0._wp
        real(wp) :: surface_mev = 0._wp
        real(wp) :: coulomb_mev = 0._wp
        real(wp) :: asymmetry_mev = 0._wp
        real(wp) :: total_mev = 0._wp
    end type

    type, public :: liquid_drop_state_t
        logical :: valid = .false.
        integer :: a = 0, z = 0
        real(wp) :: n_b_m3 = 0._wp, n_b_fm3 = 0._wp
        real(wp) :: n_nucleus_m3 = 0._wp, n_e_m3 = 0._wp
        real(wp) :: volume_mev_per_baryon = 0._wp
        real(wp) :: surface_mev_per_baryon = 0._wp
        real(wp) :: coulomb_mev_per_baryon = 0._wp
        real(wp) :: asymmetry_mev_per_baryon = 0._wp
        real(wp) :: nuclear_mev_per_baryon = 0._wp
        real(wp) :: electron_mev_per_baryon = 0._wp
        real(wp) :: total_mev_per_baryon = 0._wp
        logical :: charge_neutrality_closed = .false.
    end type

    type, public :: liquid_drop_scan_t
        logical :: valid = .false.
        integer :: candidate_count = 0
        integer :: best_a = 0, best_z = 0
        logical :: minimum_on_boundary = .false.
        logical :: stationary_point_claimed = .false.
        logical :: mesh_is_discrete = .true.
        type(liquid_drop_state_t) :: best
    end type

    public :: init_liquid_drop_model, liquid_drop_nuclear_terms
    public :: evaluate_liquid_drop_candidate, scan_liquid_drop_candidates

contains

    subroutine init_liquid_drop_model(model,status)
        type(liquid_drop_model_t),intent(out)::model
        type(status_t),intent(out)::status
        model=liquid_drop_model_t()
        model%ready=.true.
        call clear_status(status)
    end subroutine

    subroutine liquid_drop_nuclear_terms(model,a,z,terms,status)
        type(liquid_drop_model_t),intent(in)::model
        integer,intent(in)::a,z
        type(liquid_drop_terms_t),intent(out)::terms
        type(status_t),intent(out)::status
        real(wp)::ar
        terms=liquid_drop_terms_t()
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','liquid-drop model not ready');return
        end if
        if(a<=0.or.z<0.or.z>a) then
            call fail_status(status,NS_INVALID_INPUT,'physics','liquid-drop candidate requires A>0 and 0<=Z<=A');return
        end if
        ar=real(a,wp)
        terms%volume_mev=-model%a_volume_mev*ar
        terms%surface_mev=model%a_surface_mev*ar**(2._wp/3._wp)
        terms%coulomb_mev=model%a_coulomb_mev*real(z,wp)**2*ar**(-1._wp/3._wp)
        terms%asymmetry_mev=model%a_asymmetry_mev*(ar-2._wp*real(z,wp))**2/ar
        terms%total_mev=terms%volume_mev+terms%surface_mev+terms%coulomb_mev+terms%asymmetry_mev
        if(.not.all(ieee_is_finite([terms%volume_mev,terms%surface_mev,terms%coulomb_mev,&
                                    terms%asymmetry_mev,terms%total_mev]))) then
            terms=liquid_drop_terms_t()
            call fail_status(status,NS_RANGE_ERROR,'numerics','liquid-drop term overflow');return
        end if
        terms%valid=.true.
        call clear_status(status)
    end subroutine

    subroutine evaluate_liquid_drop_candidate(model,n_b_m3,a,z,state,status)
        type(liquid_drop_model_t),intent(in)::model
        real(wp),intent(in)::n_b_m3
        integer,intent(in)::a,z
        type(liquid_drop_state_t),intent(out)::state
        type(status_t),intent(out)::status
        type(liquid_drop_terms_t)::terms
        type(electron_state_t)::electron
        real(wp)::ar,n_fm3,neutrality_scale
        state=liquid_drop_state_t()
        if(.not.model%ready) then
            call fail_status(status,NS_UNINITIALIZED,'implementation','liquid-drop model not ready');return
        end if
        if(.not.ieee_is_finite(n_b_m3).or.n_b_m3<=0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','baryon density must be finite and positive');return
        end if
        call liquid_drop_nuclear_terms(model,a,z,terms,status)
        if(.not.is_ok(status))return
        n_fm3=n_b_m3/inverse_fm3_in_si
        if(n_fm3<model%n_b_min_fm3.or.n_fm3>model%n_b_max_fm3) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics',&
                'T12 toy is admitted only for the declared fully-ionized pre-neutron-drip density window');return
        end if
        ar=real(a,wp)
        state%a=a;state%z=z;state%n_b_m3=n_b_m3;state%n_b_fm3=n_fm3
        state%n_nucleus_m3=n_b_m3/ar
        state%n_e_m3=state%n_nucleus_m3*real(z,wp)
        call electron_state(state%n_e_m3,electron,status)
        if(.not.is_ok(status)) then
            state=liquid_drop_state_t();return
        end if
        state%volume_mev_per_baryon=terms%volume_mev/ar
        state%surface_mev_per_baryon=terms%surface_mev/ar
        state%coulomb_mev_per_baryon=terms%coulomb_mev/ar
        state%asymmetry_mev_per_baryon=terms%asymmetry_mev/ar
        state%nuclear_mev_per_baryon=terms%total_mev/ar
        state%electron_mev_per_baryon=(electron%energy_j_m3/n_b_m3)/mev_in_j
        state%total_mev_per_baryon=state%nuclear_mev_per_baryon+state%electron_mev_per_baryon
        neutrality_scale=max(1._wp,abs(state%n_e_m3),abs(real(z,wp)*state%n_nucleus_m3))
        state%charge_neutrality_closed=abs(state%n_e_m3-real(z,wp)*state%n_nucleus_m3)<=&
            16._wp*epsilon(1._wp)*neutrality_scale
        if(.not.all(ieee_is_finite([state%volume_mev_per_baryon,state%surface_mev_per_baryon,&
                                    state%coulomb_mev_per_baryon,state%asymmetry_mev_per_baryon,&
                                    state%nuclear_mev_per_baryon,state%electron_mev_per_baryon,&
                                    state%total_mev_per_baryon]))) then
            state=liquid_drop_state_t()
            call fail_status(status,NS_RANGE_ERROR,'numerics','liquid-drop candidate is not representable');return
        end if
        state%valid=.true.
        call clear_status(status)
    end subroutine

    subroutine scan_liquid_drop_candidates(model,n_b_m3,a_min,a_max,a_step,z_min,z_max,z_step,result,status)
        type(liquid_drop_model_t),intent(in)::model
        real(wp),intent(in)::n_b_m3
        integer,intent(in)::a_min,a_max,a_step,z_min,z_max,z_step
        type(liquid_drop_scan_t),intent(out)::result
        type(status_t),intent(out)::status
        type(liquid_drop_state_t)::trial
        integer::a,z,z_hi,last_a,last_z
        real(wp)::best_energy
        result=liquid_drop_scan_t();best_energy=huge(1._wp)
        if(a_min<=0.or.a_max<a_min.or.a_step<=0.or.z_min<0.or.z_max<z_min.or.z_step<=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid discrete liquid-drop candidate mesh');return
        end if
        ! Establish the density-domain status before scanning the candidate mesh.
        call evaluate_liquid_drop_candidate(model,n_b_m3,a_min,min(max(z_min,0),a_min),trial,status)
        if(.not.is_ok(status))return
        do a=a_min,a_max,a_step
            z_hi=min(z_max,a)
            if(z_min>z_hi)cycle
            do z=z_min,z_hi,z_step
                call evaluate_liquid_drop_candidate(model,n_b_m3,a,z,trial,status)
                if(.not.is_ok(status))return
                result%candidate_count=result%candidate_count+1
                if(trial%total_mev_per_baryon<best_energy) then
                    best_energy=trial%total_mev_per_baryon
                    result%best=trial;result%best_a=a;result%best_z=z
                end if
            end do
        end do
        if(result%candidate_count==0) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','candidate mesh contains no admissible A/Z pair');return
        end if
        last_a=a_min+((a_max-a_min)/a_step)*a_step
        z_hi=min(z_max,result%best_a)
        last_z=z_min+((z_hi-z_min)/z_step)*z_step
        result%minimum_on_boundary=result%best_a==a_min.or.result%best_a==last_a.or.&
            result%best_z==z_min.or.result%best_z==last_z
        result%stationary_point_claimed=.false.
        result%mesh_is_discrete=.true.
        result%valid=.true.
        call clear_status(status)
    end subroutine

end module ns_liquid_drop
