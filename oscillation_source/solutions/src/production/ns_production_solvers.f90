! Public production dispatch for smooth and sharp-Maxwell stellar backgrounds.
!
! The numerical kernels remain separately testable. These generic entry points
! select the physical representation from the procedure arguments / background
! type without interpolating across a density discontinuity.
module ns_production_solvers
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types, only: eos_state_t
    use ns_background_records, only: star_options_t,star_profile_t
    use ns_background_access, only: background_table_t
    use ns_tov, only: eos_pressure_callback,integrate_tov
    use ns_maxwell_tov, only: maxwell_star_t,integrate_maxwell_tov
    use ns_mode_spectrum, only: mode_options_t,mode_spectrum_t,scan_modes
    use ns_maxwell_mode_spectrum, only: maxwell_mode_spectrum_t,scan_maxwell_modes
    implicit none
    private

    interface solve_tov_star
        module procedure solve_tov_star_smooth
        module procedure solve_tov_star_maxwell
    end interface

    interface solve_cowling_spectrum
        module procedure solve_cowling_spectrum_background
        module procedure solve_cowling_spectrum_profile
        module procedure solve_cowling_spectrum_maxwell
    end interface

    public :: solve_tov_star,solve_cowling_spectrum
contains
    subroutine solve_tov_star_smooth(eos,pc,options,profile,eos_events)
        procedure(eos_pressure_callback) :: eos
        real(wp),intent(in) :: pc
        type(star_options_t),intent(in) :: options
        type(star_profile_t),intent(out) :: profile
        real(wp),intent(in),optional :: eos_events(:)
        if(present(eos_events)) then
            call integrate_tov(eos,pc,options,profile,eos_events)
        else
            call integrate_tov(eos,pc,options,profile)
        end if
    end subroutine

    subroutine solve_tov_star_maxwell(eos_high,eos_low,pc,pt,options,star)
        procedure(eos_pressure_callback) :: eos_high,eos_low
        real(wp),intent(in) :: pc,pt
        type(star_options_t),intent(in) :: options
        type(maxwell_star_t),intent(out) :: star
        call integrate_maxwell_tov(eos_high,eos_low,pc,pt,options,star)
    end subroutine

    subroutine solve_cowling_spectrum_background(background,grid,options,result)
        type(background_table_t),intent(in) :: background
        real(wp),intent(in) :: grid(:)
        type(mode_options_t),intent(in) :: options
        type(mode_spectrum_t),intent(out) :: result
        call scan_modes(background,grid,options,result)
    end subroutine

    subroutine solve_cowling_spectrum_profile(profile,grid,options,result)
        type(star_profile_t),intent(in) :: profile
        real(wp),intent(in) :: grid(:)
        type(mode_options_t),intent(in) :: options
        type(mode_spectrum_t),intent(out) :: result
        type(background_table_t) :: background
        type(status_t) :: status
        call profile_as_background(profile,background,status)
        if(.not.is_ok(status)) then
            result=mode_spectrum_t()
            result%status=status
            return
        end if
        call scan_modes(background,grid,options,result)
    end subroutine

    subroutine solve_cowling_spectrum_maxwell(star,grid,options,conversion,result)
        type(maxwell_star_t),intent(in) :: star
        real(wp),intent(in) :: grid(:)
        type(mode_options_t),intent(in) :: options
        character(len=*),intent(in) :: conversion
        type(maxwell_mode_spectrum_t),intent(out) :: result
        call scan_maxwell_modes(star,grid,conversion,options,result)
    end subroutine

    subroutine profile_as_background(profile,background,status)
        type(star_profile_t),intent(in) :: profile
        type(background_table_t),intent(out) :: background
        type(status_t),intent(out) :: status
        integer :: n,istat
        background=background_table_t()
        call fail_status(status,NS_INVALID_INPUT,'data','valid smooth TOV profile required for Cowling dispatch')
        if(.not.profile%valid) return
        if(.not.allocated(profile%q)) return
        if(.not.allocated(profile%r_m)) return
        if(.not.allocated(profile%m_kg)) return
        if(.not.allocated(profile%p_pa)) return
        if(.not.allocated(profile%energy_j_m3)) return
        if(.not.allocated(profile%n_b_m3)) return
        if(.not.allocated(profile%mu_j)) return
        if(.not.allocated(profile%chi)) return
        if(.not.allocated(profile%nu)) return
        if(.not.allocated(profile%lambda)) return
        if(.not.allocated(profile%dnu_dr)) return
        n=size(profile%q)
        if(n<2) return
        if(size(profile%r_m)/=n.or.size(profile%m_kg)/=n.or.size(profile%p_pa)/=n.or. &
           size(profile%energy_j_m3)/=n.or.size(profile%n_b_m3)/=n.or.size(profile%mu_j)/=n.or. &
           size(profile%chi)/=n.or.size(profile%nu)/=n.or.size(profile%lambda)/=n.or. &
           size(profile%dnu_dr)/=n) return
        allocate(background%q(n),background%r_m(n),background%m_kg(n),background%p_pa(n), &
            background%energy_j_m3(n),background%n_b_m3(n),background%mu_j(n),background%chi(n), &
            background%nu(n),background%lambda(n),background%dnu_dr(n),stat=istat)
        if(istat/=0) then
            call fail_status(status,NS_RANGE_ERROR,'environment','smooth Cowling background allocation failed')
            return
        end if
        background%q=profile%q
        background%r_m=profile%r_m
        background%m_kg=profile%m_kg
        background%p_pa=profile%p_pa
        background%energy_j_m3=profile%energy_j_m3
        background%n_b_m3=profile%n_b_m3
        background%mu_j=profile%mu_j
        background%chi=profile%chi
        background%nu=profile%nu
        background%lambda=profile%lambda
        background%dnu_dr=profile%dnu_dr
        background%n=n
        background%ready=.true.
        call clear_status(status)
    end subroutine
end module ns_production_solvers
