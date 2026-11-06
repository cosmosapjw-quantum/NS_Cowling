! Shared matching records only. Physics algorithms stay in the selected workspace.
module ns_matching_types
    use ns_precision, only: wp
    implicit none
    private
    type, public :: direct_join_report_t
        real(wp) :: pressure_pa=0._wp, pressure_residual_pa=0._wp
        real(wp) :: mu_low_j=0._wp, mu_high_j=0._wp, mu_mismatch_j=0._wp
        real(wp) :: n_low_m3=0._wp, n_high_m3=0._wp
        real(wp) :: energy_low_j_m3=0._wp, energy_high_j_m3=0._wp
        logical :: pressure_closed=.false., chemical_closed=.false.
        logical :: thermodynamically_consistent=.false., ready=.false.
    end type
    type, public :: transition_state_t
        real(wp) :: mu_t_j=0._wp, p_t_pa=0._wp
        real(wp) :: p_low_pa=0._wp, p_high_pa=0._wp, pressure_residual_pa=0._wp
        real(wp) :: mu_low_j=0._wp, mu_high_j=0._wp, mu_mismatch_j=0._wp
        real(wp) :: n_low_m3=0._wp, n_high_m3=0._wp, delta_n_m3=0._wp
        real(wp) :: energy_low_j_m3=0._wp, energy_high_j_m3=0._wp, delta_e_j_m3=0._wp
        integer :: root_iterations=0, root_evaluations=0
        logical :: branch_selection_valid=.false.
        logical :: thermodynamically_consistent=.false., valid=.false.
    end type
    type, public :: mu_bridge_t
        real(wp) :: mu_low_j=0._wp, mu_high_j=0._wp
        real(wp) :: p_low_pa=0._wp, p_high_pa=0._wp
        real(wp) :: n_low_m3=0._wp, n_high_m3=0._wp
        real(wp) :: a0=0._wp,a1=0._wp,a2=0._wp,a3=0._wp
        real(wp) :: min_d2p=0._wp, min_causal_margin=0._wp
        logical :: interval_extrema_checked=.false., stable_interval=.false., causal_interval=.false.
        logical :: ready=.false.
    end type
end module ns_matching_types
