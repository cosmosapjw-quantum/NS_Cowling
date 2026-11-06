! Shared state/metadata only. No TOV solver is provided to learners here.
module ns_background_records
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: ode_options_t,root_options_t,root_result_t
    implicit none
    private
    type,public :: star_options_t
        real(wp) :: r0_m=1._wp,p_cut_pa=1e22_wp,r_limit_m=1e6_wp,nu_start=0._wp
        type(ode_options_t) :: ode
        character(len=80) :: eos_id='unspecified'
        character(len=64) :: surface_policy='finite_cutoff_schwarzschild_time'
    end type
    type,public :: star_profile_t
        real(wp),allocatable :: q(:),r_m(:),m_kg(:),p_pa(:),energy_j_m3(:),n_b_m3(:),mu_j(:),chi(:)
        real(wp),allocatable :: nu(:),lambda(:),dnu_dr(:)
        integer,allocatable :: segment(:),side(:)
        real(wp) :: pc_pa=0._wp,p_start_pa=0._wp,r0_m=0._wp,p_cut_pa=0._wp
        real(wp) :: radius_m=0._wp,mass_kg=0._wp,surface_residual_pa=0._wp,lapse_residual=0._wp
        real(wp) :: lapse_shift=0._wp
        integer :: accepted=0,rejected=0,evaluations=0,stops_hit=0
        logical :: valid=.false.
        character(len=80) :: eos_id='unspecified'
        character(len=64) :: surface_policy='',coordinate='log_pressure',branch_id='single_pc'
        type(status_t) :: status
    end type
    type,public :: mass_options_t
        integer :: scan_n=9
        real(wp) :: mass_atol_kg=1e20_wp,mass_rtol=1e-8_wp
        type(root_options_t) :: root
        character(len=64) :: branch_id='sampled_increasing_pc'
    end type
    type,public :: star_result_t
        type(star_profile_t) :: profile
        type(root_result_t) :: root
        real(wp),allocatable :: scan_pc_pa(:),scan_mass_kg(:)
        real(wp) :: mass_residual_kg=0._wp,target_mass_kg=0._wp
        integer :: star_evaluations=0
        logical :: valid=.false.,scan_monotone=.false.
        type(status_t) :: status
    end type
end module
