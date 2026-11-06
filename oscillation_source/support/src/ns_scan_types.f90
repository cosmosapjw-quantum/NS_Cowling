! Shared records only. Scan/shooting algorithms belong to each workspace.
module ns_scan_types
    use ns_precision, only: wp
    use ns_status, only: status_t
    use ns_numerics_types, only: ode_options_t,root_options_t,root_result_t
    implicit none
    private
    integer,parameter,public :: CAND_SIGN=1,CAND_NEAR_BRACKET=2,CAND_SAMPLE=3
    type,public :: trial_record_t
        real(wp) :: parameter=0._wp,residual=0._wp
        integer :: orientation=1
        logical :: valid=.false.
        type(status_t) :: status
    end type
    type,public :: bracket_candidate_t
        real(wp) :: lo=0._wp,hi=0._wp
        integer :: left_index=0,right_index=0,kind=CAND_SAMPLE
    end type
    type,public :: scan_result_t
        type(trial_record_t),allocatable :: trials(:)
        type(bracket_candidate_t),allocatable :: candidates(:)
        integer :: barriers=0,sign_brackets=0,sample_candidates=0
        logical :: valid=.false.,all_trials_succeeded=.false.
        type(status_t) :: status
    end type
    type,public :: shooting_options_t
        real(wp) :: amplitude=1._wp,zero_tol=1e-10_wp,residual_tol=1e-8_wp
        character(len=16) :: method='adaptive'
        integer :: fixed_steps=256
        type(ode_options_t) :: ode
        type(root_options_t) :: root
    end type
    type,public :: string_shot_t
        real(wp) :: k_per_m=0._wp,residual=0._wp,boundary_y_m=0._wp,boundary_slope=0._wp
        real(wp),allocatable :: x_m(:),y_m(:),slope(:)
        integer :: accepted=0,rejected=0,evaluations=0
        logical :: valid=.false.
        type(status_t) :: status
    end type
    type,public :: string_root_t
        type(root_result_t) :: search
        type(string_shot_t) :: shot
        integer :: candidate_index=0
        logical :: valid=.false.
        type(status_t) :: status
    end type
    type,public :: string_spectrum_t
        type(scan_result_t) :: scan
        type(string_root_t),allocatable :: roots(:),refinements(:)
        integer :: roots_found=0,unresolved_candidates=0
        logical :: valid=.false.
        type(status_t) :: status
    end type
end module ns_scan_types
