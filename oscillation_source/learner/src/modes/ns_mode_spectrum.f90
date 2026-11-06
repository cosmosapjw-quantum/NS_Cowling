! Root convergence and smooth-mode identity are independent output states.
module ns_mode_spectrum
    use ns_precision, only: wp
    use ns_physical_constants, only: two_pi
    use ns_status
    use ns_scan_types
    use ns_numerics_types, only: root_options_t,root_result_t
    use ns_background_access
    use ns_cowling
    use ns_mode_normalization
    use ns_mode_identity
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    implicit none
    private
    type,public :: mode_options_t
        type(cowling_options_t) :: shot
        type(root_options_t) :: root
        character(len=24) :: method='adaptive'
        integer :: cell_substeps=8
        real(wp)::residual_tol=1e-5_wp,zero_tol=1e-10_wp,node_zero=1e-7_wp
        real(wp)::node_xlo=0.001_wp,node_xhi=0.999_wp
    end type
    type,public :: mode_result_t
        type(root_result_t)::search
        type(cowling_shot_t)::shot
        real(wp),allocatable::x(:),wn(:),vn(:)
        integer::w_nodes=-1,v_nodes=-1
        character(len=16)::label='unresolved'
        type(status_t)::root_status,identity_status
        logical::root_valid=.false.,diagnostics_valid=.false.
    end type
    type,public :: mode_spectrum_t
        type(scan_result_t)::scan
        type(mode_result_t),allocatable::refinements(:)
        integer::roots_found=0,f_index=0,p1_index=0
        logical::valid=.false.
        type(status_t)::status
    end type
    public::refine_mode,scan_modes,evaluate_mode_shot
contains
    subroutine evaluate_mode_shot(bg,omega,options,shot)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::omega
        type(mode_options_t),intent(in)::options
        type(cowling_shot_t),intent(out)::shot
        shot=cowling_shot_t()
        call fail_status(shot%status,NS_NOT_IMPLEMENTED,'implementation','T08 learner shot selection not implemented')
    end subroutine
    subroutine refine_mode(bg,candidate,options,result)
        type(background_table_t),intent(in)::bg
        type(bracket_candidate_t),intent(in)::candidate
        type(mode_options_t),intent(in)::options
        type(mode_result_t),intent(out)::result
        result=mode_result_t()
        call fail_status(result%root_status,NS_NOT_IMPLEMENTED,'implementation','T08 learner mode refinement not implemented')
        call fail_status(result%identity_status,NS_IDENTITY_UNRESOLVED,'physics','no computed root')
    end subroutine
    subroutine scan_modes(bg,grid,options,result)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::grid(:)
        type(mode_options_t),intent(in)::options
        type(mode_spectrum_t),intent(out)::result
        result=mode_spectrum_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','T08 learner mode scan not implemented')
    end subroutine
end module ns_mode_spectrum
