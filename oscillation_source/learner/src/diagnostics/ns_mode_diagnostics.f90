! Ch10 learner slot. Difference/resolution logic is not supplied by the solution package at runtime.
module ns_mode_diagnostics
    use ns_precision, only: wp
    use ns_status
    implicit none
    private
    type,public :: mode_difference_t
        real(wp)::delta_f_hz=0._wp,delta_p1_hz=0._wp,delta_f_relative=0._wp,delta_p1_relative=0._wp
        character(len=24)::f_interpretation='unresolved',p1_interpretation='unresolved'
        logical::f_relative_defined=.false.,p1_relative_defined=.false.,identity_compared=.false.,valid=.false.
    end type mode_difference_t
    type,public :: paired_refinement_t
        real(wp)::f_effect_hz=0._wp,p1_effect_hz=0._wp
        real(wp)::f_numerical_floor_hz=0._wp,p1_numerical_floor_hz=0._wp
        real(wp)::f_mass_floor_hz=0._wp,p1_mass_floor_hz=0._wp
        real(wp)::f_resolution_floor_hz=0._wp,p1_resolution_floor_hz=0._wp
        character(len=24)::f_interpretation='unresolved',p1_interpretation='unresolved'
        logical::f_convergent=.false.,p1_convergent=.false.,mode_identity_consistent=.false.,valid=.false.
    end type paired_refinement_t
    public :: compare_mode_frequencies,classify_paired_refinement
contains
    subroutine compare_mode_frequencies(f_ref,f_cmp,p1_ref,p1_cmp,resolution,result,status)
        real(wp),intent(in)::f_ref,f_cmp,p1_ref,p1_cmp,resolution
        type(mode_difference_t),intent(out)::result
        type(status_t),intent(out)::status
        result=mode_difference_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','Ch10 learner mode-difference diagnostic not implemented')
    end subroutine compare_mode_frequencies

    subroutine classify_paired_refinement(df0,df1,df2,dp0,dp1,dp2,absolute_resolution, &
        mass_residual_ref,mass_residual_cmp,df_dmass,dp_dmass,f_ref_label,f_cmp_label, &
        p1_ref_label,p1_cmp_label,result,status)
        real(wp),intent(in)::df0,df1,df2,dp0,dp1,dp2,absolute_resolution
        real(wp),intent(in)::mass_residual_ref,mass_residual_cmp,df_dmass,dp_dmass
        character(len=*),intent(in)::f_ref_label,f_cmp_label,p1_ref_label,p1_cmp_label
        type(paired_refinement_t),intent(out)::result
        type(status_t),intent(out)::status
        result=paired_refinement_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','Ch10 learner paired refinement not implemented')
    end subroutine classify_paired_refinement
end module ns_mode_diagnostics
