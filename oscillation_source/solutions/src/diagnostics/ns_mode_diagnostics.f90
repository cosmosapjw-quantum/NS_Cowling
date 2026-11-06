! T15 paired-frequency diagnostics. Mode identity is deliberately not inferred here.
module ns_mode_diagnostics
    use ns_precision, only: wp
    use ns_status
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
        result=mode_difference_t();call fail_status(status,NS_INVALID_INPUT,'data','valid paired frequencies required')
        if(.not.all(ieee_is_finite([f_ref,f_cmp,p1_ref,p1_cmp,resolution]))) return
        if(f_ref<=0._wp.or.f_cmp<=0._wp.or.p1_ref<=0._wp.or.p1_cmp<=0._wp.or.resolution<=0._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','positive frequencies and absolute resolution required');return
        end if
        result%delta_f_hz=f_cmp-f_ref;result%delta_p1_hz=p1_cmp-p1_ref
        result%f_relative_defined=.true.;result%p1_relative_defined=.true.
        result%delta_f_relative=result%delta_f_hz/f_ref;result%delta_p1_relative=result%delta_p1_hz/p1_ref
        if(abs(result%delta_f_hz)<=resolution) then
            result%f_interpretation='below_resolution'
        else
            result%f_interpretation='resolved'
        end if
        if(abs(result%delta_p1_hz)<=resolution) then
            result%p1_interpretation='below_resolution'
        else
            result%p1_interpretation='resolved'
        end if
        result%identity_compared=.false.;result%valid=.true.;call clear_status(status)
    end subroutine compare_mode_frequencies

    subroutine classify_paired_refinement(df0,df1,df2,dp0,dp1,dp2,absolute_resolution, &
        mass_residual_ref,mass_residual_cmp,df_dmass,dp_dmass,f_ref_label,f_cmp_label, &
        p1_ref_label,p1_cmp_label,result,status)
        real(wp),intent(in)::df0,df1,df2,dp0,dp1,dp2,absolute_resolution
        real(wp),intent(in)::mass_residual_ref,mass_residual_cmp,df_dmass,dp_dmass
        character(len=*),intent(in)::f_ref_label,f_cmp_label,p1_ref_label,p1_cmp_label
        type(paired_refinement_t),intent(out)::result
        type(status_t),intent(out)::status
        real(wp)::mass_floor_input,f_previous_change,f_latest_change,p1_previous_change,p1_latest_change
        result=paired_refinement_t()
        call fail_status(status,NS_INVALID_INPUT,'data','valid paired refinement data required')
        if(.not.all(ieee_is_finite([df0,df1,df2,dp0,dp1,dp2,absolute_resolution, &
            mass_residual_ref,mass_residual_cmp,df_dmass,dp_dmass]))) return
        if(absolute_resolution<=0._wp.or.mass_residual_ref<0._wp.or.mass_residual_cmp<0._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'numerics','nonnegative mass residuals and positive resolution required')
            return
        end if
        result%mode_identity_consistent=trim(f_ref_label)=='f'.and.trim(f_cmp_label)=='f'.and. &
            trim(p1_ref_label)=='p1'.and.trim(p1_cmp_label)=='p1'
        if(.not.result%mode_identity_consistent) then
            call fail_status(status,NS_IDENTITY_UNRESOLVED,'physics','paired comparison mode labels are inconsistent')
            return
        end if
        result%f_effect_hz=abs(df2);result%p1_effect_hz=abs(dp2)
        f_previous_change=abs(df1-df0);f_latest_change=abs(df2-df1)
        p1_previous_change=abs(dp1-dp0);p1_latest_change=abs(dp2-dp1)
        result%f_numerical_floor_hz=max(f_previous_change,f_latest_change)
        result%p1_numerical_floor_hz=max(p1_previous_change,p1_latest_change)
        result%f_convergent=f_latest_change<=f_previous_change
        result%p1_convergent=p1_latest_change<=p1_previous_change
        mass_floor_input=max(mass_residual_ref,mass_residual_cmp)
        result%f_mass_floor_hz=abs(df_dmass)*mass_floor_input
        result%p1_mass_floor_hz=abs(dp_dmass)*mass_floor_input
        result%f_resolution_floor_hz=max(absolute_resolution,result%f_numerical_floor_hz,result%f_mass_floor_hz)
        result%p1_resolution_floor_hz=max(absolute_resolution,result%p1_numerical_floor_hz,result%p1_mass_floor_hz)
        if(.not.result%f_convergent) then
            result%f_interpretation='convergence_unresolved'
        else if(result%f_effect_hz>result%f_resolution_floor_hz) then
            result%f_interpretation='resolved'
        else
            result%f_interpretation='below_resolution'
        end if
        if(.not.result%p1_convergent) then
            result%p1_interpretation='convergence_unresolved'
        else if(result%p1_effect_hz>result%p1_resolution_floor_hz) then
            result%p1_interpretation='resolved'
        else
            result%p1_interpretation='below_resolution'
        end if
        result%valid=.true.;call clear_status(status)
    end subroutine classify_paired_refinement
end module ns_mode_diagnostics
