module ns_direct_join
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use ns_polytrope
    use ns_matching_types, only: direct_join_report_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: evaluate_direct_join
contains
    pure subroutine evaluate_direct_join(low,high,pressure_pa,mu_tolerance_j,report,status)
        type(polytrope_t),intent(in)::low,high
        real(wp),intent(in)::pressure_pa,mu_tolerance_j
        type(direct_join_report_t),intent(out)::report
        type(status_t),intent(out)::status
        type(eos_state_t)::a,b
        report=direct_join_report_t()
        if(.not.ieee_is_finite(pressure_pa) .or. .not.ieee_is_finite(mu_tolerance_j)) then
            call fail_status(status,NS_NONFINITE_INPUT,'physics','nonfinite direct-join input');return
        end if
        if(pressure_pa<=0._wp .or. mu_tolerance_j<0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','direct join requires P>0 and tolerance>=0');return
        end if
        call eos_at_pressure(low,pressure_pa,SIDE_REGULAR,a,status);if(.not.is_ok(status))return
        call eos_at_pressure(high,pressure_pa,SIDE_REGULAR,b,status);if(.not.is_ok(status))return
        report%pressure_pa=pressure_pa;report%pressure_residual_pa=b%pressure_pa-a%pressure_pa
        report%mu_low_j=a%mu_j;report%mu_high_j=b%mu_j;report%mu_mismatch_j=b%mu_j-a%mu_j
        report%n_low_m3=a%n_b_m3;report%n_high_m3=b%n_b_m3
        report%energy_low_j_m3=a%energy_j_m3;report%energy_high_j_m3=b%energy_j_m3
        report%pressure_closed=report%pressure_residual_pa==0._wp
        report%chemical_closed=abs(report%mu_mismatch_j)<=mu_tolerance_j
        report%thermodynamically_consistent=report%pressure_closed .and. report%chemical_closed
        report%ready=.true.;call clear_status(status)
    end subroutine
end module ns_direct_join
