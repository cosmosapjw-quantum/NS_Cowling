module ns_direct_join
    use ns_precision,only:wp
    use ns_status
    use ns_polytrope
    use ns_matching_types,only:direct_join_report_t
    implicit none
contains
    subroutine evaluate_direct_join(low,high,p,tol,report,status)
        type(polytrope_t),intent(in)::low,high
        real(wp),intent(in)::p,tol
        type(direct_join_report_t),intent(out)::report
        type(status_t),intent(out)::status
        report=direct_join_report_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student direct-join slot not implemented')
    end subroutine
end module
