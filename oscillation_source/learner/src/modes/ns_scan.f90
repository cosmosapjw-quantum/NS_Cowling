module ns_scan
    use ns_precision, only: wp
    use ns_status
    use ns_scan_types
    implicit none
    private
    public :: build_brackets,linear_scan_grid
contains
    subroutine linear_scan_grid(lo,hi,n,grid,status)
        real(wp),intent(in)::lo,hi
        integer,intent(in)::n
        real(wp),allocatable,intent(out)::grid(:)
        type(status_t),intent(out)::status
        ! Implement the declared strictly increasing finite parameter grid.
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T06 learner scan grid')
    end subroutine
    subroutine build_brackets(trials,zero_tol,result)
        type(trial_record_t),intent(in)::trials(:)
        real(wp),intent(in)::zero_tol
        type(scan_result_t),intent(out)::result
        result=scan_result_t()
        ! Failed or nonfinite trials must be barriers, not discarded samples.
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','T06 learner bracket builder')
    end subroutine
end module ns_scan
