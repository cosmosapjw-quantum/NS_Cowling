module ns_shooting_toy
    use ns_precision, only: wp
    use ns_status
    use ns_scan_types
    implicit none
    private
    public :: shoot_string,refine_string,scan_string
contains
    subroutine shoot_string(length_m,k_per_m,options,result,fault)
        real(wp),intent(in)::length_m,k_per_m
        type(shooting_options_t),intent(in)::options
        type(string_shot_t),intent(out)::result
        character(len=*),intent(in),optional::fault
        result=string_shot_t()
        ! Use this workspace's IVP solver, not the analytic solution or solutions/.
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','T06 learner IVP shooting')
    end subroutine
    subroutine refine_string(length_m,candidate,options,result)
        real(wp),intent(in)::length_m
        type(bracket_candidate_t),intent(in)::candidate
        type(shooting_options_t),intent(in)::options
        type(string_root_t),intent(out)::result
        result=string_root_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','T06 learner root refinement')
    end subroutine
    subroutine scan_string(length_m,grid,options,result,failed_index)
        real(wp),intent(in)::length_m,grid(:)
        type(shooting_options_t),intent(in)::options
        type(string_spectrum_t),intent(out)::result
        integer,intent(in),optional::failed_index
        result=string_spectrum_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','T06 learner spectrum scan')
    end subroutine
end module ns_shooting_toy
