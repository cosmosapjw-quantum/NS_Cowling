module ns_maxwell
    use ns_precision,only:wp
    use ns_status
    use ns_polytrope
    use ns_matching_types,only:transition_state_t
    use ns_numerics_types,only:root_options_t
    implicit none
contains
    subroutine solve_maxwell_match(low,high,a,b,options,transition,status)
        type(polytrope_t),intent(in)::low,high
        real(wp),intent(in)::a,b
        type(root_options_t),intent(in)::options
        type(transition_state_t),intent(out)::transition
        type(status_t),intent(out)::status
        transition=transition_state_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student Maxwell slot not implemented')
    end subroutine
end module
