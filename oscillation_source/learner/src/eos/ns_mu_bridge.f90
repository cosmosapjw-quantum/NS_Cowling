module ns_mu_bridge
    use ns_precision,only:wp
    use ns_status
    use ns_eos_types
    use ns_polytrope
    use ns_matching_types,only:mu_bridge_t
    implicit none
contains
    subroutine construct_branch_bridge(low,high,a,b,bridge,status)
        type(polytrope_t),intent(in)::low,high
        real(wp),intent(in)::a,b
        type(mu_bridge_t),intent(out)::bridge
        type(status_t),intent(out)::status
        bridge=mu_bridge_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student bridge slot not implemented')
    end subroutine
    subroutine construct_bridge_endpoints(a,p0,n0,b,p1,n1,bridge,status)
        real(wp),intent(in)::a,p0,n0,b,p1,n1
        type(mu_bridge_t),intent(out)::bridge
        type(status_t),intent(out)::status
        bridge=mu_bridge_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student bridge slot not implemented')
    end subroutine
    subroutine evaluate_mu_bridge(bridge,mu,state,dn,status)
        type(mu_bridge_t),intent(in)::bridge
        real(wp),intent(in)::mu
        type(eos_state_t),intent(out)::state
        real(wp),intent(out)::dn
        type(status_t),intent(out)::status
        call clear_eos_state(state);dn=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student bridge slot not implemented')
    end subroutine
end module
