! Learning slot: implement here; never forwards to the solutions package.
module ns_bisection
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: scalar_callback,root_options_t,root_result_t
    implicit none
contains
    recursive subroutine bisect(f,a,b,options,result)
        procedure(scalar_callback)::f
        real(wp),intent(in)::a,b
        type(root_options_t),intent(in)::options
        type(root_result_t),intent(out)::result
        result=root_result_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','implement bisection: lab ch02')
    end subroutine
end module
