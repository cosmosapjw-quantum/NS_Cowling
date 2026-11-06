module ns_target_mass
    use ns_precision,only:wp
    use ns_status
    use ns_background_records
    use ns_tov,only:eos_pressure_callback
    implicit none
contains
    subroutine solve_target_mass(eos,m_target,pc_bracket,star_options,options,result,events)
        procedure(eos_pressure_callback)::eos
        real(wp),intent(in)::m_target,pc_bracket(2)
        type(star_options_t),intent(in)::star_options
        type(mass_options_t),intent(in)::options
        type(star_result_t),intent(out)::result
        real(wp),optional,intent(in)::events(:)
        result=star_result_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','learner target-mass residual')
    end subroutine
end module
