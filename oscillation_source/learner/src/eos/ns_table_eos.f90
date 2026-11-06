! Ch09 learner slot. Same public table-EOS contract as solutions, no fallback.
module ns_table_eos
    use ns_precision, only: wp
    use ns_status
    use ns_eos_types
    use ns_matching_types, only: mu_bridge_t
    use ns_table_io, only: normalized_table_t
    implicit none
    private
    character(len=*), parameter, public :: TABLE_INTERPOLATION_ID='learner_not_implemented'
    type, public :: table_audit_t
        integer :: row_count=0,interval_count=0,admitted_intervals=0,rejected_intervals=0
        integer :: first_rejected_interval=0
        logical :: mu_strict_increasing=.false.,pressure_strict_increasing=.false.
        logical :: fully_admitted=.false.
        real(wp) :: raw_energy_closure_max_rel=0._wp
        real(wp) :: raw_trapezoid_first_law_max_rel=0._wp
        real(wp) :: first_bad_min_d2p=0._wp,first_bad_min_causal_margin=0._wp
    end type table_audit_t
    type, public :: table_eos_t
        type(normalized_table_t) :: raw
        type(mu_bridge_t), allocatable :: cells(:)
        logical :: ready=.false.
    end type table_eos_t
    public :: audit_table_eos,init_table_eos,eos_table_at_mu,eos_table_at_pressure
contains
    subroutine audit_table_eos(raw,audit,status)
        type(normalized_table_t),intent(in)::raw
        type(table_audit_t),intent(out)::audit
        type(status_t),intent(out)::status
        audit=table_audit_t()
        audit%row_count=raw%n
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','Ch09 learner table admission not implemented')
    end subroutine audit_table_eos

    subroutine init_table_eos(raw,model,audit,status)
        type(normalized_table_t),intent(in)::raw
        type(table_eos_t),intent(out)::model
        type(table_audit_t),intent(out)::audit
        type(status_t),intent(out)::status
        model=table_eos_t();audit=table_audit_t();model%raw=raw
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','Ch09 learner table representation not implemented')
    end subroutine init_table_eos

    subroutine eos_table_at_mu(model,mu_j,side,state,dn_dmu,status)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::mu_j
        integer,intent(in)::side
        type(eos_state_t),intent(out)::state
        real(wp),intent(out)::dn_dmu
        type(status_t),intent(out)::status
        call clear_eos_state(state);dn_dmu=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','Ch09 learner P(mu) evaluation not implemented')
    end subroutine eos_table_at_mu

    subroutine eos_table_at_pressure(model,p_pa,side,state,status)
        type(table_eos_t),intent(in)::model
        real(wp),intent(in)::p_pa
        integer,intent(in)::side
        type(eos_state_t),intent(out)::state
        type(status_t),intent(out)::status
        call clear_eos_state(state)
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','Ch09 learner pressure inverse not implemented')
    end subroutine eos_table_at_pressure
end module ns_table_eos
