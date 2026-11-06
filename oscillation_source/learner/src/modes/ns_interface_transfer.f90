module ns_interface_transfer
    use ns_precision, only: wp
    use ns_status
    implicit none
    private
    public :: material_interface_transfer,lagrangian_pressure_jump,identify_interface_family
contains
    subroutine lagrangian_pressure_jump(r,omega,p,eps,nu,lambda,dnu,w,v,delta_p,qterm,status)
        real(wp),intent(in)::r,omega,p,eps,nu,lambda,dnu,w,v
        real(wp),intent(out)::delta_p,qterm
        type(status_t),intent(out)::status
        delta_p=0._wp;qterm=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','learner interface transfer is not implemented')
    end subroutine
    subroutine material_interface_transfer(r,omega,p,eps_in,eps_out,nu,lambda,dnu,w_in,v_in, &
        inner_side,outer_side,w_out,v_out,dp_in,dp_out,status)
        real(wp),intent(in)::r,omega,p,eps_in,eps_out,nu,lambda,dnu,w_in,v_in
        integer,intent(in)::inner_side,outer_side
        real(wp),intent(out)::w_out,v_out,dp_in,dp_out
        type(status_t),intent(out)::status
        w_out=0._wp;v_out=0._wp;dp_in=0._wp;dp_out=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','learner interface transfer is not implemented')
    end subroutine
    subroutine identify_interface_family(reference_label,root_valid,transfer_valid,overlap,candidate_interface_g,label,status)
        character(len=*),intent(in)::reference_label
        logical,intent(in)::root_valid,transfer_valid,candidate_interface_g
        real(wp),intent(in)::overlap
        character(len=*),intent(out)::label
        type(status_t),intent(out)::status
        label='unresolved'
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','learner interface identity is not implemented')
    end subroutine
end module ns_interface_transfer
