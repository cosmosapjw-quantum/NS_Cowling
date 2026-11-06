! Scoped to smooth barotropic ell=2 acoustic families. No frequency ranking.
module ns_mode_identity
    use ns_status
    implicit none
    private
    character(len=*),parameter,public :: MODE_IDENTITY_ID='smooth_l2_WV_nodes_v1'
    public :: identify_smooth_mode
contains
    subroutine identify_smooth_mode(root_valid,smooth,diagnostics_valid,w_nodes,v_nodes,label,status)
        logical,intent(in)::root_valid,smooth,diagnostics_valid
        integer,intent(in)::w_nodes,v_nodes
        character(len=*),intent(out)::label
        type(status_t),intent(out)::status
        label='unresolved'
        call fail_status(status,NS_IDENTITY_UNRESOLVED,'physics','requires a converged root and stable smooth W/V node evidence')
        if(.not.root_valid.or..not.smooth.or..not.diagnostics_valid) return
        if(w_nodes==0.and.v_nodes==0) then
            label='f'
        else if(w_nodes==1.and.v_nodes==1) then
            label='p1'
        else
            return
        end if
        call clear_status(status)
    end subroutine
end module ns_mode_identity
