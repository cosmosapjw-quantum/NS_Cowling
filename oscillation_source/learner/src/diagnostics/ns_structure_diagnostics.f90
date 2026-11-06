! Ch10 learner slot. Read-only diagnostics must be implemented by the learner.
module ns_structure_diagnostics
    use ns_precision, only: wp
    use ns_status
    use ns_background_access, only: background_table_t
    implicit none
    private
    character(len=*),parameter,public :: STRUCTURE_DIAGNOSTIC_ID='learner_not_implemented'
    type,public :: structure_diagnostics_t
        real(wp) :: radius_m=0._wp,mass_kg=0._wp,proper_radius_m=0._wp,acoustic_time_s=0._wp
        real(wp) :: r_split_m=0._wp,outer_areal_thickness_m=0._wp,outer_proper_thickness_m=0._wp
        real(wp) :: outer_acoustic_time_s=0._wp,sound_speed_min_m_s=0._wp,sound_speed_max_m_s=0._wp
        real(wp) :: f_hat=0._wp,p1_hat=0._wp
        logical :: valid=.false.
    end type structure_diagnostics_t
    public :: compute_structure_diagnostics
contains
    subroutine compute_structure_diagnostics(bg,split_pressure_pa,f_hz,p1_hz,result,status)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::split_pressure_pa,f_hz,p1_hz
        type(structure_diagnostics_t),intent(out)::result
        type(status_t),intent(out)::status
        result=structure_diagnostics_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','Ch10 learner structure diagnostics not implemented')
    end subroutine compute_structure_diagnostics
end module ns_structure_diagnostics
