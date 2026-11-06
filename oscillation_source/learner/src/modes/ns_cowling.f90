module ns_cowling
    use ns_precision, only: wp
    use ns_status
    use ns_numerics_types, only: ode_options_t
    use ns_background_access
    implicit none
    private
    type,public :: cowling_options_t
        integer :: ell=2
        real(wp) :: amplitude=1._wp,nu_shift=0._wp
        type(ode_options_t) :: ode
    end type
    type,public :: cowling_shot_t
        real(wp),allocatable :: q(:),r_m(:),w_m3(:),v_m2(:)
        real(wp) :: omega_rad_s=0._wp,residual=0._wp,q_over_r2=0._wp
        integer :: accepted=0,rejected=0,evaluations=0
        logical :: valid=.false.
        type(status_t) :: status
    end type
    public :: cowling_rhs_point,shoot_cowling
contains
    subroutine cowling_rhs_point(r,omega,chi,nu,lambda,dnu,y,dydr,qterm,status)
        real(wp),intent(in)::r,omega,chi,nu,lambda,dnu,y(:)
        real(wp),intent(out)::dydr(:),qterm
        type(status_t),intent(out)::status
        dydr=0._wp;qterm=0._wp
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','T07 learner Cowling RHS not implemented')
    end subroutine
    subroutine shoot_cowling(background,omega,options,result)
        type(background_table_t),intent(in)::background
        real(wp),intent(in)::omega
        type(cowling_options_t),intent(in)::options
        type(cowling_shot_t),intent(out)::result
        result=cowling_shot_t()
        call fail_status(result%status,NS_NOT_IMPLEMENTED,'implementation','T07 learner Cowling shot not implemented')
    end subroutine
end module ns_cowling
