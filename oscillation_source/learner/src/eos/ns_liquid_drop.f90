module ns_liquid_drop
    use ns_precision,only:wp
    use ns_status
    implicit none
    private
    type,public :: liquid_drop_model_t
        logical::ready=.false.
        real(wp)::a_volume_mev=15.75_wp,a_surface_mev=17.8_wp,a_coulomb_mev=0.711_wp,a_asymmetry_mev=23.7_wp
        real(wp)::n_b_min_fm3=1e-10_wp,n_b_max_fm3=2.5e-4_wp
    end type
    type,public :: liquid_drop_terms_t
        logical::valid=.false.
        real(wp)::volume_mev=0._wp,surface_mev=0._wp,coulomb_mev=0._wp
        real(wp)::asymmetry_mev=0._wp,total_mev=0._wp
    end type
    type,public :: liquid_drop_state_t
        logical::valid=.false.,charge_neutrality_closed=.false.;integer::a=0,z=0
        real(wp)::n_b_m3=0._wp,n_b_fm3=0._wp,n_nucleus_m3=0._wp,n_e_m3=0._wp
        real(wp)::volume_mev_per_baryon=0._wp,surface_mev_per_baryon=0._wp,coulomb_mev_per_baryon=0._wp
        real(wp)::asymmetry_mev_per_baryon=0._wp,nuclear_mev_per_baryon=0._wp
        real(wp)::electron_mev_per_baryon=0._wp,total_mev_per_baryon=0._wp
    end type
    type,public :: liquid_drop_scan_t
        logical::valid=.false.,minimum_on_boundary=.false.,stationary_point_claimed=.false.,mesh_is_discrete=.true.
        integer::candidate_count=0,best_a=0,best_z=0;type(liquid_drop_state_t)::best
    end type
    public :: init_liquid_drop_model,liquid_drop_nuclear_terms,evaluate_liquid_drop_candidate,scan_liquid_drop_candidates
contains
    subroutine init_liquid_drop_model(model,status)
        type(liquid_drop_model_t),intent(out)::model;type(status_t),intent(out)::status
        model=liquid_drop_model_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation','student liquid-drop slot not implemented')
    end subroutine
    subroutine liquid_drop_nuclear_terms(model,a,z,terms,status)
        type(liquid_drop_model_t),intent(in)::model;integer,intent(in)::a,z
        type(liquid_drop_terms_t),intent(out)::terms;type(status_t),intent(out)::status
        terms=liquid_drop_terms_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student liquid-drop slot not implemented')
    end subroutine
    subroutine evaluate_liquid_drop_candidate(model,n_b_m3,a,z,state,status)
        type(liquid_drop_model_t),intent(in)::model;real(wp),intent(in)::n_b_m3;integer,intent(in)::a,z
        type(liquid_drop_state_t),intent(out)::state;type(status_t),intent(out)::status
        state=liquid_drop_state_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student liquid-drop slot not implemented')
    end subroutine
    subroutine scan_liquid_drop_candidates(model,n_b_m3,a_min,a_max,a_step,z_min,z_max,z_step,result,status)
        type(liquid_drop_model_t),intent(in)::model;real(wp),intent(in)::n_b_m3
        integer,intent(in)::a_min,a_max,a_step,z_min,z_max,z_step
        type(liquid_drop_scan_t),intent(out)::result;type(status_t),intent(out)::status
        result=liquid_drop_scan_t()
        call fail_status(status,NS_NOT_IMPLEMENTED,'implementation', &
            'student liquid-drop slot not implemented')
    end subroutine
end module ns_liquid_drop
