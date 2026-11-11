module ns_interface_cases
    use, intrinsic :: iso_fortran_env, only: int64
    use ns_precision,only:wp
    use ns_status
    use ns_records
    use ns_eos_types,only:SIDE_REGULAR,SIDE_INNER,SIDE_OUTER
    use ns_segment_events
    use ns_interface_transfer
    use ns_adaptive_rk, only: adaptive_integrate
    use ns_numerics_types, only: ode_options_t,ode_result_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: is_interface_case,run_interface_case
contains
    logical function is_interface_case(name)
        character(len=*),intent(in)::name
        is_interface_case=name=='material_transfer'.or.name=='rapid_transfer'.or.name=='segment_events'.or. &
            name=='segment_restart_probe'.or.name=='interface_identity'
    end function
    subroutine run_interface_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=32)::scenario,reference_label
        real(wp)::r_m,omega_rad_s,pressure_pa,eps_in,eps_out,nu,lambda_value,dnu_dr,w_in_m3,v_in_m2,overlap
        integer::inner_side,outer_side,ios,u
        logical::candidate_interface_g,root_valid,transfer_valid
        real(wp)::w_out,v_out,dp_in,dp_out,w_rev,v_rev,d1,d2,v_zero,dp_c_in,dp_c_out
        real(wp)::xi_in,xi_out,iso_in,iso_out,pprime_in,pprime_out,q_in,q_out
        type(status_t)::status,st
        type(segment_event_t)::ev(2),z
        real(wp),allocatable::stops(:)
        type(ode_options_t)::opt
        type(ode_result_t)::seg1,seg2,seg3
        real(wp)::w_event_in,v_event_in,w_event_out,v_event_out,dp_event_in,dp_event_out
        real(wp)::chi_w_before,chi_v_before,chi_w_after,chi_v_after,q_end
        integer::nlast
        logical::crossed
        character(len=16)::label
        namelist /request/ request_token
        namelist /parameters/ scenario,reference_label,r_m,omega_rad_s,pressure_pa,eps_in,eps_out,nu,lambda_value,dnu_dr, &
            w_in_m3,v_in_m2,inner_side,outer_side,overlap,candidate_interface_g,root_valid,transfer_valid
        request_token='';scenario='default';reference_label='f';r_m=8000._wp;omega_rad_s=18849.55592153876_wp
        pressure_pa=1e33_wp;eps_in=1.1e34_wp;eps_out=7e33_wp;nu=-0.2_wp;lambda_value=0.15_wp;dnu_dr=1e-4_wp
        w_in_m3=2e11_wp;v_in_m2=-3e7_wp;inner_side=SIDE_INNER;outer_side=SIDE_OUTER;overlap=0.95_wp
        candidate_interface_g=.false.;root_valid=.true.;transfer_valid=.true.
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0)then;read(u,nml=request,iostat=ios);if(ios==0)read(u,nml=parameters,iostat=ios);close(u);end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0)then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid interface NAMELIST')
            call finish();return
        end if
        select case(case_id)
        case('material_transfer')
            call material_interface_transfer(r_m,omega_rad_s,pressure_pa,eps_in,eps_out,nu,lambda_value,dnu_dr,w_in_m3,v_in_m2, &
                 inner_side,outer_side,w_out,v_out,dp_in,dp_out,status)
            if(is_ok(status))then
                call material_interface_transfer(r_m,omega_rad_s,pressure_pa,eps_out,eps_in,nu,lambda_value,dnu_dr,w_out,v_out, &
                    SIDE_INNER,SIDE_OUTER,w_rev,v_rev,d1,d2,st)
                call lagrangian_pressure_jump(r_m,omega_rad_s,pressure_pa,eps_out,nu, &
                    lambda_value,dnu_dr,w_in_m3,v_in_m2,dp_c_out,d1,st)
                call lagrangian_pressure_jump(r_m,omega_rad_s,pressure_pa,eps_in,nu, &
                    lambda_value,dnu_dr,w_in_m3,v_in_m2,dp_c_in,d2,st)
                call material_interface_transfer(r_m,omega_rad_s,pressure_pa,eps_in,eps_in,nu,lambda_value,dnu_dr,w_in_m3,v_in_m2, &
                    SIDE_INNER,SIDE_OUTER,w_rev,v_zero,d1,d2,st)
                call emit_real('w_in_m3',w_in_m3);call emit_real('v_in_m2',v_in_m2)
                call emit_real('w_out_m3',w_out);call emit_real('v_out_m2',v_out)
                call emit_real('delta_p_in_pa',dp_in);call emit_real('delta_p_out_pa',dp_out)
                call emit_real('v_reverse_m2',v_rev)
                call emit_real('v_continuous_delta_p_mismatch_pa',dp_c_out-dp_c_in)
                call emit_bool('zero_jump_identity',abs(v_zero-v_in_m2)<=2e-13_wp*max(abs(v_in_m2),1._wp))
                call emit_bool('same_radius_two_sides',inner_side==SIDE_INNER.and.outer_side==SIDE_OUTER)
            end if
        case('rapid_transfer')
            call rapid_interface_transfer(r_m,omega_rad_s,pressure_pa,eps_in,eps_out,nu,lambda_value,dnu_dr, &
                 w_in_m3,v_in_m2,inner_side,outer_side,w_out,v_out,dp_in,dp_out,status)
            if(is_ok(status))then
                call lagrangian_pressure_jump(r_m,omega_rad_s,pressure_pa,eps_in,nu,lambda_value,dnu_dr, &
                    w_in_m3,v_in_m2,dp_c_in,q_in,st)
                call lagrangian_pressure_jump(r_m,omega_rad_s,pressure_pa,eps_out,nu,lambda_value,dnu_dr, &
                    w_out,v_out,dp_c_out,q_out,st)
                xi_in=exp(-lambda_value)*w_in_m3/(r_m*r_m)
                xi_out=exp(-lambda_value)*w_out/(r_m*r_m)
                pprime_in=-(eps_in+pressure_pa)*dnu_dr
                pprime_out=-(eps_out+pressure_pa)*dnu_dr
                iso_in=xi_in-dp_in/pprime_in
                iso_out=xi_out-dp_out/pprime_out
                call emit_real('w_in_m3',w_in_m3);call emit_real('v_in_m2',v_in_m2)
                call emit_real('w_out_m3',w_out);call emit_real('v_out_m2',v_out)
                call emit_real('delta_p_in_pa',dp_in);call emit_real('delta_p_out_pa',dp_out)
                call emit_real('isobar_displacement_in_m',iso_in);call emit_real('isobar_displacement_out_m',iso_out)
                call emit_real('delta_p_relative_jump',abs(dp_out-dp_in)/max(abs(dp_out),abs(dp_in),1._wp))
                call emit_real('rapid_condition_residual',abs(iso_out-iso_in)/max(abs(iso_out),abs(iso_in),1._wp))
            end if
        case('segment_events')
            select case(trim(scenario))
            case('default')
                call make_segment_event(1e34_wp,1e22_wp,1e33_wp,EVENT_MATERIAL,SIDE_INNER,SIDE_OUTER,ev(1),status)
                if(is_ok(status)) call make_segment_event(1e34_wp,1e22_wp,5e32_wp, &
                    EVENT_CHI_KNOT,SIDE_OUTER,SIDE_OUTER,ev(2),status)
            case('unsorted')
                call make_segment_event(1e34_wp,1e22_wp,5e32_wp,EVENT_MATERIAL,SIDE_INNER,SIDE_OUTER,ev(1),status)
                if(is_ok(status)) call make_segment_event(1e34_wp,1e22_wp,1e33_wp, &
                    EVENT_CHI_KNOT,SIDE_OUTER,SIDE_OUTER,ev(2),status)
            case('same_side_material')
                call make_segment_event(1e34_wp,1e22_wp,1e33_wp,EVENT_MATERIAL,SIDE_OUTER,SIDE_OUTER,ev(1),status)
            case('chi_side_change')
                call make_segment_event(1e34_wp,1e22_wp,1e33_wp,EVENT_CHI_KNOT,SIDE_INNER,SIDE_OUTER,ev(1),status)
            case default;call fail_status(status,NS_INVALID_INPUT,'data','unknown event scenario')
            end select
            if(is_ok(status).and.trim(scenario)/='same_side_material'.and. &
                trim(scenario)/='chi_side_change') call validate_segment_events(ev,status)
            if(is_ok(status))then
                call event_stops(ev,stops,status)
                if(is_ok(status))then
                    call emit_bool('events_valid',.true.);call emit_real('q1',stops(1));call emit_real('q2',stops(2))
                    call emit_bool('cross_stage_allowed',.false.)
                    call emit_bool('material_side_changes',ev(1)%inner_side/=ev(1)%outer_side)
                    call emit_bool('chi_knot_state_continuous',ev(2)%inner_side==ev(2)%outer_side)
                end if
            end if
        case('segment_restart_probe')
            call make_segment_event(1e34_wp,1e22_wp,1e33_wp,EVENT_MATERIAL,SIDE_INNER,SIDE_OUTER,ev(1),status)
            if(is_ok(status)) call make_segment_event(1e34_wp,1e22_wp,5e32_wp, &
                EVENT_CHI_KNOT,SIDE_OUTER,SIDE_OUTER,ev(2),status)
            if(is_ok(status)) call event_stops(ev,stops,status)
            if(.not.is_ok(status)) then
                call finish();return
            end if
            opt%atol=[1e-6_wp,1e-6_wp];opt%rtol=1e-10_wp
            opt%h_initial=1._wp;opt%h_min=1e-12_wp;opt%h_max=10._wp
            opt%max_steps=64;opt%record_trials=.true.
            call adaptive_integrate(zero_rhs,0._wp,stops(1),[w_in_m3,v_in_m2],opt,seg1)
            if(.not.seg1%valid.or..not.is_ok(seg1%status)) then
                status=seg1%status;call finish();return
            end if
            nlast=size(seg1%x);w_event_in=seg1%y(1,nlast);v_event_in=seg1%y(2,nlast)
            call material_interface_transfer(r_m,omega_rad_s,pressure_pa,eps_in,eps_out,nu,lambda_value,dnu_dr, &
                w_event_in,v_event_in,SIDE_INNER,SIDE_OUTER,w_event_out,v_event_out,dp_event_in,dp_event_out,status)
            if(.not.is_ok(status)) then;call finish();return;end if
            call adaptive_integrate(zero_rhs,stops(1),stops(2),[w_event_out,v_event_out],opt,seg2)
            if(.not.seg2%valid.or..not.is_ok(seg2%status)) then
                status=seg2%status;call finish();return
            end if
            nlast=size(seg2%x);chi_w_before=seg2%y(1,nlast);chi_v_before=seg2%y(2,nlast)
            ! A chi-only knot changes the derivative law, not the state itself.
            chi_w_after=chi_w_before;chi_v_after=chi_v_before
            q_end=stops(2)+1._wp
            call adaptive_integrate(zero_rhs,stops(2),q_end,[chi_w_after,chi_v_after],opt,seg3)
            if(.not.seg3%valid.or..not.is_ok(seg3%status)) then
                status=seg3%status;call finish();return
            end if
            crossed=segment_crossed(seg1,stops(1)).or.segment_crossed(seg2,stops(2))
            if(seg1%x(size(seg1%x))/=stops(1).or.seg2%x(size(seg2%x))/=stops(2).or.seg3%x(size(seg3%x))/=q_end) then
                call fail_status(status,NS_STAGNATION,'numerics','segmented restart did not land exactly on an endpoint')
                call finish();return
            end if
            call emit_int('segments_completed',3_int64)
            call emit_bool('material_event_hit',.true.);call emit_bool('chi_event_hit',.true.)
            call emit_bool('any_segment_crossed_event',crossed)
            call emit_real('material_q_inner',stops(1));call emit_real('material_q_outer',stops(1))
            call emit_real('material_r_inner_m',r_m);call emit_real('material_r_outer_m',r_m)
            call emit_int('material_side_inner',int(SIDE_INNER,int64));call emit_int('material_side_outer',int(SIDE_OUTER,int64))
            call emit_real('material_w_inner_m3',w_event_in);call emit_real('material_w_outer_m3',w_event_out)
            call emit_real('material_v_inner_m2',v_event_in);call emit_real('material_v_outer_m2',v_event_out)
            call emit_real('material_delta_p_inner_pa',dp_event_in);call emit_real('material_delta_p_outer_pa',dp_event_out)
            call emit_real('chi_q_before',stops(2));call emit_real('chi_q_after',stops(2))
            call emit_real('chi_w_before_m3',chi_w_before);call emit_real('chi_w_after_m3',chi_w_after)
            call emit_real('chi_v_before_m2',chi_v_before);call emit_real('chi_v_after_m2',chi_v_after)
            call clear_status(status)
        case('interface_identity')
            call identify_interface_family(trim(reference_label),root_valid,transfer_valid, &
                overlap,candidate_interface_g,label,status)
            call emit_text('mode_label',trim(label))
            if(is_ok(status))then
                call emit_text('identity_status','identified')
            else
                call emit_text('identity_status','unresolved')
            end if
        end select
        call finish()
    contains
        subroutine zero_rhs(q,state,dstate,status_local)
            real(wp),intent(in)::q,state(:)
            real(wp),intent(out)::dstate(:)
            type(status_t),intent(out)::status_local
            if(.not.ieee_is_finite(q).or..not.all(ieee_is_finite(state)).or.size(dstate)/=size(state)) then
                dstate=0._wp;call fail_status(status_local,NS_INVALID_INPUT,'numerics','invalid segmented probe state');return
            end if
            dstate=0._wp;call clear_status(status_local)
        end subroutine
        logical function segment_crossed(seg,event_q)
            type(ode_result_t),intent(in)::seg
            real(wp),intent(in)::event_q
            integer::j
            segment_crossed=.false.
            if(.not.allocated(seg%trial_x).or..not.allocated(seg%trial_h)) return
            do j=1,size(seg%trial_x)
                if((seg%trial_x(j)<event_q.and.seg%trial_x(j)+seg%trial_h(j)>event_q).or. &
                   (seg%trial_x(j)>event_q.and.seg%trial_x(j)+seg%trial_h(j)<event_q)) then
                    segment_crossed=.true.;return
                end if
            end do
        end function
        subroutine finish()
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            call emit_bool('result_valid',is_ok(status))
        end subroutine
    end subroutine
end module ns_interface_cases
