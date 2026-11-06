module ns_cowling_cases
    use ns_precision, only: wp
    use ns_status
    use ns_records
    use ns_background_access
    use ns_center_state
    use ns_surface_residual
    use ns_cowling
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    private
    public :: is_cowling_case,run_cowling_case
contains
    logical function is_cowling_case(name)
        character(len=*),intent(in)::name
        is_cowling_case=name=='cowling_center'.or.name=='cowling_rhs_point'.or.name=='background_probe'.or. &
                        name=='surface_residual'.or.name=='cowling_shot'
    end function
    subroutine run_cowling_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token
        character(len=4096)::background_path
        character(len=64)::special_input
        real(wp)::r0_m,amplitude,r_m,r_probe_m,omega_rad_s,w_m3,v_m2,chi,nu,lambda_value,dnu_dr,nu_shift
        real(wp)::ode_rtol,atol_w_m3,atol_v_m2,h_initial,h_min,h_max,dw,dv,qterm,residual,q_over_r2
        integer::ell,max_steps,ios,u,i
        logical::emit_profile
        type(status_t)::status
        type(background_table_t)::bg
        type(background_point_t)::point
        type(cowling_options_t)::options
        type(cowling_shot_t)::shot
        real(wp)::y(2),dydr(2),w,v
        character(len=32)::suffix
        namelist /request/ request_token
        namelist /parameters/ background_path,special_input,r0_m,amplitude,ell,r_m,r_probe_m,omega_rad_s,w_m3,v_m2,chi,nu, &
            lambda_value,dnu_dr,nu_shift,ode_rtol,atol_w_m3,atol_v_m2,h_initial,h_min,h_max,max_steps,emit_profile
        request_token='';background_path='';special_input='none';r0_m=1._wp;amplitude=1._wp;ell=2
        r_m=1000._wp;r_probe_m=1._wp;omega_rad_s=6283.185307179586_wp;w_m3=1e9_wp;v_m2=-1e6_wp
        chi=0.1_wp;nu=-0.1_wp;lambda_value=0.1_wp;dnu_dr=1e-5_wp;nu_shift=0._wp
        ode_rtol=1e-8_wp;atol_w_m3=1e-8_wp;atol_v_m2=1e-8_wp
        h_initial=1e-4_wp;h_min=1e-14_wp;h_max=0.1_wp;max_steps=50000;emit_profile=.false.
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then;read(u,nml=request,iostat=ios);if(ios==0) read(u,nml=parameters,iostat=ios);close(u);end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then;call fail_status(status,NS_INVALID_INPUT,'data','invalid Cowling NAMELIST');call finish();return;end if
        select case(trim(special_input))
        case('none')
        case('nan_r');r0_m=ieee_value(0._wp,ieee_quiet_nan)
        case('nan_w');w_m3=ieee_value(0._wp,ieee_quiet_nan)
        case default;call fail_status(status,NS_INVALID_INPUT,'data','unknown Cowling special input');call finish();return
        end select
        select case(case_id)
        case('cowling_center')
            call cowling_center_state(ell,r0_m,amplitude,w,v,status)
            if(is_ok(status)) then
                call emit_real('w_m3',w);call emit_real('v_m2',v)
                call emit_real('w_over_r3',w/r0_m**3);call emit_real('v_over_r2',v/r0_m**2)
            end if
        case('cowling_rhs_point')
            y=[w_m3,v_m2];call cowling_rhs_point(r_m,omega_rad_s,chi,nu,lambda_value,dnu_dr,y,dydr,qterm,status)
            if(is_ok(status)) then
                call emit_real('dw_dr_m2',dydr(1));call emit_real('dv_dr_m',dydr(2));call emit_real('q_m2',qterm)
            end if
        case('surface_residual')
            call cowling_surface_residual(r_m,omega_rad_s,nu,lambda_value,dnu_dr,w_m3,v_m2,residual,q_over_r2,status)
            if(is_ok(status)) then;call emit_real('residual',residual);call emit_real('q_over_r2',q_over_r2);end if
        case('background_probe')
            call load_background_csv(trim(background_path),bg,status)
            if(is_ok(status)) call background_at_r(bg,r_probe_m,point,status,nu_shift)
            if(is_ok(status)) then
                call emit_real('pressure_pa',point%p_pa);call emit_real('energy_j_m3',point%energy_j_m3)
                call emit_real('chi',point%chi)
                call emit_real('nu',point%nu);call emit_real('lambda_value',point%lambda);call emit_real('dnu_dr',point%dnu_dr)
                call emit_text('interpolation_id',BACKGROUND_INTERPOLATION_ID)
            end if
        case('cowling_shot')
            call load_background_csv(trim(background_path),bg,status)
            if(is_ok(status)) then
                options%ell=ell;options%amplitude=amplitude;options%nu_shift=nu_shift;options%ode%rtol=ode_rtol
                options%ode%atol=[atol_w_m3,atol_v_m2]
                options%ode%h_initial=h_initial;options%ode%h_min=h_min;options%ode%h_max=h_max;options%ode%max_steps=max_steps
                call shoot_cowling(bg,omega_rad_s,options,shot);status=shot%status
            end if
            if(shot%valid.and.is_ok(status)) then
                call emit_real('omega_rad_s',shot%omega_rad_s);call emit_real('residual',shot%residual)
                call emit_real('surface_w_m3',shot%w_m3(size(shot%w_m3)));call emit_real('surface_v_m2',shot%v_m2(size(shot%v_m2)))
                call emit_real('r_start_m',shot%r_m(1));call emit_real('r_end_m',shot%r_m(size(shot%r_m)))
                call emit_int('accepted',int(shot%accepted,int64));call emit_int('rejected',int(shot%rejected,int64))
                call emit_int('evaluations',int(shot%evaluations,int64))
                call emit_text('surface_policy','finite_cutoff_schwarzschild_time')
                call emit_text('interpolation_id',BACKGROUND_INTERPOLATION_ID)
                call emit_bool('mode_calibrated',.false.);call emit_text('identity_status','not_attempted')
                if(emit_profile) then
                    call emit_int('profile_n',int(size(shot%q),int64))
                    do i=1,size(shot%q)
                        write(suffix,'(I5.5)') i
                        call emit_real('q_'//trim(suffix),shot%q(i));call emit_real('r_'//trim(suffix),shot%r_m(i))
                        call emit_real('w_'//trim(suffix),shot%w_m3(i));call emit_real('v_'//trim(suffix),shot%v_m2(i))
                    end do
                end if
            end if
        end select
        call finish()
    contains
        subroutine finish()
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            if(case_id=='cowling_shot') then
                call emit_bool('result_valid',shot%valid.and.is_ok(status))
            else
                call emit_bool('result_valid',is_ok(status))
            end if
        end subroutine
    end subroutine
end module ns_cowling_cases
