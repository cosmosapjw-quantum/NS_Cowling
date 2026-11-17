module ns_mode_cases
    use ns_precision, only: wp
    use ns_physical_constants, only: two_pi
    use ns_status
    use ns_records
    use ns_background_access
    use ns_cowling, only: cowling_shot_t
    use ns_scan, only: linear_scan_grid
    use ns_scan_types, only: bracket_candidate_t,CAND_SIGN,CAND_SAMPLE
    use ns_mode_spectrum
    use ns_production_solvers, only: solve_cowling_spectrum
    use ns_mode_normalization
    use ns_mode_identity
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private
    public::is_mode_case,run_mode_case
contains
    logical function is_mode_case(name)
        character(len=*),intent(in)::name
        is_mode_case=name=='smooth_spectrum'.or.name=='mode_refine'.or.name=='mode_trial'.or. &
            name=='mode_diagnostics'.or.name=='mode_identity_probe'
    end function
    subroutine run_mode_case(case_id,path)
        character(len=*),intent(in)::case_id,path
        character(len=32)::request_token,suffix
        character(len=4096)::background_path
        character(len=24)::method
        character(len=16)::label
        character(len=48)::diagnostic_action
        real(wp)::frequency_min_hz,frequency_max_hz,frequency_hz,ode_rtol,amplitude,nu_shift, &
            root_xtol,root_rtol,residual_tol,node_zero,node_xlo,node_xhi,h_max, &
            norm2,op,on,r(101),w(101),v(101)
        real(wp),allocatable::grid(:),x(:),wn(:),vn(:),xp(:),wpn(:),vpn(:),xn(:),wnn(:),vnn(:)
        integer::scan_n,cell_substeps,max_steps,max_iter,u,ios,i,n,j,wnodes,vnodes,candidate_kind
        logical::emit_profile,root_valid,smooth,diagnostics_valid
        type(status_t)::status,st
        type(background_table_t)::bg
        type(mode_options_t)::opt
        type(mode_spectrum_t)::spectrum
        type(mode_result_t)::rr
        type(cowling_shot_t)::shot
        type(bracket_candidate_t)::candidate
        namelist /request/ request_token
        namelist /parameters/ background_path,method,frequency_min_hz,frequency_max_hz,frequency_hz,scan_n, &
            ode_rtol,amplitude,nu_shift,root_xtol,root_rtol,residual_tol,node_zero,node_xlo,node_xhi, &
            h_max,cell_substeps,max_steps,max_iter,emit_profile,diagnostic_action
        request_token='';background_path='';method='adaptive';diagnostic_action='regular'
        frequency_min_hz=400._wp;frequency_max_hz=9000._wp;frequency_hz=1500._wp;scan_n=65
        ode_rtol=1e-9_wp;amplitude=1._wp;nu_shift=0._wp;root_xtol=1e-10_wp;root_rtol=0._wp
        residual_tol=1e-5_wp;node_zero=1e-7_wp;node_xlo=0.001_wp;node_xhi=0.999_wp
        h_max=0.05_wp;cell_substeps=8;max_steps=200000;max_iter=100;emit_profile=.false.
        call emit_text('case_id',case_id)
        open(newunit=u,file=path,status='old',action='read',iostat=ios)
        if(ios==0) then
            read(u,nml=request,iostat=ios)
            if(ios==0) read(u,nml=parameters,iostat=ios)
            close(u)
        end if
        call emit_text('request_token',trim(request_token))
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','invalid T08 NAMELIST');call finish(.false.);return
        end if
        if(case_id=='mode_identity_probe') then
            call clear_status(status);call emit_text('root_status','ok')
            call identify_smooth_mode(.true.,.true.,.true.,0,0,label,st)
            call emit_text('zero_nodes_label',label)
            call identify_smooth_mode(.true.,.true.,.true.,1,1,label,st)
            call emit_text('one_nodes_label',label)
            call identify_smooth_mode(.true.,.true.,.true.,1,0,label,st)
            call emit_text('ambiguous_identity',trim(status_name(st%code)))
            call identify_smooth_mode(.true.,.false.,.true.,0,0,label,st)
            call emit_text('jump_identity',trim(status_name(st%code)))
            call identify_smooth_mode(.false.,.true.,.true.,0,0,label,st)
            call emit_text('unconverged_identity',trim(status_name(st%code)))
            call finish(.true.);return
        end if
        if(case_id=='mode_diagnostics') then
            do i=1,101
                r(i)=0.001_wp+real(i-1,wp)*0.00999_wp
                w(i)=r(i)**3*(1._wp-2._wp*r(i));v(i)=r(i)**2*(1._wp-2._wp*r(i))
            end do
            select case(trim(diagnostic_action))
            case('regular')
            case('zero');w=0._wp;v=0._wp
            case('duplicate_r');r(51)=r(50)
            case('duplicate_identical');r(51)=r(50);w(51)=w(50);v(51)=v(50)
            case('extra_node');w=w*(r-0.75_wp)
            case default
                call fail_status(status,NS_INVALID_INPUT,'data','unknown diagnostic action');call finish(.false.);return
            end select
            call normalize_mode(r,w,v,x,wn,vn,status)
            if(.not.is_ok(status)) then;call finish(.false.);return;end if
            call normalize_mode(r,1e40_wp*w,1e40_wp*v,xp,wpn,vpn,st)
            if(.not.is_ok(st)) then;status=st;call finish(.false.);return;end if
            call normalize_mode(r,-w,-v,xn,wnn,vnn,st)
            norm2=0._wp;op=0._wp;on=0._wp
            do i=2,size(x)
                norm2=norm2+0.5_wp*(x(i)-x(i-1))*(wn(i)**2+wn(i-1)**2+6._wp*(vn(i)**2+vn(i-1)**2))
                op=op+0.5_wp*(x(i)-x(i-1))*(wn(i)*wpn(i)+wn(i-1)*wpn(i-1)+ &
                    6._wp*(vn(i)*vpn(i)+vn(i-1)*vpn(i-1)))
                on=on+0.5_wp*(x(i)-x(i-1))*(wn(i)*wnn(i)+wn(i-1)*wnn(i-1)+ &
                    6._wp*(vn(i)*vnn(i)+vn(i-1)*vnn(i-1)))
            end do
            call count_smooth_nodes(x,wn,node_xlo,node_xhi,node_zero,wnodes,status)
            call count_smooth_nodes(x,vn,node_xlo,node_xhi,node_zero,vnodes,st)
            if(.not.is_ok(st)) status=st
            call emit_real('norm_squared',norm2);call emit_real('positive_overlap',op);call emit_real('negative_overlap',on)
            call emit_int('w_nodes',int(wnodes,int64));call emit_int('v_nodes',int(vnodes,int64))
            call finish(is_ok(status));return
        end if
        call load_background_csv(trim(background_path),bg,status)
        if(.not.is_ok(status)) then;call finish(.false.);return;end if
        opt%method=method;opt%cell_substeps=cell_substeps;opt%residual_tol=residual_tol
        opt%root%xtol=root_xtol;opt%root%rtol=root_rtol;opt%root%ftol=0._wp;opt%root%max_iter=max_iter
        opt%node_zero=node_zero;opt%node_xlo=node_xlo;opt%node_xhi=node_xhi
        opt%shot%amplitude=amplitude;opt%shot%nu_shift=nu_shift
        opt%shot%ode%rtol=ode_rtol;opt%shot%ode%atol=[1e-10_wp,1e-10_wp]
        opt%shot%ode%h_initial=1e-6_wp;opt%shot%ode%h_min=1e-16_wp;opt%shot%ode%h_max=h_max
        opt%shot%ode%max_steps=max_steps
        call emit_text('interpolation_id',BACKGROUND_INTERPOLATION_ID)
        call emit_text('normalization_id',MODE_NORM_ID);call emit_text('identity_rule',MODE_IDENTITY_ID)
        call emit_text('integrator',trim(method));call emit_bool('exact_vacuum_surface',.false.)
        call emit_real('r0_m',bg%r_m(1));call emit_real('radius_m',bg%r_m(bg%n))
        call emit_real('mass_kg',bg%m_kg(bg%n));call emit_real('p_cut_pa',bg%p_pa(bg%n))
        select case(case_id)
        case('mode_trial')
            call evaluate_mode_shot(bg,frequency_hz*two_pi,opt,shot);status=shot%status
            if(shot%valid) then
                call emit_real('residual',shot%residual);call emit_int('accepted',int(shot%accepted,int64))
                call emit_int('evaluations',int(shot%evaluations,int64))
            end if
            call finish(shot%valid)
        case('mode_refine')
            candidate=bracket_candidate_t(frequency_min_hz*two_pi,frequency_max_hz*two_pi,1,2,CAND_SIGN)
            if(trim(diagnostic_action)=='sample') candidate%kind=CAND_SAMPLE
            call refine_mode(bg,candidate,opt,rr);status=rr%root_status
            call write_root(rr,'mode_');call finish(rr%root_valid)
        case('smooth_spectrum')
            call linear_scan_grid(frequency_min_hz*two_pi,frequency_max_hz*two_pi,scan_n,grid,status)
            if(.not.is_ok(status)) then;call finish(.false.);return;end if
            call solve_cowling_spectrum(bg,grid,opt,spectrum);status=spectrum%status
            call emit_int('roots_found',int(spectrum%roots_found,int64))
            call emit_int('f_index',int(spectrum%f_index,int64));call emit_int('p1_index',int(spectrum%p1_index,int64))
            if(allocated(spectrum%scan%trials)) then
                call emit_int('scan_n',int(size(spectrum%scan%trials),int64))
                call emit_int('barriers',int(spectrum%scan%barriers,int64))
                call emit_int('sign_brackets',int(spectrum%scan%sign_brackets,int64))
                do i=1,size(spectrum%scan%trials)
                    write(suffix,'(I5.5)')i
                    call emit_real('scan_hz_'//trim(suffix),spectrum%scan%trials(i)%parameter/two_pi)
                    call emit_text('scan_status_'//trim(suffix),trim(status_name(spectrum%scan%trials(i)%status%code)))
                    if(spectrum%scan%trials(i)%valid) call emit_real('scan_residual_'//trim(suffix), &
                spectrum%scan%trials(i)%residual)
                end do
            end if
            if(allocated(spectrum%refinements)) then
                do j=1,size(spectrum%refinements)
                    write(suffix,'(I5.5)')j
                    call write_root(spectrum%refinements(j),'candidate_'//trim(suffix)//'_')
                end do
            end if
            if(spectrum%f_index>0) call write_root(spectrum%refinements(spectrum%f_index),'f_')
            if(spectrum%p1_index>0) call write_root(spectrum%refinements(spectrum%p1_index),'p1_')
            call finish(spectrum%valid)
        end select
    contains
        subroutine write_root(m,prefix)
            type(mode_result_t),intent(in)::m
            character(len=*),intent(in)::prefix
            integer::idx,step,nm
            character(len=32)::tag
            call emit_text(prefix//'root_status',trim(status_name(m%root_status%code)))
            call emit_text(prefix//'label',trim(m%label))
            if(is_ok(m%identity_status)) then
                call emit_text(prefix//'identity','identified')
            else
                call emit_text(prefix//'identity',trim(status_name(m%identity_status%code)))
            end if
            call emit_bool(prefix//'root_valid',m%root_valid)
            if(m%search%valid) then
                call emit_real(prefix//'hz',m%search%root/two_pi)
                call emit_real(prefix//'bracket_width_hz',(m%search%hi-m%search%lo)/two_pi)
                call emit_int(prefix//'root_iterations',int(m%search%iterations,int64))
            end if
            if(m%shot%valid) then
                call emit_real(prefix//'residual',m%shot%residual)
                nm=size(m%shot%r_m)
                call emit_real(prefix//'min_dr_m',minval(m%shot%r_m(2:)-m%shot%r_m(:nm-1)))
                call emit_int(prefix//'duplicate_r_count',int(count(m%shot%r_m(2:)==m%shot%r_m(:nm-1)),int64))
            end if
            call emit_int(prefix//'w_nodes',int(m%w_nodes,int64));call emit_int(prefix//'v_nodes',int(m%v_nodes,int64))
            if(emit_profile.and.m%diagnostics_valid) then
                ! Full profile, no hidden thinning or interpolation.
                nm=size(m%x);call emit_int(prefix//'profile_n',int(nm,int64))
                do idx=1,nm
                    write(tag,'(I7.7)')idx
                    call emit_real(prefix//'x_'//trim(tag),m%x(idx))
                    call emit_real(prefix//'wn_'//trim(tag),m%wn(idx))
                    call emit_real(prefix//'vn_'//trim(tag),m%vn(idx))
                end do
            end if
        end subroutine
        subroutine finish(valid)
            logical,intent(in)::valid
            call emit_text('status',trim(status_name(status%code)))
            call emit_text('category',trim(status%category));call emit_text('detail',trim(status%detail))
            call emit_bool('result_valid',valid.and.is_ok(status))
        end subroutine
    end subroutine
end module ns_mode_cases
