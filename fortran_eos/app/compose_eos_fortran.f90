program compose_eos_fortran
  use compose_eos_types
  use compose_eos_component, only: extract_core_rows, extract_crust_rows, verify_core_roundtrip
  use compose_eos_io
  use compose_eos_verify
  use compose_eos_reconstruct, only: fill_density_gap, reconstruct_master
  use compose_eos_causality, only: causal_complete, causal_completion_report_t
  use compose_eos_matching, only: match_crust_core, match_report_t
  use compose_eos_phase_transition, only: construct_maxwell_transition, maxwell_transition_report_t
  implicit none
  type(eos_data_t)::eos,out,crust,core,phase_low,phase_high,readback
  type(verify_report_t)::report
  type(causal_completion_report_t)::causal_report
  type(match_report_t)::match_report
  type(maxwell_transition_report_t)::phase_report
  character(len=32)::action
  character(len=2048)::arg1,arg2,arg3
  character(len=128)::arg4,arg5,arg6,arg7
  character(len=512)::status
  integer::narg,ios,n_fill,fit_points,gap_phase,n_points,n_bridge
  integer::core_first,core_last,expected_source_n,expected_phase
  integer::crust_outer_last,expected_outer_phase,expected_inner_phase
  real(wp)::cs2_extension,n_max_requested,smooth_sigma,dissolution_tol
  real(wp)::min_ratio_local,best_sel
  integer::igap_sel
  logical::has_gap_phase,drop_composition,local_neutrality_confirmed
  logical::smooth_master_flag,fix_dissolution,recompute_extra_cs2,outer_crust_corrected
  logical::has_min_ratio
  character(len=256)::why
  action='';arg1='';arg2='';arg3='';arg4='';arg5='';arg6='';arg7=''
  narg=command_argument_count()
  if(narg<2)then
    write(*,'(a)')'usage: compose_eos_fortran {read|verify|verify-csv|normalize|write|extract-core|extract-crust|fill-gap|reconstruct-master|causal-complete|match-core|maxwell-transition} INPUT OUTPUT [OPTIONS]'
    stop 2
  endif
  call get_command_argument(1,action,status=ios)
  call get_command_argument(2,arg1,status=ios)
  if(narg>=3)call get_command_argument(3,arg2,status=ios)
  if(narg>=4)call get_command_argument(4,arg3,status=ios)
  if(narg>=5)call get_command_argument(5,arg4,status=ios)
  if(narg>=6)call get_command_argument(6,arg5,status=ios)
  if(narg>=7)call get_command_argument(7,arg6,status=ios)
  if(narg>=8)call get_command_argument(8,arg7,status=ios)
  select case(trim(action))
  case('read')
    if(narg/=3)then;call emit_failure('read requires input CompOSE directory and output directory');stop 2;endif
    call read_compose(trim(arg1),eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_internal_csv(trim(arg2)//'/eos_data.csv',eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_metadata(trim(arg2)//'/eos_data.nml',eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_normalized_csv(trim(arg2)//'/normalized_eos.csv',eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(*,'(a)')'{"status":"ok","operation":"read","n":'//integer_text(eos%n)//'}'
  case('verify')
    if(narg/=3)then;call emit_failure('verify requires CompOSE directory and JSON output path');stop 2;endif
    call read_compose(trim(arg1),eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call check_consistency(eos,report)
    call write_verify_json(trim(arg2),eos,report)
    if(report%overall)then
      write(*,'(a)')'{"status":"ok","operation":"verify","overall":"PASS"}'
    else
      write(*,'(a)')'{"status":"ok","operation":"verify","overall":"FAIL","failed_checks":"'//trim(report%failed_checks)//'"}'
    endif
  case('verify-csv')
    if(narg/=3)then;call emit_failure('verify-csv requires internal EOS directory and JSON output');stop 2;endif
    call read_internal_csv(trim(arg1)//'/eos_data.csv',trim(arg1)//'/eos_data.nml',eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call check_consistency(eos,report)
    call write_verify_json(trim(arg2),eos,report)
    if(report%overall)then
      write(*,'(a)')'{"status":"ok","operation":"verify-csv","overall":"PASS"}'
    else
      write(*,'(a)')'{"status":"ok","operation":"verify-csv","overall":"FAIL","failed_checks":"'//trim(report%failed_checks)//'"}'
    endif
  case('normalize')
    if(narg/=3)then;call emit_failure('normalize requires CompOSE directory and normalized CSV path');stop 2;endif
    call read_compose(trim(arg1),eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_normalized_csv(trim(arg2),eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(*,'(a)')'{"status":"ok","operation":"normalize","n":'//integer_text(eos%n)//'}'
  case('write')
    if(narg/=3)then;call emit_failure('write requires internal EOS directory and output CompOSE directory');stop 2;endif
    call read_internal_csv(trim(arg1)//'/eos_data.csv',trim(arg1)//'/eos_data.nml',eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_compose(trim(arg2),eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(*,'(a)')'{"status":"ok","operation":"write","n":'//integer_text(eos%n)//'}'
  case('fill-gap')
    if(narg<5.or.narg>7)then
      call emit_failure('fill-gap requires input internal EOS directory, output directory, n_fill, fit_points, and optional gap_phase and min_ratio');stop 2
    endif
    read(arg3,*,iostat=ios)n_fill
    if(ios/=0)then;call emit_failure('invalid n_fill');stop 2;endif
    read(arg4,*,iostat=ios)fit_points
    if(ios/=0)then;call emit_failure('invalid fit_points');stop 2;endif
    has_gap_phase=.false.;gap_phase=0
    min_ratio_local=5._wp
    has_min_ratio=.false.
    ! remaining args: gap_phase, then min_ratio (or just min_ratio)
    if(narg>=6)then
      ! arg5 could be gap_phase (int) or min_ratio (real); disambiguate by presence of a '.'
      if(index(arg5,'.')>0.or.index(arg5,'e')>0.or.index(arg5,'E')>0)then
        read(arg5,*,iostat=ios)min_ratio_local
        if(ios/=0)then;call emit_failure('invalid min_ratio');stop 2;endif
        has_min_ratio=.true.
        if(narg==7)then;call emit_failure('fill-gap: gap_phase and min_ratio both given but gap_phase must come first');stop 2;endif
      else
        has_gap_phase=.true.
        read(arg5,*,iostat=ios)gap_phase
        if(ios/=0)then;call emit_failure('invalid gap_phase');stop 2;endif
        if(narg==7)then
          read(arg6,*,iostat=ios)min_ratio_local
          if(ios/=0)then;call emit_failure('invalid min_ratio');stop 2;endif
          has_min_ratio=.true.
        endif
      endif
    endif
    call read_internal_csv(trim(arg1)//'/eos_data.csv',trim(arg1)//'/eos_data.nml',eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    if(has_min_ratio)then
      call fill_density_gap(eos,n_fill,fit_points,gap_phase,has_gap_phase,out,status,min_ratio_local,igap_sel,best_sel)
    else
      call fill_density_gap(eos,n_fill,fit_points,gap_phase,has_gap_phase,out,status,igap_out=igap_sel,best_out=best_sel)
    endif
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_internal_csv(trim(arg2)//'/eos_data.csv',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_metadata(trim(arg2)//'/eos_data.nml',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(why,'(a,i0,a,i0,a,g0.6,a)')'{"status":"ok","operation":"fill-gap","n":',out%n,&
      ',"igap":',igap_sel,',"best_ratio":',best_sel,'}'
    write(*,'(a)')trim(why)
  case('causal-complete')
    if(narg/=6)then
      call emit_failure('causal-complete requires input internal directory, output directory, cs2_extension, n_max (zero selects source maximum), and n_points');stop 2
    endif
    read(arg3,*,iostat=ios)cs2_extension
    if(ios/=0)then;call emit_failure('invalid cs2_extension');stop 2;endif
    read(arg4,*,iostat=ios)n_max_requested
    if(ios/=0)then;call emit_failure('invalid n_max');stop 2;endif
    read(arg5,*,iostat=ios)n_points
    if(ios/=0)then;call emit_failure('invalid n_points');stop 2;endif
    call read_internal_csv(trim(arg1)//'/eos_data.csv',trim(arg1)//'/eos_data.nml',eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call causal_complete(eos,cs2_extension,n_max_requested,n_points,out,causal_report,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_internal_csv(trim(arg2)//'/eos_data.csv',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_metadata(trim(arg2)//'/eos_data.nml',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(*,'(a)')'{"status":"ok","operation":"causal-complete","n":'//integer_text(out%n)// &
      ',"crossing_left_index":'//integer_text(causal_report%crossing_left_index)// &
      ',"n_match":'//real_text(causal_report%n_match)// &
      ',"pressure_match":'//real_text(causal_report%pressure_match)// &
      ',"energy_match":'//real_text(causal_report%energy_match)// &
      ',"mu_match":'//real_text(causal_report%mu_match)// &
      ',"source_mu_match":'//real_text(causal_report%source_mu_match)// &
      ',"mu_match_relative_residual":'//real_text(causal_report%mu_match_relative_residual)// &
      ',"source_max_cs2":'//real_text(causal_report%source_max_cs2)// &
      ',"cs2_extension":'//real_text(causal_report%cs2_extension)// &
      ',"n_max":'//real_text(causal_report%n_max)// &
      ',"n_tail_points":'//integer_text(causal_report%n_tail_points)// &
      ',"composition_omitted":'//logical_text(causal_report%composition_omitted)// &
      ',"sound_speed_continuous":'//logical_text(causal_report%sound_speed_continuous)//'}'
  case('match-core')
    if(narg/=7)then
      call emit_failure('match-core requires CRUST_INTERNAL CORE_INTERNAL OUTPUT METHOD N_BRIDGE DROP_COMPOSITION');stop 2
    endif
    read(arg5,*,iostat=ios)n_bridge
    if(ios/=0)then;call emit_failure('invalid n_bridge');stop 2;endif
    select case(trim(adjustl(arg6)))
    case('true','TRUE','.true.','1')
      drop_composition=.true.
    case('false','FALSE','.false.','0')
      drop_composition=.false.
    case default
      call emit_failure('DROP_COMPOSITION must be true or false');stop 2
    end select
    call read_internal_csv(trim(arg1)//'/eos_data.csv',trim(arg1)//'/eos_data.nml',crust,status)
    if(status/='ok')then;call emit_failure('invalid crust input: '//trim(status));stop 3;endif
    call read_internal_csv(trim(arg2)//'/eos_data.csv',trim(arg2)//'/eos_data.nml',core,status)
    if(status/='ok')then;call emit_failure('invalid core input: '//trim(status));stop 3;endif
    call match_crust_core(crust,core,trim(arg4),n_bridge,drop_composition,out,match_report,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_internal_csv(trim(arg3)//'/eos_data.csv',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_metadata(trim(arg3)//'/eos_data.nml',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(*,'(a)')'{"status":"ok","operation":"match-core","n":'//integer_text(out%n)// &
      ',"method":"'//trim(match_report%method)//'","construction_kind":"'//trim(match_report%construction_kind)// &
      '","spline_mode":"'//trim(match_report%spline_mode)//'","i1":'//integer_text(match_report%i1)// &
      ',"i2":'//integer_text(match_report%i2)//',"n_bridge":'//integer_text(match_report%n_bridge)// &
      ',"n1":'//real_text(match_report%n1)//',"n2":'//real_text(match_report%n2)// &
      ',"P1":'//real_text(match_report%p1)//',"P2":'//real_text(match_report%p2)// &
      ',"F1":'//real_text(match_report%f1)//',"F2":'//real_text(match_report%f2)// &
      ',"mu1":'//real_text(match_report%mu1)//',"mu2":'//real_text(match_report%mu2)// &
      ',"cs2_source_left":'//real_text(match_report%cs2_source_left)// &
      ',"cs2_source_right":'//real_text(match_report%cs2_source_right)// &
      ',"cs2_left":'//real_text(match_report%cs2_left)//',"cs2_right":'//real_text(match_report%cs2_right)// &
      ',"cs2_min":'//real_text(match_report%cs2_min)//',"cs2_max":'//real_text(match_report%cs2_max)// &
      ',"pressure_span":'//real_text(match_report%pressure_span)// &
      ',"first_law_relative_residual":'//real_text(match_report%first_law_relative_residual)// &
      ',"composition_dropped":'//logical_text(match_report%composition_dropped)// &
      ',"is_exact_maxwell":'//logical_text(match_report%is_exact_maxwell)// &
      ',"seidov_applicable":'//logical_text(match_report%seidov_applicable)//'}'
  case('reconstruct-master')
    if(narg/=8)then
      call emit_failure('reconstruct-master requires INPUT_INTERNAL OUTPUT_INTERNAL SMOOTH SMOOTH_SIGMA FIX_DISSOLUTION DISSOLUTION_TOL RECOMPUTE_EXTRA_CS2');stop 2
    endif
    select case(trim(adjustl(arg3)))
    case('true','TRUE','.true.','1')
      smooth_master_flag=.true.
    case('false','FALSE','.false.','0')
      smooth_master_flag=.false.
    case default
      call emit_failure('SMOOTH must be true or false');stop 2
    end select
    read(arg4,*,iostat=ios)smooth_sigma
    if(ios/=0)then;call emit_failure('invalid SMOOTH_SIGMA');stop 2;endif
    select case(trim(adjustl(arg5)))
    case('true','TRUE','.true.','1')
      fix_dissolution=.true.
    case('false','FALSE','.false.','0')
      fix_dissolution=.false.
    case default
      call emit_failure('FIX_DISSOLUTION must be true or false');stop 2
    end select
    read(arg6,*,iostat=ios)dissolution_tol
    if(ios/=0)then;call emit_failure('invalid DISSOLUTION_TOL');stop 2;endif
    select case(trim(adjustl(arg7)))
    case('true','TRUE','.true.','1')
      recompute_extra_cs2=.true.
    case('false','FALSE','.false.','0')
      recompute_extra_cs2=.false.
    case default
      call emit_failure('RECOMPUTE_EXTRA_CS2 must be true or false');stop 2
    end select
    call read_internal_csv(trim(arg1)//'/eos_data.csv',trim(arg1)//'/eos_data.nml',eos,status)
    if(status/='ok')then;call emit_failure('invalid EOS input: '//trim(status));stop 3;endif
    call reconstruct_master(eos,smooth_master_flag,smooth_sigma,fix_dissolution, &
      dissolution_tol,recompute_extra_cs2,out,outer_crust_corrected,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_internal_csv(trim(arg2)//'/eos_data.csv',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_metadata(trim(arg2)//'/eos_data.nml',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(*,'(a)')'{"status":"ok","operation":"reconstruct-master","n":'//integer_text(out%n)// &
      ',"composition_reconstructed":'//logical_text(out%composition_available)// &
      ',"outer_crust_correction_applied":'//logical_text(outer_crust_corrected)// &
      ',"q_extra_recomputed_as_cs2":'//logical_text(recompute_extra_cs2)//'}'
  case('maxwell-transition')
    if(narg/=5)then
      call emit_failure('maxwell-transition requires LOW_PHASE_INTERNAL HIGH_PHASE_INTERNAL OUTPUT_INTERNAL LOCAL_NEUTRALITY_CONFIRMED');stop 2
    endif
    select case(trim(adjustl(arg4)))
    case('true','TRUE','.true.','1')
      local_neutrality_confirmed=.true.
    case('false','FALSE','.false.','0')
      local_neutrality_confirmed=.false.
    case default
      call emit_failure('LOCAL_NEUTRALITY_CONFIRMED must be true or false');stop 2
    end select
    call read_internal_csv(trim(arg1)//'/eos_data.csv',trim(arg1)//'/eos_data.nml',phase_low,status)
    if(status/='ok')then;call emit_failure('invalid low-density phase: '//trim(status));stop 3;endif
    call read_internal_csv(trim(arg2)//'/eos_data.csv',trim(arg2)//'/eos_data.nml',phase_high,status)
    if(status/='ok')then;call emit_failure('invalid high-density phase: '//trim(status));stop 3;endif
    call construct_maxwell_transition(phase_low,phase_high,local_neutrality_confirmed,out,phase_report,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_internal_csv(trim(arg3)//'/eos_data.csv',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_metadata(trim(arg3)//'/eos_data.nml',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(*,'(a)')'{"status":"ok","operation":"maxwell-transition","n":'//integer_text(out%n)// &
      ',"pressure_crossings":'//integer_text(phase_report%pressure_crossings)// &
      ',"low_phase_rows":'//integer_text(phase_report%low_phase_rows)// &
      ',"high_phase_rows":'//integer_text(phase_report%high_phase_rows)// &
      ',"mu_transition":'//real_text(phase_report%mu_transition)// &
      ',"pressure_transition":'//real_text(phase_report%pressure_transition)// &
      ',"n_low":'//real_text(phase_report%n_low)//',"n_high":'//real_text(phase_report%n_high)// &
      ',"epsilon_low":'//real_text(phase_report%epsilon_low)//',"epsilon_high":'//real_text(phase_report%epsilon_high)// &
      ',"density_jump":'//real_text(phase_report%density_jump)// &
      ',"latent_energy_density":'//real_text(phase_report%latent_energy_density)// &
      ',"maxwell_identity_residual":'//real_text(phase_report%maxwell_identity_residual)// &
      ',"max_euler_residual_low":'//real_text(phase_report%max_euler_residual_low)// &
      ',"max_euler_residual_high":'//real_text(phase_report%max_euler_residual_high)// &
      ',"min_cs2_low":'//real_text(phase_report%min_cs2_low)// &
      ',"max_cs2_low":'//real_text(phase_report%max_cs2_low)// &
      ',"min_cs2_high":'//real_text(phase_report%min_cs2_high)// &
      ',"max_cs2_high":'//real_text(phase_report%max_cs2_high)// &
      ',"min_cs2":'//real_text(min(phase_report%min_cs2_low,phase_report%min_cs2_high))// &
      ',"max_cs2":'//real_text(max(phase_report%max_cs2_low,phase_report%max_cs2_high))//'}'
  case('extract-core')
    if(narg/=7)then
      call emit_failure('extract-core requires SOURCE_COMPOSE OUTPUT_COMPOSE FIRST_ROW LAST_ROW EXPECTED_SOURCE_N EXPECTED_Iphase_OR_MINUS1');stop 2
    endif
    read(arg3,*,iostat=ios)core_first
    if(ios/=0)then;call emit_failure('invalid FIRST_ROW');stop 2;endif
    read(arg4,*,iostat=ios)core_last
    if(ios/=0)then;call emit_failure('invalid LAST_ROW');stop 2;endif
    read(arg5,*,iostat=ios)expected_source_n
    if(ios/=0)then;call emit_failure('invalid EXPECTED_SOURCE_N');stop 2;endif
    read(arg6,*,iostat=ios)expected_phase
    if(ios/=0)then;call emit_failure('invalid EXPECTED_Iphase');stop 2;endif
    call read_compose(trim(arg1),eos,status)
    if(status/='ok')then;call emit_failure('invalid source CompOSE table: '//trim(status));stop 3;endif
    call extract_core_rows(eos,core_first,core_last,expected_source_n,expected_phase,out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_compose(trim(arg2),out,status)
    if(status/='ok')then;call emit_failure('native core-only CompOSE write failed: '//trim(status));stop 3;endif
    call read_compose(trim(arg2),readback,status)
    if(status/='ok')then;call emit_failure('native core-only CompOSE readback failed: '//trim(status));stop 3;endif
    call verify_core_roundtrip(out,readback,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 4;endif
    write(*,'(a)')'{"status":"ok","operation":"extract-core","n":'//integer_text(out%n)// &
      ',"first_source_row":'//integer_text(core_first)//',"last_source_row":'//integer_text(core_last)// &
      ',"source_inb_start":'//integer_text(eos%inb_start)//',"output_inb_start":'//integer_text(out%inb_start)// &
      ',"expected_phase":'//integer_text(expected_phase)//',"native_readback":"PASS"}'
  case('extract-crust')
    if(narg/=7)then
      call emit_failure('extract-crust requires SOURCE_COMPOSE OUTPUT_COMPOSE EXPECTED_SOURCE_N OUTER_LAST EXPECTED_OUTER_Iphase EXPECTED_INNER_Iphase');stop 2
    endif
    read(arg3,*,iostat=ios)expected_source_n
    if(ios/=0)then;call emit_failure('invalid EXPECTED_SOURCE_N');stop 2;endif
    read(arg4,*,iostat=ios)crust_outer_last
    if(ios/=0)then;call emit_failure('invalid OUTER_LAST');stop 2;endif
    read(arg5,*,iostat=ios)expected_outer_phase
    if(ios/=0)then;call emit_failure('invalid EXPECTED_OUTER_Iphase');stop 2;endif
    read(arg6,*,iostat=ios)expected_inner_phase
    if(ios/=0)then;call emit_failure('invalid EXPECTED_INNER_Iphase');stop 2;endif
    call read_compose(trim(arg1),eos,status)
    if(status/='ok')then;call emit_failure('invalid source CompOSE table: '//trim(status));stop 3;endif
    call extract_crust_rows(eos,expected_source_n,crust_outer_last,expected_outer_phase, &
      expected_inner_phase,out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_compose(trim(arg2),out,status)
    if(status/='ok')then;call emit_failure('native crust-only CompOSE write failed: '//trim(status));stop 3;endif
    call read_compose(trim(arg2),readback,status)
    if(status/='ok')then;call emit_failure('native crust-only CompOSE readback failed: '//trim(status));stop 3;endif
    call verify_core_roundtrip(out,readback,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 4;endif
    write(*,'(a)')'{"status":"ok","operation":"extract-crust","n":'//integer_text(out%n)// &
      ',"outer_last_source_row":'//integer_text(crust_outer_last)// &
      ',"inner_first_source_row":'//integer_text(crust_outer_last+1)// &
      ',"outer_neutron_drip_nb":'//real_text(out%nb(crust_outer_last))// &
      ',"native_readback":"PASS"}'
  case default
    call emit_failure('unknown operation: '//trim(action));stop 2
  end select
contains
  subroutine emit_failure(message)
    character(len=*),intent(in)::message
    write(*,'(a)')'{"status":"error","message":"'//trim(message)//'"}'
  end subroutine
  function integer_text(i) result(s)
    integer,intent(in)::i
    character(len=32)::s
    write(s,'(i0)')i
  end function
  function real_text(x) result(s)
    real(wp),intent(in)::x
    character(len=40)::s
    write(s,'(es24.16e3)')x
  end function
  function logical_text(value) result(s)
    logical,intent(in)::value
    character(len=5)::s
    if(value)then
      s='true '
    else
      s='false'
    endif
  end function
end program compose_eos_fortran
