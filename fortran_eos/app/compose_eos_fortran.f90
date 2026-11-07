program compose_eos_fortran
  use compose_eos_types
  use compose_eos_io
  use compose_eos_verify
  use compose_eos_reconstruct, only: fill_density_gap
  use compose_eos_causality, only: causal_complete, causal_completion_report_t
  implicit none
  type(eos_data_t)::eos,out
  type(verify_report_t)::report
  type(causal_completion_report_t)::causal_report
  character(len=32)::action
  character(len=2048)::arg1,arg2
  character(len=128)::arg3,arg4,arg5
  character(len=512)::status
  integer::narg,ios,n_fill,fit_points,gap_phase,n_points
  real(wp)::cs2_extension,n_max_requested
  logical::has_gap_phase
  action='';arg1='';arg2='';arg3='';arg4='';arg5=''
  narg=command_argument_count()
  if(narg<2)then
    write(*,'(a)')'usage: compose_eos_fortran {read|verify|verify-csv|normalize|write|fill-gap|causal-complete} INPUT OUTPUT [OPTIONS]'
    stop 2
  endif
  call get_command_argument(1,action,status=ios)
  call get_command_argument(2,arg1,status=ios)
  if(narg>=3)call get_command_argument(3,arg2,status=ios)
  if(narg>=4)call get_command_argument(4,arg3,status=ios)
  if(narg>=5)call get_command_argument(5,arg4,status=ios)
  if(narg>=6)call get_command_argument(6,arg5,status=ios)
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
    if(narg/=5.and.narg/=6)then
      call emit_failure('fill-gap requires input internal EOS directory, output directory, n_fill, fit_points, and optional gap_phase');stop 2
    endif
    read(arg3,*,iostat=ios)n_fill
    if(ios/=0)then;call emit_failure('invalid n_fill');stop 2;endif
    read(arg4,*,iostat=ios)fit_points
    if(ios/=0)then;call emit_failure('invalid fit_points');stop 2;endif
    has_gap_phase=narg==6;gap_phase=0
    if(has_gap_phase)then
      read(arg5,*,iostat=ios)gap_phase
      if(ios/=0)then;call emit_failure('invalid gap_phase');stop 2;endif
    endif
    call read_internal_csv(trim(arg1)//'/eos_data.csv',trim(arg1)//'/eos_data.nml',eos,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call fill_density_gap(eos,n_fill,fit_points,gap_phase,has_gap_phase,out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_internal_csv(trim(arg2)//'/eos_data.csv',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    call write_metadata(trim(arg2)//'/eos_data.nml',out,status)
    if(status/='ok')then;call emit_failure(trim(status));stop 3;endif
    write(*,'(a)')'{"status":"ok","operation":"fill-gap","n":'//integer_text(out%n)//'}'
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
