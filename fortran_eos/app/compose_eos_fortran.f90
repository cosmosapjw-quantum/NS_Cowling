program compose_eos_fortran
  use compose_eos_types
  use compose_eos_io
  use compose_eos_verify
  implicit none
  type(eos_data_t)::eos
  type(verify_report_t)::report
  character(len=32)::action
  character(len=2048)::arg1,arg2
  character(len=512)::status
  integer::narg,ios
  action='';arg1='';arg2=''
  narg=command_argument_count()
  if(narg<2)then
    write(*,'(a)')'usage: compose_eos_fortran {read|verify|verify-csv|normalize|write} INPUT OUTPUT'
    stop 2
  endif
  call get_command_argument(1,action,status=ios)
  call get_command_argument(2,arg1,status=ios)
  if(narg>=3)call get_command_argument(3,arg2,status=ios)
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
end program compose_eos_fortran
