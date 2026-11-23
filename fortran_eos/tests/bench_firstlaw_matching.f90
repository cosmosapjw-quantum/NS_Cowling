program bench_firstlaw_matching
  use iso_fortran_env, only:int64
  use compose_eos_types
  use compose_eos_io, only:read_internal_csv
  use compose_eos_matching, only:match_crust_core, match_report_t
  implicit none
  character(400)::directory,buffer,status,output
  type(eos_data_t)::cr,co,matched
  type(match_report_t)::report
  integer::repeats,query_size,rep,unit
  integer(int64)::start,stop,rate
  if(command_argument_count()/=5)error stop 'usage: bench CRUST_DIR CORE_DIR REPETITIONS BRIDGE_NODES OUTPUT_DUMP'
  call get_command_argument(1,directory)
  call read_internal_csv(trim(directory)//'/eos_data.csv',trim(directory)//'/eos_data.nml',cr,status)
  if(status/='ok')error stop 'cannot read crust'
  call get_command_argument(2,directory)
  call read_internal_csv(trim(directory)//'/eos_data.csv',trim(directory)//'/eos_data.nml',co,status)
  if(status/='ok')error stop 'cannot read core'
  call get_command_argument(3,buffer);read(buffer,*)repeats
  call get_command_argument(4,buffer);read(buffer,*)query_size
  call get_command_argument(5,output)
  do rep=0,repeats
    call system_clock(start,rate)
    call match_crust_core(cr,co,'first_law_constrained_hermite_pressure_bridge',query_size,.true.,matched,report,status)
    call system_clock(stop)
    if(status/='ok')then;print *,trim(status);error stop 'no match';endif
    if(rep>0)print '(es20.12)',real(stop-start,wp)/real(rate,wp)
  enddo
  open(newunit=unit,file=output,form='unformatted',access='stream',status='replace')
  write(unit)report,matched%nb,matched%q,matched%iphase
  close(unit)
end program
