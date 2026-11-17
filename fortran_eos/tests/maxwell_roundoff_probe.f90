program maxwell_roundoff_probe
 use compose_eos_types, only: wp,eos_data_t,eos_pressure,eos_energy
 use compose_eos_io, only: read_internal_csv
 use compose_eos_phase_transition, only: construct_maxwell_transition,maxwell_transition_report_t
 implicit none
 type(eos_data_t)::low,high,out
 type(maxwell_transition_report_t)::report
 character(len=2048)::low_path,high_path
 character(len=512)::status
 real(wp),allocatable::p(:),e(:)
 integer::i
 call get_command_argument(1,low_path)
 call get_command_argument(2,high_path)
 call read_internal_csv(trim(low_path)//'/eos_data.csv',trim(low_path)//'/eos_data.nml',low,status)
 if(status/='ok')stop 2
 call read_internal_csv(trim(high_path)//'/eos_data.csv',trim(high_path)//'/eos_data.nml',high,status)
 if(status/='ok')stop 2
 call construct_maxwell_transition(low,high,.true.,out,report,status)
 write(*,'(a)')'CONSTRUCTION_STATUS='//trim(status)
 if(.not.allocated(out%nb))stop 3
 p=eos_pressure(out);e=eos_energy(out)
 write(*,'(a,i0)')'OUTPUT_ROWS=',out%n
 write(*,'(a,i0)')'NONINCREASING_DENSITIES=',count(out%nb(2:)<=out%nb(:out%n-1))
 write(*,'(a,i0)')'NONINCREASING_ENERGIES=',count(e(2:)<=e(:out%n-1))
 write(*,'(a,es24.16e3)')'MIN_WITHIN_PHASE_DELTA_P=', &
   minval(p(2:)-p(:out%n-1),mask=out%iphase(2:)==out%iphase(:out%n-1))
 do i=1,out%n-1
  if(out%iphase(i)==out%iphase(i+1).and.p(i+1)>=p(i))cycle
  write(*,'(a,2i8)')'INTERVAL_ZERO_BASED=',i-1,i
  write(*,'(a,2i8)')'PHASES=',out%iphase(i),out%iphase(i+1)
  write(*,'(a,2es24.16e3)')'DENSITIES=',out%nb(i),out%nb(i+1)
  write(*,'(a,2es24.16e3)')'Q1=',out%q(1,i),out%q(1,i+1)
  write(*,'(a,2es24.16e3)')'PRESSURES=',p(i),p(i+1)
  write(*,'(a,es24.16e3)')'DELTA_P=',p(i+1)-p(i)
  write(*,'(a,es24.16e3)')'PRESSURE_ULP=',spacing(max(p(i),p(i+1)))
  write(*,'(a,es24.16e3)')'DELTA_P_ULPS=',(p(i+1)-p(i))/spacing(max(p(i),p(i+1)))
 enddo
end program

