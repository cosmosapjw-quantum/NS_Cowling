program core_native_probe
 use compose_eos_types, only: wp,eos_data_t
 use compose_eos_core_native, only: core_state_t,solve_core_state,generate_core_native,core_nuclear_state
 implicit none
 character(len=80)::model,arg
 character(len=240)::status
 real(wp)::n,lo,hi,e,p,mn,mp
 integer::i,count,ios,branch
 type(core_state_t)::s
 type(core_state_t),allocatable::states(:)
 real(wp),allocatable::grid(:)
 type(eos_data_t)::eos
 call get_command_argument(1,model)
 call get_command_argument(2,arg)
 if(trim(arg)=='nuclear')then
  call get_command_argument(3,arg);read(arg,*)lo
  call get_command_argument(4,arg);read(arg,*)hi
  call get_command_argument(5,arg);read(arg,*)branch
  call core_nuclear_state(trim(model),lo,hi,branch,e,p,mn,mp,status)
  if(trim(status)/='ok')then
   write(*,'(a)')trim(status)
   stop 5
  endif
  write(*,'(4(es25.16e3,1x))')e,p,mn,mp
  stop
 endif
 read(arg,*,iostat=ios)n
 if(ios/=0)stop 2
 if(command_argument_count()==2)then
  call solve_core_state(trim(model),n,s,status)
  if(trim(status)/='ok')then
   write(*,'(a)')trim(status)
   stop 3
  endif
  call emit(s)
 else
  count=n
  call get_command_argument(3,arg);read(arg,*)lo
  call get_command_argument(4,arg);read(arg,*)hi
  allocate(grid(count))
  do i=1,count
   grid(i)=lo*(hi/lo)**(real(i-1,wp)/real(count-1,wp))
  enddo
  call generate_core_native(trim(model),grid,eos,states,status)
  if(trim(status)/='ok')then
   write(*,'(a)')trim(status)
   stop 4
  endif
  do i=1,size(states)
   call emit(states(i))
  enddo
 endif
contains
 subroutine emit(s)
  type(core_state_t),intent(in)::s
  write(*,'(16(es25.16e3,1x),i0)')s%n_b,s%epsilon,s%pressure,s%mu_n,s%mu_p,s%mu_e, &
    s%yn,s%yp,s%ye,s%ymu,s%chi,s%n_low,s%n_high,s%yp_low,s%yp_high,s%gibbs_residual,s%phase
 end subroutine
end program
