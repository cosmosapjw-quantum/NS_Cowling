program native_endpoint
 use compose_eos_types,only:wp
 use compose_eos_core_native,only:core_state_t,solve_core_state
 use compose_eos_core_domain,only:core_sound_speed_analytic,core_sound_speed_fd
 implicit none
 type(core_state_t)::s
 real(wp)::n,rel,cs,frozen,hyy,beta,fd,cm,fl
 character(256)::arg,status
 call get_command_argument(1,arg);read(arg,*)n
 call get_command_argument(2,arg);read(arg,*)rel
 call solve_core_state('RG_SLY4',n,s,status)
 if(status/='ok')then;write(*,*)trim(status);stop 1;endif
 call core_sound_speed_analytic('RG_SLY4',n,cs,frozen,hyy,beta,status)
 if(status/='ok')then;write(*,*)trim(status);stop 2;endif
 call core_sound_speed_fd('RG_SLY4',n,rel,fd,cm,fl,status)
 if(status/='ok')then;write(*,*)trim(status);stop 3;endif
 write(*,'(*(es25.16e3,:,","))')n,s%epsilon,s%pressure,s%mu_n,cs,fd,cm,fl,frozen,hyy,beta
end program
