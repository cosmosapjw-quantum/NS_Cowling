program nv_hf_probe
 use compose_eos_types,only:wp
 use compose_eos_nv_original,only:nv_original_source_t,load_nv_original_source
 use compose_eos_nv_hf
 implicit none
 type(nv_original_source_t)::src
 type(nv_hf_config_t)::cfg
 type(nv_hf_result_t)::res
 character(1024)::source,output,arg,status
 integer::i,u,ios,c0,c1,rate
 if(command_argument_count()<2)then
  write(*,'(a)')'usage: nv_hf_probe SOURCE OUTPUT [N Z nb mesh maxiter mixing seed smoothing coefficient_mode]';stop 2
 endif
 call get_command_argument(1,source);call get_command_argument(2,output)
 cfg%log_unit=6
 do i=3,command_argument_count()
  call get_command_argument(i,arg)
  select case(i)
  case(3);read(arg,*,iostat=ios)cfg%neutrons
  case(4);read(arg,*,iostat=ios)cfg%protons
  case(5);read(arg,*,iostat=ios)cfg%nb
  case(6);read(arg,*,iostat=ios)cfg%mesh
  case(7);read(arg,*,iostat=ios)cfg%max_iterations
  case(8);read(arg,*,iostat=ios)cfg%mixing
  case(9);read(arg,*,iostat=ios)cfg%seed_kind
  case(10);read(arg,*,iostat=ios)cfg%outer_smoothing
  case(11);read(arg,*,iostat=ios)cfg%coefficient_mode
  case default;stop 2
  end select
  if(ios/=0)stop 2
 enddo
 call load_nv_original_source(trim(source),src,status)
 if(status/='ok')then;write(*,'(a)')trim(status);stop 2;endif
 call system_clock(c0,rate)
 call nv_hf_solve(src,cfg,res,status)
 call system_clock(c1)
 open(newunit=u,file=trim(output),status='replace',action='write',iostat=ios)
 if(ios/=0)stop 2
 write(u,'(a)')'# NV1973 orbital DME HF fixed-composition calculation'
 write(u,'(a)')'status='//trim(status)
 write(u,'(a,l1)')'converged=',res%converged
 write(u,'(a,l1)')'energy_valid=',res%energy_valid
 write(u,'(a,l1)')'spectral_complete=',res%spectral_complete
 write(u,'(a,l1)')'outer_smoothing=',res%outer_smoothing
 write(u,'(a,i0)')'neutrons=',cfg%neutrons
 write(u,'(a,i0)')'protons=',cfg%protons
 write(u,'(a,i0)')'mesh=',cfg%mesh
 write(u,'(a,i0)')'seed_kind=',cfg%seed_kind
 write(u,'(a,i0)')'coefficient_mode=',cfg%coefficient_mode
 write(u,'(a,i0)')'iterations=',res%iterations
 write(u,'(a,i0)')'lmax=',res%lmax_used
 write(u,'(a,i0)')'radial_states=',res%radial_states_used
 call value('nb_fm3',cfg%nb);call value('radius_fm',res%radius)
 call value('energy_per_baryon_MeV',res%energy_per_baryon)
 call value('pressure_MeV_fm3',res%pressure)
 call value('mu_n_kin_MeV',res%mu_n);call value('mu_p_kin_MeV',res%mu_p)
 call value('mu_e_total_MeV',res%mu_e);call value('beta_residual_MeV',res%beta_residual)
 call value('mu_n_addition_kin_MeV',res%mu_n_plus);call value('mu_p_addition_kin_MeV',res%mu_p_plus)
 call value('beta_interval_lower_MeV',res%beta_interval_lower)
 call value('beta_interval_upper_MeV',res%beta_interval_upper)
 call value('density_residual',res%density_residual);call value('energy_residual_MeV',res%energy_residual)
 call value('tau_residual',res%tau_residual);call value('spin_residual',res%spin_residual)
 call value('selfconsistent_residual',res%selfconsistent_residual)
 call value('maximum_raw_neutron_density_fm3',res%maximum_raw_density(1))
 call value('maximum_raw_proton_density_fm3',res%maximum_raw_density(2))
 call value('last_admissible_energy_per_baryon_MeV',res%last_admissible_energy_per_baryon)
 call value('particle_error',res%particle_error);call value('dilation_pressure_error_MeV_fm3',res%dilation_pressure_error)
 call value('nuclear_energy_MeV',res%nuclear_energy);call value('coulomb_direct_MeV',res%coulomb_direct)
 call value('coulomb_exchange_MeV',res%coulomb_exchange);call value('electron_energy_MeV',res%electron_energy)
 call value('boundary_neutron_density_fm3',res%boundary_neutron_density)
 call value('minimum_effective_mass_coefficient',res%minimum_effective_mass_coefficient)
 call value('wall_seconds',real(c1-c0,wp)/rate)
 if(allocated(res%r))then
  write(u,'(a)')'# r_fm weight_fm3 rho_n rho_p tau_n tau_p J_p'
  do i=1,size(res%r)
   write(u,'(7(es25.16e3,1x))')res%r(i),res%weight(i),res%rho(i,:),res%tau(i,:),res%jp(i)
  enddo
 endif
 close(u)
 write(*,'(a,1x,l1,1x,i0,5(1x,es14.6))')trim(status),res%converged,res%iterations,res%energy_per_baryon,res%pressure, &
 res%density_residual,res%energy_residual,real(c1-c0,wp)/rate
 if(status/='ok'.or..not.res%converged)stop 3
contains
 subroutine value(name,x)
  character(*),intent(in)::name
  real(wp),intent(in)::x
  write(u,'(a,a,es25.16e3)')name,'=',x
 end subroutine
end program
