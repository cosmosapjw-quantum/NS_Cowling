program nv_hf_basis_extension
 use compose_eos_types,only:wp
 use compose_eos_nv_original,only:nv_original_source_t,load_nv_original_source
 use compose_eos_nv_hf
 implicit none
 type(nv_original_source_t)::src
 type(nv_hf_config_t)::cfg
 type(nv_hf_result_t)::base,expanded
 character(1024)::source,status
 real(wp)::energy_difference,pressure_difference
 call get_command_argument(1,source);call load_nv_original_source(trim(source),src,status)
 if(status/='ok')error stop 'source failed'
 cfg%neutrons=130;cfg%protons=50;cfg%nb=2.79e-4_wp;cfg%mesh=160
 cfg%coefficient_mode=1;cfg%mixing=.08_wp;cfg%max_iterations=800
 call nv_hf_solve(src,cfg,base,status)
 if(status/='ok'.or..not.base%converged)error stop 'base SCF failed'
 cfg%lmax=base%lmax_used+2;cfg%radial_states=base%radial_states_used+2
 call nv_hf_solve(src,cfg,expanded,status,base)
 if(status/='ok'.or..not.expanded%converged)error stop 'expanded SCF failed'
 energy_difference=abs(expanded%energy_per_baryon-base%energy_per_baryon)
 pressure_difference=abs(expanded%pressure-base%pressure)
 if(energy_difference>1e-8_wp.or.pressure_difference>2e-8_wp)error stop 'basis extension changed result'
 write(*,'(a,2(i0,1x),a,2(i0,1x))')'BASIS base=',base%lmax_used,base%radial_states_used, &
 ' expanded=',expanded%lmax_used,expanded%radial_states_used
 write(*,'(a,es25.16)')'energy_difference_MeV_per_baryon=',energy_difference
 write(*,'(a,es25.16)')'pressure_difference_MeV_fm3=',pressure_difference
 write(*,'(a,es25.16)')'base_energy_MeV_per_baryon=',base%energy_per_baryon
 write(*,'(a,es25.16)')'expanded_energy_MeV_per_baryon=',expanded%energy_per_baryon
 write(*,'(a,es25.16)')'expanded_selfconsistent_residual=',expanded%selfconsistent_residual
 write(*,'(a)')'PASS occupied-state truncation extension at stated cell'
end program
