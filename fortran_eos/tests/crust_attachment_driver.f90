program crust_attachment_driver
 use compose_eos_types,only:wp
 use compose_eos_lower_crust_native,only:lower_crust_state_t
 use compose_eos_crust_attachment
 use, intrinsic :: ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(lower_crust_state_t)::s
 type(crust_integral_budget_t)::b
 type(crust_attachment_t)::d
 character(len=256)::arg,status,action
 real(wp)::n,n0,aerr,r,tol,left(4),right(4),boundary
 logical::ext
 call get_command_argument(1,action)
 if(action=='seam')then
  left=[.002_wp,.002_wp*941._wp-.002_wp,.002_wp,941._wp]
  right=[.003_wp,.003_wp*941.3_wp-.0025_wp,.0025_wp,941.3_wp]
  call source_attachment_diagnostics(left,right,.002_wp,d,status)
  if(status/='ok'.or..not.d%source_prescription_valid)error stop 'source seam rejected'
  if(d%equilibrium_established.or.d%cross_interval_supported)error stop 'invented phase equivalence'
  if(source_attachment_supports(d,.0025_wp))error stop 'invented gap interpolation'
  if(.not.source_attachment_supports(d,.002_wp))error stop 'lost endpoint'
  write(*,'(4(es25.16e3,1x))')d%delta_pressure,d%delta_mu,d%left_euler_residual,d%right_euler_residual
  right(1)=.001_wp
  call source_attachment_diagnostics(left,right,.002_wp,d,status)
  if(status/='ok'.or..not.d%density_overlap)error stop 'source overlap not preserved'
  if(d%cross_interval_supported)error stop 'source overlap made interpolable'
  call source_attachment_diagnostics(left,right,.004_wp,d,status)
  if(status=='ok'.or.d%source_prescription_valid)error stop 'unbracketed boundary accepted'
  right(1)=.003_wp;left(2)=ieee_value(0._wp,ieee_quiet_nan)
  call source_attachment_diagnostics(left,right,.002_wp,d,status)
  if(status=='ok'.or.d%source_prescription_valid)error stop 'NaN endpoint accepted'
  stop
 endif
 call get_command_argument(2,arg);read(arg,*)n
 call get_command_argument(3,arg);read(arg,*)n0
 call get_command_argument(4,arg);read(arg,*)aerr
 call get_command_argument(5,arg);read(arg,*)r
 call get_command_argument(6,arg);read(arg,*)tol
 call get_command_argument(7,arg);ext=trim(arg)=='extended'
 call bsk22_reference_integral(n0,n,aerr,r,s,b,status,tol,ext)
 if(status/='ok')then;write(*,'(a)')trim(status);stop 2;endif
 write(*,'(21(es25.16e3,1x))')s%nb,s%epsilon,s%pressure,s%mu_b, &
 s%source_epsilon,s%source_pressure,s%source_thermo_residual,s%quadrature_error, &
 s%gamma,s%sound2,s%composition_mu_residual,b%anchor_e,b%pressure_integral, &
 b%energy_per_baryon_model_bound,b%mu_model_bound,b%epsilon_model_bound, &
 b%epsilon_quadrature_estimate,b%mu_quadrature_estimate,b%source_delta_energy_per_baryon, &
 s%yn+s%yp+s%yh-1,s%ye-s%yh*s%zav/s%aav
end program
