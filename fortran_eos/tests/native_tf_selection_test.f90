program native_tf_selection_test
 use compose_eos_types,only:wp
 use compose_eos_spatial_tf,only:tf_state_t,select_tf_candidates_native
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(tf_state_t)::s(3)
 integer::i,selected,unresolved,lower
 integer::dims(2)=[3,2]
 character(16)::seeds(2)=[character(16)::'droplet','hole']
 character(160)::status
 do i=1,3
  s(i)%n_b=.03_wp;s(i)%baryon_constraint_error=0;s(i)%charge_constraint_error=0
  s(i)%scalar_pde_residual=0;s(i)%vector_pde_residual=0;s(i)%electron_pde_residual=0
  s(i)%epsilon=30;s(i)%energy_per_baryon=1000+i
  s(i)%dimension=3;s(i)%seed='droplet';s(i)%numerical_status='STATIONARY'
 enddo
 s(2)%dimension=2;s(2)%seed='hole';s(2)%energy_per_baryon=999
 s(3)%numerical_status='REQUIRES_REFINEMENT'
 call check('LOWEST_EVALUATED_STATIONARY_CANDIDATE',2,0,0)
 s(3)%energy_per_baryon=998
 call check('LOWER_UNCONVERGED_ATTEMPT',0,0,1)
 s(3)%energy_per_baryon=1003;s(2)%seed='droplet'
 call check('HOLD_UNRESOLVED_GEOMETRY',0,1,0)
 s(2)%seed='hole';s(2)%energy_per_baryon=ieee_value(1._wp,ieee_quiet_nan)
 call check('HOLD_UNRESOLVED_GEOMETRY',0,1,0)
 print '(a)','4 native finite-family selection checks PASS'
contains
 subroutine check(expected,index,wanted_unresolved,wanted_lower)
  character(*),intent(in)::expected
  integer,intent(in)::index,wanted_unresolved,wanted_lower
  call select_tf_candidates_native(s,dims,seeds,.03_wp,selected,status,unresolved,lower)
  if(status/=expected.or.selected/=index.or.unresolved/=wanted_unresolved.or.lower/=wanted_lower)then
   print *,trim(status),selected,unresolved,lower
   stop 1
  endif
 end subroutine
end program
