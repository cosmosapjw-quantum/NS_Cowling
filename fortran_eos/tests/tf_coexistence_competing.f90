! Finite competition check for the TF common-tangent candidates.
! Every result is an actual native solve; this does not certify global minima.
program tf_coexistence_competing
 use compose_eos_types, only: wp
 use compose_eos_spatial_tf, only: tf_state_t, solve_tf_family_native
 implicit none
 type(tf_state_t) :: selected
 type(tf_state_t), allocatable :: attempts(:)
 integer :: n, j, k, dims(6)
 real(wp) :: nb
 character(32) :: model, arg, status
 character(16) :: seeds(6), family(6)
 call get_command_argument(1,model)
 call get_command_argument(2,arg); read(arg,*) n
 call get_command_argument(3,arg); read(arg,*) nb
 dims=[1,2,2,3,3,1]
 seeds=[character(16)::'droplet','droplet','hole','droplet','hole','uniform']
 ! Canonical labels derive from the converged profile orientation, not seed names.
 family=[character(16)::'slab','rod','tube','droplet','bubble','uniform']
 print '(a,1x,a,1x,i0,1x,es25.16e3)','CAMPAIGN',trim(model),n,nb
 do j=1,6
  call solve_tf_family_native(trim(model),nb,n,dims(j),trim(seeds(j)),selected,status,attempts)
  print '(a,1x,i0,1x,a,1x,a,1x,i0)','FAMILY',j,trim(family(j)),trim(status),size(attempts)
  do k=1,size(attempts)
   print '(a,2(1x,i0),1x,a,1x,a,8(1x,es25.16e3))','ATTEMPT',j,k,trim(attempts(k)%evaluation_status), &
    trim(attempts(k)%numerical_status),attempts(k)%n_b,attempts(k)%epsilon,attempts(k)%pressure,attempts(k)%mu_B, &
    attempts(k)%energy_per_baryon,attempts(k)%radius_fm,attempts(k)%max_chemical_KKT_residual_MeV, &
    attempts(k)%radius_stationarity_MeV
  enddo
  if(status=='ok')then
   print '(a,1x,i0,1x,a,1x,a,6(1x,es25.16e3))','SELECTED',j,trim(family(j)),trim(selected%geometry), &
    selected%epsilon,selected%pressure,selected%mu_B,selected%energy_per_baryon,selected%max_chemical_KKT_residual_MeV, &
    selected%radius_stationarity_MeV
  else
   print '(a,1x,i0,1x,a,1x,a)','NO_ACCEPTED_FAMILY',j,trim(family(j)),trim(status)
  endif
 enddo
 print '(a)','FINITE_NATIVE_FAMILY_COMPARISON_NOT_GLOBAL_MINIMUM'
end program
