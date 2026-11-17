! Retained-profile continuation adapter. Physics/acceptance stays in native module.
program tf_continuation_20261117
 use compose_eos_types,only:wp
 use compose_eos_spatial_tf
 implicit none
 type(tf_state_t)::s
 integer::n,d,i
 character(160)::model,seed,status
 real(wp)::nb,r
 real(wp),allocatable::nn(:),np(:)
 read(*,*)model,n,d,seed,nb,r
 allocate(nn(n),np(n));read(*,*)nn;read(*,*)np
 call solve_spatial_tf(trim(model),nb,n,d,trim(seed),r,1800,1e-14_wp,s,status,initial_nn=nn,initial_np=np)
 print '(a)',trim(status)
 if(status/='ok')stop 2
 print '(a,1x,a)',trim(s%numerical_status),trim(s%geometry)
 print '(*(es25.16e3,1x))',s%epsilon,s%pressure,s%mu_B,s%mu_Q_lepton, &
 s%scalar_pde_residual,s%vector_pde_residual,s%electron_pde_residual, &
 s%d_epsilon_d_log_radius,s%max_chemical_KKT_residual_MeV,s%radius_stationarity_MeV,s%Y_p,s%radius_fm, &
 real(s%iterations,wp),real(s%evaluations,wp),real(s%optimizer_mode,wp), &
 s%baryon_constraint_error,s%charge_constraint_error, &
 merge(1._wp,0._wp,s%optimizer_success),merge(1._wp,0._wp,s%radius_at_search_bound), &
 merge(1._wp,0._wp,s%density_at_search_bound),merge(1._wp,0._wp,s%muon_threshold_crossed)
 do i=1,n
 print '(*(es25.16e3,1x))',s%nn(i),s%np(i),s%ne(i),s%mu_n(i),s%mu_p(i),s%mu_e(i),s%phi(i),s%omega(i),s%rho(i),s%coulomb(i)
 enddo
end program
