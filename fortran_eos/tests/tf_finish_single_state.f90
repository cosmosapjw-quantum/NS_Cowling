! One retained branch-state solve with unchanged native acceptance.
program tf_finish_single_state
 use compose_eos_types,only:wp
 use compose_eos_spatial_tf
 use, intrinsic :: ieee_exceptions,only:ieee_invalid,ieee_get_flag,ieee_set_flag
 implicit none
 type(tf_state_t)::s
 type(tf_grid_t)::g
 integer::n,d,i
 character(160)::model,arg,seed,geom,status
 real(wp)::nb,r,target,meanb
 real(wp),allocatable::nn(:),np(:)
 logical::invalid
 call get_command_argument(1,model);call get_command_argument(2,arg);read(arg,*)n
 call get_command_argument(3,arg);read(arg,*)target
 allocate(nn(n),np(n))
 read(*,*)d,seed,nb,r,geom
 read(*,*)nn;read(*,*)np
 if(target/=nb)then
  nn=nn*target/nb;np=np*target/nb
 endif
 call make_tf_grid(n,d,g,status)
 if(status/='ok')stop 1
 meanb=sum(g%weights*(nn+np))
 print '(a,5(1x,es25.16e3))','INITIAL',nb,target,r,meanb,meanb-target
 call ieee_set_flag(ieee_invalid,.false.)
 call solve_spatial_tf(trim(model),target,n,d,trim(seed),r,1800,1e-14_wp,s,status,initial_nn=nn,initial_np=np)
 call ieee_get_flag(ieee_invalid,invalid)
 print '(a,1x,a)','SOLVE_STATUS',trim(status)
 if(status/='ok')stop 2
 print '(a,1x,a,1x,a,9(1x,es25.16e3))','STATE',trim(s%numerical_status),trim(s%geometry), &
  s%n_b,s%epsilon,s%pressure,s%mu_B,s%max_chemical_KKT_residual_MeV,s%radius_stationarity_MeV,s%radius_fm, &
  s%baryon_constraint_error,s%charge_constraint_error
 print '(a,3(1x,i0),4(1x,l1))','OPTIMIZER',s%optimizer_mode,s%iterations,s%evaluations, &
  s%optimizer_success,s%radius_at_search_bound,s%density_at_search_bound,invalid
 do i=1,n
  print '(a,1x,i0,3(1x,es25.16e3))','PROFILE',i,s%nn(i),s%np(i),s%ne(i)
 enddo
 if(index(s%numerical_status,'STATIONARY')==0.or.trim(s%geometry)/=trim(geom))stop 3
 print '(a)','PASS_RETAINED_NATIVE_BRANCH'
end program
