! Finite native branch continuation from independently retained stationary seeds.
! Root parameters are reused only after separate root and h/2 evidence checks.
program tf_finish_chord
 use compose_eos_types,only:wp
 use compose_eos_spatial_tf
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use,intrinsic::ieee_exceptions,only:ieee_invalid,ieee_get_flag,ieee_set_flag
 implicit none
 integer::n,d(2),j,k,seq=0,u
 character(160)::model,arg,status,seed(2),geom(2),file
 real(wp)::na,nb,mu_t,p_t,seed_n(2),r(2),target,chord,gaps(2),gapmin
 real(wp),allocatable::nn(:,:),np(:,:)
 type(tf_state_t)::s
 logical::invalid
 call get_command_argument(1,model);call get_command_argument(2,arg);read(arg,*)n
 allocate(nn(n,2),np(n,2))
 read(*,*)na,nb,mu_t,p_t
 if(.not.all(ieee_is_finite([na,nb,mu_t,p_t])).or.na<=0.or.nb<=na.or.p_t<=0)stop 1
 do j=1,2
  read(*,*)d(j),seed(j),seed_n(j),r(j),geom(j)
  read(*,*)nn(:,j);read(*,*)np(:,j)
 enddo
 print '(a,4(1x,es25.16e3))','REUSED_NATIVE_ROOT',na,nb,mu_t,p_t
 gapmin=huge(1._wp)
 do k=0,4
  target=na+(nb-na)*real(k,wp)/4
  chord=mu_t*target-p_t
  do j=1,2
   print '(a,2(1x,i0),3(1x,es25.16e3))','WARM_SEED',k,j,seed_n(j),r(j),target
   call ieee_set_flag(ieee_invalid,.false.)
   call solve_spatial_tf(trim(model),target,n,d(j),trim(seed(j)),r(j),1800,1e-14_wp,s,status, &
     initial_nn=nn(:,j)*target/seed_n(j),initial_np=np(:,j)*target/seed_n(j))
   call ieee_get_flag(ieee_invalid,invalid)
   seq=seq+1
   print '(a,2(1x,i0),1x,a)','SOLVE_STATUS',seq,j,trim(status)
   if(status/='ok')stop 2
   print '(a,2(1x,i0),1x,a,1x,a,9(1x,es25.16e3))','STATE',seq,j,trim(s%numerical_status),trim(s%geometry), &
     s%n_b,s%epsilon,s%pressure,s%mu_B,s%max_chemical_KKT_residual_MeV,s%radius_stationarity_MeV,s%radius_fm, &
     s%baryon_constraint_error,s%charge_constraint_error
   print '(a,4(1x,i0),4(1x,l1))','OPTIMIZER',seq,s%optimizer_mode,s%iterations,s%evaluations, &
     s%optimizer_success,s%radius_at_search_bound,s%density_at_search_bound,invalid
   write(file,'(a,i3.3,a)')'state_',seq,'.txt'
   open(newunit=u,file=trim(file),status='new',action='write')
   write(u,*)d(j),trim(seed(j)),s%n_b,s%radius_fm,trim(s%geometry)
   write(u,*)s%nn;write(u,*)s%np;close(u)
   flush(6)
   if(s%numerical_status/='STATIONARY'.or.trim(s%geometry)/=trim(geom(j)))then
    print '(a)','HOLD_BRANCH_CHANGE_OR_UNRESOLVED';stop 3
   endif
   gaps(j)=s%epsilon-chord
   ! Native continuation uses the actual previous converged branch profile.
   nn(:,j)=s%nn;np(:,j)=s%np;r(j)=s%radius_fm;seed_n(j)=target
  enddo
  print '(a,4(1x,es25.16e3))','CHORD_CHECK',target,gaps,minval(gaps)
  if(k>0.and.k<4)gapmin=min(gapmin,minval(gaps))
 enddo
 print '(a,1x,es25.16e3)','MIN_INTERIOR_SAMPLED_CHORD_GAP',gapmin
 if(gapmin<=0)stop 4
 print '(a)','CONVERGED_FINITE_CHORD_NOT_GLOBAL_PHASE_APPROVAL'
end program
