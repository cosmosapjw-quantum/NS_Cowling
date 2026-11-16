! Exploratory, two-neutral-branch T=0 common-tangent solve.
! Every branch value is an actual native TF minimization; no table extrapolation.
program blocker_tf_common_tangent
 use compose_eos_types,only:wp
 use compose_eos_spatial_tf
 implicit none
 type(tf_state_t)::s(2),tmp(2),minus,plus
 type(tf_grid_t)::grid(2)
 integer::n,d(2),j,k,it,seq=0
 character(160)::model,arg,status,seed(2),geom(2)
 real(wp),allocatable::nn(:,:),np(:,:)
 real(wp)::nb(2),r(2),mu(2),p(2),curv(2),fmu,fp,da,db,step,dx(2),h,lower,upper
 real(wp)::startnb(2),startR(2),target,mu_t,p_t,chord,energygap(2),mixgapmin
 logical::converged
 call get_command_argument(1,model);call get_command_argument(2,arg);read(arg,*)n
 allocate(nn(n,2),np(n,2))
 do j=1,2
  read(*,*)d(j),seed(j),nb(j),r(j),geom(j)
  read(*,*)nn(:,j);read(*,*)np(:,j)
  call make_tf_grid(n,d(j),grid(j),status)
  call native_solve(j,nb(j),r(j),nn(:,j),np(:,j),s(j))
 enddo
 startnb=nb;startR=r;lower=minval(nb)-.003_wp;upper=maxval(nb)+.003_wp
 h=2e-5_wp;converged=.false.
 print '(a,2(1x,es25.16e3))','BOUNDED_ACTUAL_EVALUATION_DOMAIN',lower,upper
 do it=1,16
  do j=1,2
   call local_branch_derivative(j,s(j),h,mu(j),p(j),curv(j))
  enddo
  fmu=mu(1)-mu(2);fp=p(1)-p(2)
  print '(a,1x,i0,10(1x,es25.16e3))','ITERATION',it,nb,mu,p,curv,fmu,fp
  if(abs(fmu)<2e-6_wp.and.abs(fp)<2e-7_wp)then
   converged=.true.;exit
  endif
  if(any(curv<=0).or.nb(2)<=nb(1))then
   print '(a)','HOLD_NONPOSITIVE_CURVATURE_OR_ORDER';stop 4
  endif
  da=(-fp+nb(2)*fmu)/((nb(1)-nb(2))*curv(1))
  db=(fmu+curv(1)*da)/curv(2);dx=[da,db]
  step=min(1._wp,.001_wp/max(maxval(abs(dx)),tiny(1._wp)))
  do while(any(nb+step*dx<=lower).or.any(nb+step*dx>=upper).or.nb(1)+step*da>=nb(2)+step*db)
   step=step/2
   if(step<1e-4_wp)then
    print '(a)','HOLD_OUTSIDE_BOUNDED_DOMAIN';stop 5
   endif
  enddo
  do j=1,2
   target=nb(j)+step*dx(j)
   call native_solve(j,target,s(j)%radius_fm,s(j)%nn*target/nb(j),s(j)%np*target/nb(j),tmp(j))
  enddo
  s=tmp;nb=nb+step*dx
 enddo
 if(.not.converged)then
  print '(a)','HOLD_ROOT_NOT_CONVERGED';stop 6
 endif
 mu_t=sum(mu)/2;p_t=sum(p)/2
 print '(a,8(1x,es25.16e3))','COEXISTENCE',nb,mu_t,p_t,fmu,fp,s(2)%epsilon-s(1)%epsilon-mu_t*(nb(2)-nb(1)), &
  (s(2)%epsilon/nb(2)-s(1)%epsilon/nb(1))
 do j=1,2
  call local_branch_derivative(j,s(j),h/2,mu(j),p(j),curv(j))
 enddo
 print '(a,6(1x,es25.16e3))','HALF_STEP_CHECK',mu,p,mu(1)-mu(2),p(1)-p(2)
 mixgapmin=huge(1._wp)
 do k=0,4
  target=nb(1)+(nb(2)-nb(1))*real(k,wp)/4
  chord=mu_t*target-p_t
  do j=1,2
   call native_solve(j,target,s(j)%radius_fm,s(j)%nn*target/nb(j),s(j)%np*target/nb(j),tmp(j))
   energygap(j)=tmp(j)%epsilon-chord
  enddo
  if(k>0.and.k<4)mixgapmin=min(mixgapmin,minval(energygap))
  print '(a,4(1x,es25.16e3))','CHORD_CHECK',target,energygap,minval(energygap)
 enddo
 print '(a,1x,es25.16e3)','MIN_INTERIOR_SAMPLED_CHORD_GAP',mixgapmin
 do j=1,2
  do k=1,n
   print '(a,2(1x,i0),3(1x,es25.16e3))','PROFILE',j,k,s(j)%nn(k),s(j)%np(k),s(j)%ne(k)
  enddo
 enddo
 print '(a)','CONVERGED_TWO_BRANCH_CANDIDATE_NOT_GLOBAL_PHASE_APPROVAL'
 contains
 subroutine native_solve(branch,density,radius,inn,inp,out)
 integer,intent(in)::branch
 real(wp),intent(in)::density,radius,inn(:),inp(:)
 type(tf_state_t),intent(out)::out
 call solve_spatial_tf(trim(model),density,n,d(branch),trim(seed(branch)),radius,1800,1e-13_wp,out,status,initial_nn=inn,initial_np=inp)
 seq=seq+1
 print '(a,2(1x,i0),1x,a)','SOLVE_STATUS',seq,branch,trim(status)
 if(status/='ok')stop 2
 print '(a,2(1x,i0),1x,a,1x,a,9(1x,es25.16e3))','STATE',seq,branch,trim(out%numerical_status),trim(out%geometry), &
  out%n_b,out%epsilon,out%pressure,out%mu_B,out%max_chemical_KKT_residual_MeV,out%radius_stationarity_MeV,out%radius_fm, &
  out%baryon_constraint_error,out%charge_constraint_error
 if(index(out%numerical_status,'STATIONARY')==0.or.trim(out%geometry)/=trim(geom(branch)))then
  print '(a)','HOLD_BRANCH_CHANGE_OR_UNRESOLVED';stop 3
 endif
 end subroutine
 subroutine local_branch_derivative(branch,base,delta,mub,press,curvature)
 integer,intent(in)::branch
 type(tf_state_t),intent(in)::base
 real(wp),intent(in)::delta
 real(wp),intent(out)::mub,press,curvature
 real(wp)::n0
 n0=base%n_b
 call native_solve(branch,n0-delta,base%radius_fm,base%nn*(n0-delta)/n0,base%np*(n0-delta)/n0,minus)
 call native_solve(branch,n0+delta,base%radius_fm,base%nn*(n0+delta)/n0,base%np*(n0+delta)/n0,plus)
 mub=(plus%epsilon-minus%epsilon)/(2*delta);press=n0*mub-base%epsilon
 curvature=(plus%epsilon-2*base%epsilon+minus%epsilon)/delta**2
 end subroutine
end program
