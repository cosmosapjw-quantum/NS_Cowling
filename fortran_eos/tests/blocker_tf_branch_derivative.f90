! Local native branch continuation and independent optimized-energy differences.
! No pressure fitting, phase mixing or density removal.
program blocker_tf_branch_derivative
 use compose_eos_types, only:wp
 use compose_eos_spatial_tf
 implicit none
 type(tf_state_t)::center,local(4)
 type(tf_grid_t)::g
 character(160)::model,arg,seed,status
 real(wp),allocatable::nn(:),np(:)
 real(wp)::nb,radius,e0,p0,mu0,h,offset(4),target,shape,rdiff,d1,d2,p1,p2
 integer::n,d,i,j
 call get_command_argument(1,model)
 call get_command_argument(2,arg);read(arg,*)n
 call get_command_argument(3,arg);read(arg,*)d
 call get_command_argument(4,arg);read(arg,*)nb
 call get_command_argument(5,arg);read(arg,*)radius
 call get_command_argument(6,seed)
 allocate(nn(n),np(n))
 read(*,*)e0,p0,mu0
 read(*,*)nn
 read(*,*)np
 call make_tf_grid(n,d,g,status)
 if(status/='ok')stop 1
 call solve_spatial_tf(trim(model),nb,n,d,trim(seed),radius,1800,1e-13_wp,center,status,initial_nn=nn,initial_np=np)
 print '(a,1x,a)','CENTER_STATUS',trim(status)
 if(status/='ok')stop 2
 call print_state('CENTER',center,0._wp,0._wp,0._wp)
 print '(a,3(1x,es25.16e3))','REFINEMENT_DELTA',center%epsilon-e0,center%pressure-p0,center%mu_B-mu0
 h=min(1e-5_wp,1e-3_wp*nb);offset=[-h,h,-h/2,h/2]
 do j=1,4
  target=nb+offset(j)
  call solve_spatial_tf(trim(model),target,n,d,trim(seed),center%radius_fm,1800,1e-13_wp,local(j),status, &
    initial_nn=center%nn*target/nb,initial_np=center%np*target/nb)
  print '(a,1x,i0,1x,a)','LOCAL_STATUS',j,trim(status)
  if(status/='ok')stop 3
  shape=sqrt(sum(g%weights*((local(j)%nn*nb/target-center%nn)**2+(local(j)%np*nb/target-center%np)**2)))/nb
  rdiff=local(j)%radius_fm/center%radius_fm-1
  call print_state('LOCAL',local(j),offset(j),shape,rdiff)
 enddo
 d1=(local(2)%epsilon-local(1)%epsilon)/(2*h)
 d2=(local(4)%epsilon-local(3)%epsilon)/h
 p1=nb*d1-center%epsilon;p2=nb*d2-center%epsilon
 print '(a,9(1x,es25.16e3))','DERIVATIVE',h,d1,d2,p1,p2,center%pressure,p1-center%pressure,p2-center%pressure,p2-p1
 ! Profiles bind the actual branch, not only a phase or geometry label.
 do j=0,4
  do i=1,n
   if(j==0)then
    print '(a,2(1x,i0),3(1x,es25.16e3))','PROFILE',j,i,center%nn(i),center%np(i),center%ne(i)
   else
    print '(a,2(1x,i0),3(1x,es25.16e3))','PROFILE',j,i,local(j)%nn(i),local(j)%np(i),local(j)%ne(i)
   endif
  enddo
 enddo
 contains
 subroutine print_state(label,s,delta,dshape,dradius)
 character(*),intent(in)::label
 type(tf_state_t),intent(in)::s
 real(wp),intent(in)::delta,dshape,dradius
 print '(a,1x,a,1x,a,13(1x,es25.16e3),3(1x,i0))',label,trim(s%numerical_status),trim(s%geometry),s%n_b,delta,s%epsilon,s%pressure,s%mu_B, &
  s%max_chemical_KKT_residual_MeV,s%radius_stationarity_MeV,s%radius_fm,dshape,dradius,s%scalar_pde_residual,s%vector_pde_residual,s%electron_pde_residual, &
  s%iterations,s%evaluations,s%optimizer_mode
 end subroutine
end program
