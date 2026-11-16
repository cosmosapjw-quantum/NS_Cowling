program native_tf_driver
  use compose_eos_types, only: wp
  use compose_eos_spatial_tf
  implicit none
  type(tf_grid_t) :: g
  type(tf_state_t) :: s
  type(tf_state_t),allocatable :: attempts(:)
  character(128) :: arg,model,seed,status,mode
  real(wp),allocatable :: nn(:),np(:),ne(:),v(:)
  real(wp) :: radius,nbar,mu,error,ftol
  integer :: n,d,i,maxiter
  call get_command_argument(1,mode)
  if(trim(mode)=='electron') then
    call get_command_argument(2,arg);read(arg,*) n
    call get_command_argument(3,arg);read(arg,*) d
    call get_command_argument(4,arg);read(arg,*) radius
    call make_tf_grid(n,d,g,status)
    if(trim(status)/='ok') then
      print '(a)',trim(status);stop 1
    endif
    allocate(np(n),ne(n),v(n));read(*,*) np
    call screened_electrons_native(g,np,radius,ne,v,mu,error,status)
    print '(a)',trim(status)
    if(trim(status)/='ok') stop 1
    print '(*(es25.16e3,1x))',mu,error
    do i=1,n
      print '(*(es25.16e3,1x))',ne(i),v(i)
    enddo
  else
    call get_command_argument(2,model)
    call get_command_argument(3,arg);read(arg,*) n
    call get_command_argument(4,arg);read(arg,*) d
    call get_command_argument(5,arg);read(arg,*) nbar
    call get_command_argument(6,arg);read(arg,*) radius
    call make_tf_grid(n,d,g,status)
    if(trim(status)/='ok') then
      print '(a)',trim(status);stop 1
    endif
    if(trim(mode)=='evaluate') then
      allocate(nn(n),np(n),ne(n));read(*,*) nn;read(*,*) np;read(*,*) ne
      call evaluate_spatial_tf_profile(trim(model),nbar,g,nn,np,ne,radius,s,status)
    else if(trim(mode)=='family') then
      call get_command_argument(7,seed)
      call solve_tf_family_native(trim(model),nbar,n,d,trim(seed),s,status,attempts)
      do i=1,size(attempts)
        write(0,'(a,i0,1x,a,1x,a,1x,4(es18.10e3,1x))') 'ATTEMPT ',i,trim(attempts(i)%numerical_status), &
          trim(attempts(i)%evaluation_status),attempts(i)%energy_per_baryon,attempts(i)%max_chemical_KKT_residual_MeV, &
          attempts(i)%radius_stationarity_MeV,attempts(i)%radius_fm
      enddo
    else if(trim(mode)=='solve') then
      call get_command_argument(7,seed)
      call get_command_argument(8,arg);read(arg,*) maxiter
      call get_command_argument(9,arg);read(arg,*) ftol
      call solve_spatial_tf(trim(model),nbar,n,d,trim(seed),radius,maxiter,ftol,s,status)
    else
      print '(a)','bad_driver_mode';stop 1
    endif
    print '(a)',trim(status)
    if(trim(status)/='ok') stop 1
    print '(*(es25.16e3,1x))',s%epsilon,s%pressure,s%mu_B,s%mu_Q_lepton, &
      s%scalar_pde_residual,s%vector_pde_residual,s%electron_pde_residual, &
      s%d_epsilon_d_log_radius,s%max_chemical_KKT_residual_MeV, &
      s%radius_stationarity_MeV,s%Y_p,s%radius_fm, &
      real(s%iterations,wp),real(s%evaluations,wp),real(s%optimizer_mode,wp), &
      merge(1._wp,0._wp,index(s%numerical_status,'STATIONARY')>0)
    do i=1,n
      print '(*(es25.16e3,1x))',s%nn(i),s%np(i),s%ne(i),s%mu_n(i),s%mu_p(i), &
        s%mu_e(i),s%phi(i),s%omega(i),s%rho(i),s%coulomb(i)
    enddo
  endif
end program native_tf_driver
