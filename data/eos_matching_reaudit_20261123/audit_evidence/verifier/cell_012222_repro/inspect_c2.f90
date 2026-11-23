program inspect_c2
  use compose_eos_types
  use compose_eos_io, only: read_internal_csv
  implicit none
  type(eos_data_t)::e
  character(len=2048)::directory
  character(len=512)::status
  real(wp),allocatable::mu_source(:),mu_euler(:),relative(:)
  real(wp)::pressure,energy
  integer::i,j
  call get_command_argument(1,directory)
  call read_internal_csv(trim(directory)//'/eos_data.csv',trim(directory)//'/eos_data.nml',e,status)
  if(trim(status)/='ok')then
    print '(a)',trim(status)
    stop 1
  endif
  allocate(mu_source(e%n),mu_euler(e%n),relative(e%n))
  mu_source=(e%q(3,:)+1._wp)*e%mn
  if(e%il==0)mu_source=mu_source+e%charge_fraction_grid*e%q(4,:)*e%mn
  mu_euler=(e%q(6,:)+1._wp)*e%mn+e%q(1,:)
  relative=abs(mu_source-mu_euler)/max(abs(mu_source),abs(mu_euler),tiny(1._wp))
  j=maxloc(relative,1)
  print '(a)',trim(directory)
  print '(a,i0,a,i0,a,es25.16e3)','n=',e%n,' worst_row_1based=',j,' max_relative=',relative(j)
  print '(a)','row_1based,n_fm-3,P_MeV_fm-3,epsilon_MeV_fm-3,mu_source_MeV,mu_Euler_MeV,delta_mu_MeV,relative'
  do i=1,e%n
    if(relative(i)<=2.e-6_wp.and.i/=j)cycle
    pressure=e%nb(i)*e%q(1,i)
    energy=e%nb(i)*(1._wp+e%q(7,i))*e%mn
    print '(i0,7(a,es25.16e3))',i,',',e%nb(i),',',pressure,',',energy,',',mu_source(i),',',mu_euler(i),',',mu_source(i)-mu_euler(i),',',relative(i)
  enddo
end program inspect_c2
