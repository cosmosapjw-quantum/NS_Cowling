program bbp_original_benchmark
  use compose_eos_types,only:wp
  use compose_eos_bbp_original,only:bbp_original_state_t,bbp_original_solve
  implicit none
  type(bbp_original_state_t)::published,variational
  real(wp),parameter::n(9)=[2.79e-4_wp,5.25e-4_wp,8.79e-4_wp,1.59e-3_wp,3.73e-3_wp, &
    8.91e-3_wp,.0204_wp,.0475_wp,.0789_wp]
  real(wp),parameter::kr(9)=[1.32_wp,1.32_wp,1.32_wp,1.31_wp,1.30_wp,1.27_wp,1.24_wp,1.20_wp,1.21_wp]
  real(wp),parameter::xr(9)=[.313_wp,.307_wp,.302_wp,.295_wp,.280_wp,.256_wp,.220_wp,.163_wp,.121_wp]
  real(wp),parameter::ar(9)=[127._wp,134._wp,140._wp,149._wp,170._wp,211._wp,296._wp,548._wp,990._wp]
  real(wp),parameter::zr(9)=[40._wp,41._wp,42._wp,44._wp,48._wp,54._wp,65._wp,89._wp,120._wp]
  real(wp),parameter::pr(9)=[5e-4_wp,6.42e-4_wp,8.73e-4_wp,1.47e-3_wp,4.04e-3_wp, &
    1.38e-2_wp,.0477_wp,.176_wp,.384_wp]
  real(wp),parameter::rhor(9)=[4.66e11_wp,8.79e11_wp,1.47e12_wp,2.67e12_wp,6.25e12_wp, &
    1.50e13_wp,3.44e13_wp,8.01e13_wp,1.30e14_wp]
  real(wp),parameter::rho_per_energy=1.602176634e33_wp/(2.99792458e10_wp)**2
  integer::i
  character(160)::status
  print '(a)','nb_ref,k_ref,x_ref,A_ref,Z_ref,P_ref,rho_ref,k_published,x_published,A_published,'// &
    'Z_published,P_published,rho_published,P_variational,relative_P_published_minus_ref'
  do i=1,9
    call bbp_original_solve(n(i),.true.,published,status)
    if(trim(status)/='ok')error stop 'published benchmark failed'
    call bbp_original_solve(n(i),.false.,variational,status)
    if(trim(status)/='ok')error stop 'variational benchmark failed'
    print '(15(es24.16e3,:,","))',n(i),kr(i),xr(i),ar(i),zr(i),pr(i),rhor(i), &
      published%k,published%x,published%a,published%z,published%pressure, &
      published%energy*rho_per_energy,variational%pressure,published%pressure/pr(i)-1
  enddo
end program
