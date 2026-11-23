! Local high-gradient diagnostic, not a proof of an equilibrium crust phase.
program nv_dme_gradient_stability
 use compose_eos_types,only:wp
 use compose_eos_nv_original,only:nv_original_source_t,load_nv_original_source
 use compose_eos_nv_dme,only:nv_dme_coeff_t,nv_dme_coefficients,nv_dme_h2m
 implicit none
 type(nv_original_source_t)::src
 type(nv_dme_coeff_t)::co
 real(wp),parameter::points(2,8)=reshape([ &
  .001_wp,.0001_wp,.01_wp,.001_wp,.05_wp,.01_wp,.08_wp,.02_wp, &
  .1_wp,.04_wp,.08_wp,.08_wp,.12_wp,.02_wp,.14_wp,.04_wp],[2,8])
 real(wp)::n,p,an,ap,eig,det
 character(1024)::path,status
 integer::i,mode
 call get_command_argument(1,path)
 call load_nv_original_source(trim(path),src,status)
 if(status/='ok')then;print*,trim(status);error stop 'source load';endif
 print '(a)','mode,n_n,n_p,B_n,B_p,C_n,C_p,D,Weizsaecker_plus_C_n,Weizsaecker_plus_C_p,min_eigenvalue,determinant'
 do i=1,size(points,2)
  n=points(1,i);p=points(2,i)
  do mode=0,1
   call nv_dme_coefficients(src,n,p,co,status,coefficient_mode=mode)
   if(status/='ok')then;print*,trim(status);error stop 'coefficient evaluation';endif
   an=co%c(1)+(nv_dme_h2m+co%b(1))/(4*n)
   ap=co%c(2)+(nv_dme_h2m+co%b(2))/(4*p)
   eig=.5_wp*(an+ap-sqrt((an-ap)**2+co%d**2));det=an*ap-co%d**2/4
   print '(i0,11(",",es25.16e3))',mode,n,p,co%b,co%c,co%d,an,ap,eig,det
  enddo
 enddo
end program
