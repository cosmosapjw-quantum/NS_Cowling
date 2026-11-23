program bbp_original_density_identity
  use, intrinsic::iso_fortran_env,only:real64
  implicit none
  integer,parameter::wp=real64
  real(wp),parameter::pi=acos(-1._wp)
  real(wp)::a,k,kg,nn,ni,ng,nb,lo,hi
  integer::ia,ik,ig,in
  lo=huge(lo);hi=-huge(hi)
  print '(a)','A,k,kg,nN,nb_from_identity'
  do ia=-1,1,2
    do ik=-1,1,2
      do ig=-1,1,2
        do in=-1,1,2
          a=990._wp+ia*.5_wp;k=1.21_wp+ik*.005_wp
          kg=1.01_wp+ig*.005_wp;nn=17.8e-6_wp+in*.05e-6_wp
          ni=k**3/(1.5_wp*pi*pi);ng=kg**3/(1.5_wp*pi*pi)
          nb=a*nn+(1-a*nn/ni)*ng
          lo=min(lo,nb);hi=max(hi,nb)
          print '(5(es25.16e3,:,","))',a,k,kg,nn,nb
        enddo
      enddo
    enddo
  enddo
  a=990;k=1.21_wp;kg=1.01_wp;nn=17.8e-6_wp
  ni=k**3/(1.5_wp*pi*pi);ng=kg**3/(1.5_wp*pi*pi)
  nb=a*nn+(1-a*nn/ni)*ng
  print '(a,es25.16e3)','nominal=',nb
  print '(a,2es25.16e3)','corner_interval=',lo,hi
  print '(a,2es25.16e3)','Table2_printed_interval=',.0789_wp-.00005_wp,.0789_wp+.00005_wp
  if(hi>=.07885_wp)error stop 'intervals unexpectedly overlap'
end program
