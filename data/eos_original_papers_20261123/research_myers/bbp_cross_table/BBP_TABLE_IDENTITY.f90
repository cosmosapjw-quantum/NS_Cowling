! Independent of all EOS production modules.
! BBP1971 Eq.(2.5), p.230 and Table1 definitions, p.251.
program bbp_table_identity
 use, intrinsic::iso_fortran_env,only:real64
 implicit none
 integer,parameter::wp=real64
 real(wp),parameter::pi=acos(-1._wp),rho_conv=1.7826619216278976e12_wp
 real(wp)::amin,amax,kmin,kmax,kgmin,kgmax,nnmin,nnmax,lo,hi,wide_lo,wide_hi
 real(wp)::nominal,printed_lo,printed_hi,gap,awide
 integer::ia,ik,ig,ic,pass
 real(wp)::a,k,kg,nn,ni,ng,u,nb
 amin=989.5_wp;amax=990.5_wp
 kmin=1.205_wp;kmax=1.215_wp;kgmin=1.005_wp;kgmax=1.015_wp
 nnmin=17.75e-6_wp;nnmax=17.85e-6_wp
 nominal=nb_from_table(990._wp,1.21_wp,1.01_wp,17.8e-6_wp)
 printed_lo=.07885_wp;printed_hi=.07895_wp
 do pass=1,2
  if(pass==2)then;amin=985._wp;amax=995._wp;endif
  lo=huge(1._wp);hi=-huge(1._wp)
  do ia=0,1
   a=amin+(amax-amin)*ia
   do ik=0,1
    k=kmin+(kmax-kmin)*ik
    do ig=0,1
     kg=kgmin+(kgmax-kgmin)*ig
     do ic=0,1
      nn=nnmin+(nnmax-nnmin)*ic
      ni=k**3/(1.5_wp*pi*pi);ng=kg**3/(1.5_wp*pi*pi);u=a*nn/ni
      ! All partial derivatives of nb(a,k,kg,nn) are positive in this box.
      if(ni<=ng.or.u<=0.or.u>=1)error stop 'monotonicity domain failed'
      nb=a*nn+(1-u)*ng
      lo=min(lo,nb);hi=max(hi,nb)
     enddo
    enddo
   enddo
  enddo
  if(pass==1)then
   wide_lo=lo;wide_hi=hi
  else
   ! Even A rounded to the nearest ten does not make the printed nb feasible.
   awide=hi
  endif
 enddo
 gap=printed_lo-awide
 if(gap<=1e-12_wp)error stop 'claimed table inconsistency not established'
 print '(a)','{'
 print '(a,es25.16e3,a)',' "Table1_nominal_nb":',nominal,','
 print '(a,es25.16e3,a)',' "Table1_half_last_digit_nb_min":',wide_lo,','
 print '(a,es25.16e3,a)',' "Table1_half_last_digit_nb_max":',wide_hi,','
 print '(a,es25.16e3,a)',' "Table1_A_plusminus5_nb_min":',lo,','
 print '(a,es25.16e3,a)',' "Table1_A_plusminus5_nb_max":',hi,','
 print '(a,es25.16e3,a)',' "Table2_printed_nb_min":',printed_lo,','
 print '(a,es25.16e3,a)',' "Table2_printed_nb_max":',printed_hi,','
 print '(a,es25.16e3,a)',' "minimum_disjoint_gap_fm3":',gap,','
 print '(a,es25.16e3,a)',' "Table2_nominal_energy_per_baryon_MeV":',1.30e14_wp/(.0789_wp*rho_conv),','
 print '(a)',' "original_nb_retained":0.0789,'
 print '(a)',' "possible_typo_nb_inference_only":0.0769,'
 print '(a)',' "rounding_boxes_disjoint":true,'
 print '(a)',' "production_module_used":false}'
contains
 pure real(wp) function nb_from_table(a,k,kg,nn) result(nb)
  real(wp),intent(in)::a,k,kg,nn
  real(wp)::ni,ng
  ni=k**3/(1.5_wp*pi*pi);ng=kg**3/(1.5_wp*pi*pi)
  nb=a*nn+(1-a*nn/ni)*ng
 end function
end program
