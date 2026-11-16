! Derived one-charge, separately neutral, T=0 two-branch research chord.
! The charge chemical potential and nuclear composition are deliberately absent.
program blocker_tf_export_chord
 use compose_eos_types,only:wp
 implicit none
 character(512)::out,path,header
 real(wp)::na,nb,mu,p,n(3),e(3),fraction(3),nr,er,pr,mur,fr
 integer::u,v,i,ios
 call get_command_argument(1,out);read(*,*)na,nb,mu,p
 if(na<=0.or.nb<=na.or.p<=0)stop 1
 n=[na,(na+nb)/2,nb];e=mu*n-p;fraction=[0._wp,.5_wp,1._wp]
 path=trim(out)//'/candidate_thermo.csv'
 open(newunit=u,file=trim(path),status='new')
 write(u,'(a)')'n_b_fm3,epsilon_MeV_fm3,pressure_MeV_fm3,mu_B_MeV,high_density_volume_fraction'
 do i=1,3
  write(u,'(5(es25.16e3,:,","))')n(i),e(i),p,mu,fraction(i)
 enddo
 close(u)
 open(newunit=v,file=trim(path),status='old');read(v,'(a)')header
 do i=1,3
  read(v,*,iostat=ios)nr,er,pr,mur,fr
  if(ios/=0)stop 2
  if(maxval(abs([nr-n(i),er-e(i),pr-p,mur-mu,fr-fraction(i)]))>1e-12_wp*maxval(abs([e(i),mu])))stop 3
  if(abs(nr*mur-er-pr)>1e-11_wp)stop 4
 enddo
 read(v,*,iostat=ios)nr
 if(ios>=0)stop 5
 close(v)
 print '(a)','PASS_NATIVE_CHORD_WRITE_READBACK_3_ROWS'
 print '(a)','RESEARCH_CANDIDATE; COMPOSITION_UNAVAILABLE; MU_Q_UNDEFINED_FOR_MACROSCOPIC_MIXTURE'
end program
