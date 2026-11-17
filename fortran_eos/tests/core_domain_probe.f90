program core_domain_probe
 use compose_eos_types, only: wp
 use compose_eos_core_native, only: core_state_t,core_model_bounds,solve_apr_gibbs_fraction
 use compose_eos_core_domain
 implicit none
 character(len=12),parameter::models(5)=[character(len=12)::'APR_APR','RG_KDE0v1','RG_SLY4','BBB_BHF2','GM_GM1']
 real(wp),parameter::steps(6)=[2.e-3_wp,1.e-3_wp,5.e-4_wp,2.e-4_wp,1.e-4_wp,5.e-5_wp]
 real(wp),parameter::lb(3)=[.85_wp,1.39_wp,1.24_wp],ub(3)=[.91_wp,1.44_wp,1.29_wp]
 real(wp)::lo,hi,mn,mp,n,cs,cm,fl,ca,fr,hh,br,cl,ch,mm,pa,minimum,maximum
 character(len=240)::status
 type(core_state_t)::s0,s1
 integer::i,j,k,unit,fail
 fail=0
 open(newunit=unit,file='core_domain_samples.csv',status='replace')
 write(unit,'(a)')'model,n_fm3,relstep,cs2_fd,cs2_mu,firstlaw_rel,cs2_analytic,cs2_frozen,hessian_yy,beta_residual'
 do i=1,5
  call core_model_bounds(trim(models(i)),lo,hi,mn,mp,status);call require(status)
  minimum=huge(1._wp);maximum=-huge(1._wp)
  do j=0,40
   n=(lo+1.e-3_wp)+(hi-lo-2.e-3_wp)*real(j,wp)/40
   call core_sound_speed_fd(trim(models(i)),n,2.e-5_wp,cs,cm,fl,status)
   if(status/='ok')then
    write(*,'(a,a,es18.9,1x,a)')'SKIP ',trim(models(i)),n,trim(status)
    cycle
   endif
   minimum=min(minimum,cs);maximum=max(maximum,cs)
   if(abs(cs-cm)>2.e-7_wp.or.fl>2.e-7_wp)fail=fail+1
   ca=-1;fr=-1;hh=-1;br=-1
   if(i<=3)then
    call core_sound_speed_analytic(trim(models(i)),n,ca,fr,hh,br,status)
    if(status=='ok')then
     if(abs(ca-cs)>2.e-7_wp.or.abs(br)>2.e-8_wp.or.fr<ca)fail=fail+1
    else
     ! APR two-charge Gibbs is tested by direct equilibrium path below.
     ca=-1;fr=-1;hh=-1;br=-1
     if(index(status,'mixed-phase')/=1)then
      write(*,'(a)')trim(status);fail=fail+1
     endif
    endif
   endif
   write(unit,'(a,9(",",es25.16e3))')trim(models(i)),n,2.e-5_wp,cs,cm,fl,ca,fr,hh,br
  enddo
  write(*,'(a,1x,a,2(es25.16e3,1x))')'SCAN',trim(models(i)),minimum,maximum
 enddo
 close(unit)
 open(newunit=unit,file='core_domain_crossings.csv',status='replace')
 write(unit,'(a)')'model,relstep,lo_fm3,hi_fm3,cs2_fd,cs2_mu,firstlaw_rel,cs2_analytic,cs2_below,cs2_above'
 do i=1,3
  do k=1,size(steps)
   call core_causal_bracket(trim(models(i)),lb(i),ub(i),steps(k),2.e-10_wp,lo,hi,status);call require(status)
   n=(lo+hi)/2
   call core_sound_speed_fd(trim(models(i)),n,steps(k),cs,cm,fl,status);call require(status)
   call core_sound_speed_analytic(trim(models(i)),n,ca,fr,hh,br,status);call require(status)
   call core_sound_speed_fd(trim(models(i)),n-1.e-5_wp,steps(k),cl,mm,pa,status);call require(status)
   call core_sound_speed_fd(trim(models(i)),n+1.e-5_wp,steps(k),ch,mm,pa,status);call require(status)
   if(abs(cs-ca)>2.e-7_wp.or.cl>=1.or.ch<=1)fail=fail+1
   write(unit,'(a,9(",",es25.16e3))')trim(models(i)),steps(k),lo,hi,cs,cm,fl,ca,cl,ch
   write(*,'(a,1x,a,4(es25.16e3,1x))')'ROOT',trim(models(i)),steps(k),lo,hi,ca
  enddo
 enddo
 close(unit)
 call solve_apr_gibbs_fraction(0._wp,s0,status);call require(status)
 call solve_apr_gibbs_fraction(1._wp,s1,status);call require(status)
 write(*,'(a,2es25.16e3)')'APR_GIBBS_ENDPOINTS ',s0%n_b,s1%n_b
 open(newunit=unit,file='core_domain_gibbs.csv',status='replace')
 write(unit,'(a)')'n_fm3,relstep,cs2_fd,cs2_mu,firstlaw_rel,cs2_gibbs_analytic'
 do j=1,3
  n=s0%n_b+(s1%n_b-s0%n_b)*j/4
  do k=3,6
   call core_sound_speed_fd('APR_APR',n,steps(k),cs,cm,fl,status);call require(status)
   call core_sound_speed_gibbs(n,ca,status);call require(status)
   if(abs(cs-cm)>1.e-6_wp.or.fl>1.e-6_wp.or.abs(cs-ca)>1.e-6_wp)fail=fail+1
   write(unit,'(6(es25.16e3,:,","))')n,steps(k),cs,cm,fl,ca
  enddo
 enddo
 close(unit)
 open(newunit=unit,file='core_domain_muon_convergence.csv',status='replace')
 write(unit,'(a)')'n_fm3,relstep,cs2_fd,cs2_analytic,absolute_difference,status'
 n=.117325_wp
 call core_sound_speed_analytic('APR_APR',n,ca,fr,hh,br,status);call require(status)
 do k=1,5
  mm=1.e-4_wp/2**(k-1)
  call core_sound_speed_fd('APR_APR',n,mm,cs,cm,fl,status)
  write(unit,'(5(es25.16e3,","),a)')n,mm,cs,ca,abs(cs-ca),trim(status)
 enddo
 if(abs(cs-ca)>2.e-9_wp)fail=fail+1
 close(unit)
 call core_sound_speed_fd('BBB_BHF2',.4_wp,1.e-4_wp,cs,cm,fl,status)
 if(index(status,'stencil crosses BBB')/=1)fail=fail+1
 call core_sound_speed_fd('APR_APR',s0%n_b,1.e-4_wp,cs,cm,fl,status)
 if(index(status,'stencil crosses')/=1)fail=fail+1
 call core_sound_speed_fd('RG_SLY4',1._wp,0._wp,cs,cm,fl,status)
 if(status=='ok')fail=fail+1
 call core_causal_bracket('RG_SLY4',.3_wp,.4_wp,1.e-4_wp,2.e-10_wp,lo,hi,status)
 if(status=='ok')fail=fail+1
 write(*,'(a,i0)')'FAILURES=',fail
 if(fail/=0)stop 1
contains
 subroutine require(status)
  character(len=*),intent(in)::status
  if(status/='ok')then
   write(*,'(a)')trim(status)
   stop 2
  endif
 end subroutine
end program
