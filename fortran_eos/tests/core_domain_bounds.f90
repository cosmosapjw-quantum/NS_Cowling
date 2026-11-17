program core_domain_bounds
 use compose_eos_types, only: wp
 use compose_eos_core_domain
 use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 character(len=12),parameter::models(3)=[character(len=12)::'RG_SLY4','GM_GM1','APR_APR']
 real(wp),parameter::outside(3)=[1.5_wp,1.1_wp,1.5_wp]
 real(wp)::cs,cm,fl,fr,hh,br,lo,hi,nan,inf
 character(len=240)::status
 integer::i,fail,count
 fail=0;count=0
 do i=1,3
  call core_sound_speed_fd(trim(models(i)),outside(i),1.e-4_wp,cs,cm,fl,status)
  call expect(trim(models(i))//' FD outside',status,'density outside declared native core model bounds')
  call core_sound_speed_analytic(trim(models(i)),outside(i),cs,fr,hh,br,status)
  call expect(trim(models(i))//' analytic outside',status,'density outside declared native core model bounds')
 enddo
 call core_sound_speed_fd('RG_SLY4',1.3999_wp,1.e-4_wp,cs,cm,fl,status)
 call expect('upper stencil overrun',status,'finite-difference stencil outside declared native core model bounds')
 call core_sound_speed_fd('RG_SLY4',.020001_wp,1.e-4_wp,cs,cm,fl,status)
 call expect('lower stencil overrun',status,'finite-difference stencil outside declared native core model bounds')
 call core_sound_speed_analytic('RG_SLY4',1.4_wp,cs,fr,hh,br,status)
 call expect('analytic inclusive upper endpoint',status,'ok')
 call core_sound_speed_analytic('RG_SLY4',.02_wp,cs,fr,hh,br,status)
 call expect('analytic inclusive lower endpoint',status,'ok')
 call core_sound_speed_fd('RG_SLY4',1.399_wp,2.e-5_wp,cs,cm,fl,status)
 call expect('interior full stencil',status,'ok')
 nan=ieee_value(0._wp,ieee_quiet_nan);inf=ieee_value(0._wp,ieee_positive_inf)
 call core_sound_speed_analytic('RG_SLY4',nan,cs,fr,hh,br,status)
 call expect('analytic NaN',status,'density must be finite')
 call core_sound_speed_analytic('APR_APR',inf,cs,fr,hh,br,status)
 call expect('analytic infinity',status,'density must be finite')
 call core_sound_speed_gibbs(nan,cs,status)
 call expect('Gibbs NaN',status,'density must be finite')
 call core_sound_speed_gibbs(1.5_wp,cs,status)
 call expect('Gibbs outside',status,'density outside declared native core model bounds')
 call core_sound_speed_gibbs(.22_wp,cs,status)
 call expect('Gibbs interior',status,'ok')
 call core_causal_bracket('RG_SLY4',1.5_wp,1.6_wp,1.e-4_wp,2.e-10_wp,lo,hi,status)
 call expect('bracket outside',status,'density outside declared native core model bounds')
 ! One interior derivative and one original crossing ensure guards do not alter physics.
 call core_sound_speed_fd('APR_APR',.9_wp,1.e-4_wp,cs,cm,fl,status)
 call expect('APR interior derivative',status,'ok')
 call core_sound_speed_analytic('APR_APR',.9_wp,fr,hh,br,lo,status)
 call expect('APR analytic interior',status,'ok')
 if(abs(cs-fr)>2.e-7_wp.or.abs(cs-cm)>2.e-7_wp)fail=fail+1
 call core_causal_bracket('RG_SLY4',1.24_wp,1.29_wp,1.e-4_wp,2.e-10_wp,lo,hi,status)
 call expect('original SLy4 crossing',status,'ok')
 if(lo>1.2650056715193643_wp.or.hi<1.2650056715193643_wp)fail=fail+1
 write(*,'(a,i0,a,i0)')'BOUNDARY_CHECKS=',count,' FAILURES=',fail
 if(fail/=0)stop 1
contains
 subroutine expect(label,actual,wanted)
  character(len=*),intent(in)::label,actual,wanted
  count=count+1
  write(*,'(a,2a)')trim(label),': ',trim(actual)
  if(actual/=wanted)fail=fail+1
 end subroutine
end program
