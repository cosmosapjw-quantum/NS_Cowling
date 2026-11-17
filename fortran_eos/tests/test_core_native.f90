program test_core_native
 use compose_eos_types, only: wp,eos_data_t
 use compose_eos_core_native, only: core_state_t,solve_core_state,generate_core_native,solve_apr_gibbs_fraction, &
   core_nuclear_state,core_free_fermion
 implicit none
 character(len=16),parameter::models(5)=[character(len=16)::'APR_APR','BBB_BHF2','GM_GM1','RG_KDE0v1','RG_SLY4']
 type(core_state_t)::s,sp,sm
 type(core_state_t),allocatable::states(:)
 type(eos_data_t)::eos
 character(len=240)::status
 integer::i,j,branch
 real(wp)::n,h,deriv,grid(4),nn,np,e,p,mn,mp,ep,em,dummy1,dummy2,dummy3,x,expected
 do i=1,5
  do j=1,4
   n=.12_wp+.137_wp*real(j,wp);h=n*1.e-5_wp
   call solve_core_state(trim(models(i)),n,s,status)
   call require(trim(status)=='ok','state: '//trim(status))
   call require(abs(s%yn+s%yp-1)<1.e-12_wp,'baryon closure')
   call require(abs(s%yp-s%ye-s%ymu)<1.e-10_wp,'charge closure')
   call require(abs(s%epsilon+s%pressure-n*s%mu_n)<1.e-8_wp,'Euler identity')
   call require(abs(s%mu_n-s%mu_p-s%mu_e)<1.e-8_wp,'beta equilibrium')
   call solve_core_state(trim(models(i)),n+h,sp,status);call require(trim(status)=='ok','plus state')
   call solve_core_state(trim(models(i)),n-h,sm,status);call require(trim(status)=='ok','minus state')
   deriv=(sp%epsilon-sm%epsilon)/(2*h)
   call require(abs(deriv-s%mu_n)<2.e-5_wp,'independent first law')
  enddo
  ! Derivatives at fixed species density independently test chemical potentials,
  ! beyond identities imposed algebraically by the beta-equilibrium solver.
  nn=.37_wp*.88_wp;np=.37_wp*.12_wp;h=1.e-6_wp
  do branch=1,merge(2,1,i==1)
   call core_nuclear_state(trim(models(i)),nn,np,branch,e,p,mn,mp,status)
   call require(trim(status)=='ok','asymmetric nuclear state')
   call core_nuclear_state(trim(models(i)),nn+h,np,branch,ep,dummy1,dummy2,dummy3,status)
   call require(trim(status)=='ok','neutron derivative plus')
   call core_nuclear_state(trim(models(i)),nn-h,np,branch,em,dummy1,dummy2,dummy3,status)
   call require(trim(status)=='ok','neutron derivative minus')
   call require(abs((ep-em)/(2*h)-mn)<5.e-6_wp,'independent neutron chemical derivative')
   call core_nuclear_state(trim(models(i)),nn,np+h,branch,ep,dummy1,dummy2,dummy3,status)
   call require(trim(status)=='ok','proton derivative plus')
   call core_nuclear_state(trim(models(i)),nn,np-h,branch,em,dummy1,dummy2,dummy3,status)
   call require(trim(status)=='ok','proton derivative minus')
   call require(abs((ep-em)/(2*h)-mp)<5.e-6_wp,'independent proton chemical derivative')
  enddo
  write(*,'(a)')'PASS independent core '//trim(models(i))
 enddo
 do j=0,4
  call solve_apr_gibbs_fraction(real(j,wp)/4,s,status)
  call require(trim(status)=='ok','APR Gibbs')
  call require(s%gibbs_residual<1.e-7_wp,'Gibbs equations')
  call require(s%n_high>s%n_low,'ordered phases')
  call require(abs(s%yp-s%ye-s%ymu)<1.e-10_wp,'Gibbs charge')
 enddo
 ! Nonrelativistic limit checks cancellation-free Fermi-gas pressure.
 x=1.e-4_wp;n=(x*939._wp/197.3269804_wp)**3/(3*acos(-1._wp)**2)
 call core_free_fermion(n,939._wp,e,p,mn)
 expected=n*939._wp*x*x/5
 call require(abs(p/expected-1)<1.e-8_wp,'nonrelativistic Fermi pressure')
 call require(abs((e+p)/(n*mn)-1)<1.e-14_wp,'Fermi-gas Euler identity')
 ! Source endpoints are inside support; a genuinely outside point is rejected.
 call solve_core_state('BBB_BHF2',.08_wp,s,status)
 call require(trim(status)=='ok','exact BBB lower endpoint')
 call solve_core_state('BBB_BHF2',1.5_wp,s,status)
 call require(trim(status)=='ok','exact BBB upper endpoint')
 call solve_core_state('BBB_BHF2',.079999_wp,s,status)
 call require(trim(status)/='ok','BBB outside source support rejected')
 grid=[.5_wp,.53_wp,.57_wp,.6_wp]
 call generate_core_native('APR_APR',grid,eos,states,status)
 call require(trim(status)=='ok','APR generation')
 call require(eos%n==4.and.all(eos%iphase==7),'no outside APR endpoints')
 call require(all(eos%c3_applicable).and.all(eos%c4_applicable),'native composition masks')
 call solve_core_state('NO_SUCH_MODEL',.16_wp,s,status)
 call require(trim(status)/='ok','unknown model rejected')
 grid=[.08_wp,.16_wp,.12_wp,.4_wp]
 call generate_core_native('GM_GM1',grid,eos,states,status)
 call require(trim(status)/='ok','unordered grid rejected')
 write(*,'(a)')'PASS native core functional, Gibbs, composition and generation tests'
contains
 subroutine require(condition,label)
  logical,intent(in)::condition
  character(len=*),intent(in)::label
  if(.not.condition)then
   write(*,'(a)')'FAIL '//trim(label)
   stop 1
  endif
 end subroutine
end program
