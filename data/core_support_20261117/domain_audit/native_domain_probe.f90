
program core_support_audit
 use compose_eos_types,only:wp
 use compose_eos_core_native,only:core_state_t,solve_core_state
 use compose_eos_core_domain,only:core_sound_speed_analytic
 implicit none
 type(core_state_t)::s,t(5)
 character(len=32)::model
 character(len=256)::status,ast
 real(wp)::n,h,de,dp,dm,cs,csm,fl,ca,fr,hy,be,w(5)
 integer::ios,k,j,shift,attempt
 logical::same
 write(*,'(a)') 'model,n_b,epsilon,pressure,mu_n,mu_p,mu_e,Yp,Ye,Ymu,cs2_fd,cs2_mu,firstlaw_relative,fd_h,cs2_analytic,fixed_proton_fraction_lepton_equilibrated_curvature_ratio,hyy,beta_gradient,phase,fd_status,analytic_status'
 do
  read(*,*,iostat=ios) model,n
  if(ios/=0)exit
  call solve_core_state(trim(model),n,s,status)
  if(status/='ok')then
   write(0,*) model,n,trim(status)
   stop 21
  endif
  h=n*1.e-4_wp
  same=.false.
  do attempt=1,24
   shift=-2;w=[1._wp,-8._wp,0._wp,8._wp,-1._wp]
   if(n-2*h<.02_wp)then
    shift=0;w=[-25._wp,48._wp,-36._wp,16._wp,-3._wp]
   else if(n+2*h>1._wp)then
    shift=-4;w=[3._wp,-16._wp,36._wp,-48._wp,25._wp]
   endif
   do k=1,5
    j=shift+k-1
    call solve_core_state(trim(model),n+j*h,t(k),status)
    if(status/='ok')then
     write(0,*) model,n+j*h,trim(status)
     stop 22
    endif
   enddo
   same=all(t%phase==s%phase).and.all((t%ymu>0).eqv.(s%ymu>0))
   if(same)exit
   h=h/2
  enddo
  if(.not.same)then
   write(0,*) 'could not keep derivative stencil within a composition branch',model,n
   stop 23
  endif
  de=sum(w*t%epsilon)/(12*h);dp=sum(w*t%pressure)/(12*h);dm=sum(w*t%mu_n)/(12*h)
  cs=dp/de;csm=n*dm/s%mu_n;fl=abs(de/s%mu_n-1)
  call core_sound_speed_analytic(trim(model),n,ca,fr,hy,be,ast)
  write(*,'(a,",",17(es24.16e3,","),i0,",ok,",a)') trim(model),n,s%epsilon,s%pressure,s%mu_n,s%mu_p,s%mu_e,s%yp,s%ye,s%ymu,cs,csm,fl,h,ca,fr,hy,be,s%phase,trim(ast)
 enddo
end program
