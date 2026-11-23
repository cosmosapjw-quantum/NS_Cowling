program density_probe
 use compose_eos_types,only:wp
 use compose_eos_bbp_original,only:bbp_original_state_t,bbp_original_solve
 implicit none
 type(bbp_original_state_t)::st
 character(160)::status
 real(wp),parameter::rho_conv=1.7826619216278976e12_wp
 real(wp)::n(4)=[.076_wp,.0769_wp,.077_wp,.0789_wp],lo,hi,mid
 integer::i
 print '(a)','control,nb,rho,P,A,Z,x,k,kg,u'
 do i=1,4
  call bbp_original_solve(n(i),.true.,st,status)
  if(status/='ok')then;print *,trim(status);stop 2;endif
  print '(a,9(es25.16e3,:,","))','fixed_nb,',st%nb,st%energy*rho_conv,st%pressure,st%a,st%z,st%x,st%k,st%kg,st%u
 enddo
 lo=.075_wp;hi=.0789_wp
 do i=1,36
  mid=(lo+hi)/2
  call bbp_original_solve(mid,.true.,st,status)
  if(status/='ok')stop 2
  if(st%energy*rho_conv<1.30e14_wp)then;lo=mid;else;hi=mid;endif
 enddo
 call bbp_original_solve((lo+hi)/2,.true.,st,status)
 if(status/='ok')stop 2
 print '(a,9(es25.16e3,:,","))','fixed_rho,',st%nb,st%energy*rho_conv,st%pressure,st%a,st%z,st%x,st%k,st%kg,st%u
end program
