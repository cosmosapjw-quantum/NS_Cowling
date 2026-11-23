program bbp_rounding_feasibility
 use iso_fortran_env,only: real64
 implicit none
 integer,parameter::wp=real64,k=26
 real(wp),parameter::c=29979245800._wp,conv=1e-39_wp/1.602176634e-6_wp
 real(wp)::rho(k),p(k),n(k),er(k),nr(k),eps(k),inte(k),lnfac(k),mu(k),g,h,s,a,q,lo,hi,plo,phi,de_min,de_max
 integer::i,u,ios,z,aa
 character(1024)::path,line
 character(40)::rt,pt,nt
 call get_command_argument(1,path)
 open(newunit=u,file=trim(path),status='old')
 i=0
 do
 read(u,'(a)',iostat=ios)line
 if(ios/=0)exit
 if(index(adjustl(line),'#')==1.or.len_trim(line)==0)cycle
 i=i+1
 read(line,*)rt,pt,nt,z,aa,g
 read(rt,*)rho(i);read(pt,*)p(i);read(nt,*)n(i)
 er(i)=half_decimal(rt)*c*c*conv;nr(i)=half_decimal(nt)*1e-39_wp
 eps(i)=rho(i)*c*c*conv;p(i)=p(i)*conv;n(i)=n(i)*1e-39_wp
 if(i==k)exit
 enddo
 close(u)
 inte=0;lnfac=0;mu=(eps+p)/n
 do i=2,k
  h=log(n(i)/n(i-1));g=log(p(i)/p(i-1))/h;s=g-1
  inte(i)=inte(i-1)+(p(i-1)/n(i-1))*h*exprel(s*h)
  h=log(eps(i)/eps(i-1));g=log(p(i)/p(i-1))/h;s=g-1;a=p(i-1)/eps(i-1)
  q=a*expm1_stable(s*h)/(1+a)
  if(abs(s)>1e-8_wp)then
   lnfac(i)=lnfac(i-1)+h-log1p_stable(q)/s
  else
   lnfac(i)=lnfac(i-1)+h/(1+a)
  endif
 enddo
 lo=maxval((eps-er)/n-inte);hi=minval((eps+er)/n-inte)
 print '(a,3es25.16)','fixed n, rho rounding interval energy anchor: ',lo,hi,hi-lo
 lo=maxval((n-nr)*exp(-lnfac));hi=minval((n+nr)*exp(-lnfac))
 print '(a,3es25.16)','fixed rho, n rounding interval density anchor: ',lo,hi,hi-lo
 print '(a,i0)','raw negative dmu steps: ',count(mu(2:)<=mu(:k-1))
 print '(a)','simultaneous rounding necessary cold first-law intervals: row, de_min,de_max, integral_min,integral_max'
 do i=2,k
  de_min=(eps(i)-er(i))/(n(i)+nr(i))-(eps(i-1)+er(i-1))/(n(i-1)-nr(i-1))
  de_max=(eps(i)+er(i))/(n(i)-nr(i))-(eps(i-1)-er(i-1))/(n(i-1)+nr(i-1))
  plo=(p(i-1)-.0005_wp*10._wp**floor(log10(p(i-1)/conv))*conv)* &
   (1/(n(i-1)+nr(i-1))-1/(n(i)-nr(i)))
  phi=(p(i)+.0005_wp*10._wp**floor(log10(p(i)/conv))*conv)* &
   (1/(n(i-1)-nr(i-1))-1/(n(i)+nr(i)))
  if(de_min>phi.or.de_max<plo)print '(i0,4(",",es25.16))',i-1,de_min,de_max,plo,phi
 enddo
 print '(a)','row, epsilon/n-integral, energy_half_rounding, n_interval_anchor_low, n_interval_anchor_high'
 do i=1,k
  print '(i0,4(",",es25.16))',i,eps(i)/n(i)-inte(i),er(i)/n(i),(n(i)-nr(i))*exp(-lnfac(i)), &
   (n(i)+nr(i))*exp(-lnfac(i))
 enddo
contains
 real(wp)function half_decimal(t)result(v)
 character(*),intent(in)::t
 integer::j,e,pos
 pos=index(t,'E');read(t(pos+1:),*)e;j=pos-index(t,'.')-1
 v=.5_wp*10._wp**(e-j)
 end function
 real(wp)function exprel(t)result(v)
 real(wp),intent(in)::t
 if(abs(t)<1e-5_wp)then
 v=1+t*(.5_wp+t*(1._wp/6+t*(1._wp/24+t/120)))
 else
 v=(exp(t)-1)/t
 endif
 end function
 real(wp)function expm1_stable(t)result(v)
 real(wp),intent(in)::t
 v=t*exprel(t)
 end function
 real(wp)function log1p_stable(t)result(v)
 real(wp),intent(in)::t
 if(abs(t)<1e-5_wp)then
 v=t*(1-t*(.5_wp-t*(1._wp/3-t*(.25_wp-t*.2_wp))))
 else
 v=log(1+t)
 endif
 end function
end program
