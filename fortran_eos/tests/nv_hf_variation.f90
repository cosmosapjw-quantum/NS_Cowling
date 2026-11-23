program nv_hf_variation
 use compose_eos_types,only:wp
 use compose_eos_nv_original,only:nv_original_source_t,load_nv_original_source
 use compose_eos_nv_hf,only:nv_hf_cell_energy,nv_hf_fields
 implicit none
 integer,parameter::n=128
 real(wp),parameter::pi=acos(-1._wp),radius=24._wp,h=radius/n,eps=1e-6_wp
 real(wp)::r(n),w(n),rr(n,2),dr(n),rho(n,2),tau(n,2),jp(n),up(n,2),bm(n,2),so(n),hr(n)
 real(wp)::rp(n,2),ec,ex,en,p,muec,eplus,eminus,deriv,fd,rel,edge,fac,vol,ptrial
 type(nv_original_source_t)::src
 character(1024)::path,status
 integer::i,s,trial,l(2),ncheck,mode
 call get_command_argument(1,path);call load_nv_original_source(trim(path),src,status)
 if(status/='ok')error stop 'source failed'
 do i=1,n
  r(i)=(i-.5_wp)*h;w(i)=((i*h)**3-((i-1)*h)**3)/3
 enddo
 vol=sum(4*pi*w);ncheck=0
 do mode=0,1
 do trial=1,4
  l=[mod(trial,2),mod(trial+1,2)]
  rr(:,1)=r**l(1)*exp(-r*r/(18+trial))
  rr(:,2)=r**l(2)*exp(-r*r/(15+trial))
  do s=1,2
   rr(:,s)=rr(:,s)/sqrt(sum(w*rr(:,s)**2))
  enddo
  call densities(rr,l,rho,tau,jp)
  call nv_hf_fields(src,r,w,h,rho,tau,jp,2/vol,up,bm,so,en,ec,ex,p,muec,status,mode)
  if(status/='ok')error stop 'field failed'
  do s=1,2
   dr=rr(:,s)*cos(r*(.3_wp+trial*.1_wp))
   ! Tangent to normalization, so rest masses and uniform electrons stay fixed
   ! to first order; second-order normalization is included explicitly below.
   dr=dr-rr(:,s)*sum(w*rr(:,s)*dr)
   hr=(up(:,s)+bm(:,s)*real(l(s)*(l(s)+1),wp)*h/w)*rr(:,s)
   if(s==2)hr=hr+so*real(l(s),wp)*rr(:,s)/r
   do i=1,n-1
    edge=(i*h)**2*(bm(i,s)+bm(i+1,s))/(2*h)*(rr(i,s)-rr(i+1,s))
    hr(i)=hr(i)+edge/w(i);hr(i+1)=hr(i+1)-edge/w(i+1)
   enddo
   if(mod(l(s),2)==0)hr(n)=hr(n)+2*radius**2*bm(n,s)*rr(n,s)/(h*w(n))
   deriv=4*sum(w*dr*hr)
   rp=rr;rp(:,s)=rr(:,s)+eps*dr;rp(:,s)=rp(:,s)/sqrt(sum(w*rp(:,s)**2))
   call densities(rp,l,rho,tau,jp)
   call nv_hf_cell_energy(src,r,w,h,rho,tau,jp,2,eplus,p,status,mode)
   if(status/='ok')error stop 'plus energy failed'
   rp=rr;rp(:,s)=rr(:,s)-eps*dr;rp(:,s)=rp(:,s)/sqrt(sum(w*rp(:,s)**2))
   call densities(rp,l,rho,tau,jp)
   call nv_hf_cell_energy(src,r,w,h,rho,tau,jp,2,eminus,p,status,mode)
   if(status/='ok')error stop 'minus energy failed'
   fd=(eplus-eminus)/(2*eps);rel=abs(fd-deriv)/max(1._wp,abs(deriv))
   write(*,'(a,i0,a,i0,a,i0,3(1x,es20.12))')'VARIATION mode=',mode,' trial=',trial,' species=',s,deriv,fd,rel
   if(rel>2e-7_wp)error stop 'orbital functional derivative failure'
   ncheck=ncheck+1
  enddo
  call densities(rr,l,rho,tau,jp)
  call nv_hf_cell_energy(src,r,w,h,rho,tau,jp,2,en,p,status,mode)
  if(status/='ok')error stop 'base dilation energy failed'
  fac=exp(eps)
  call nv_hf_cell_energy(src,r*fac,w*fac**3,h*fac,rho/fac**3,tau/fac**5,jp/fac**4,2,eplus,ptrial,status,mode)
  if(status/='ok')error stop 'expanded dilation energy failed'
  fac=exp(-eps)
  call nv_hf_cell_energy(src,r*fac,w*fac**3,h*fac,rho/fac**3,tau/fac**5,jp/fac**4,2,eminus,ptrial,status,mode)
  if(status/='ok')error stop 'contracted dilation energy failed'
  ptrial=-(eplus-eminus)/(6*eps*vol)
  if(abs(ptrial-p)>2e-7_wp*max(abs(p),1e-6_wp))error stop 'dilation pressure failure'
  ncheck=ncheck+1
 enddo
 enddo
 write(*,'(a,i0)')'PASS orbital_variation_checks=',ncheck
contains
 subroutine densities(orb,ell,rho,tau,jp)
  real(wp),intent(in)::orb(n,2)
  integer,intent(in)::ell(2)
  real(wp),intent(out)::rho(n,2),tau(n,2),jp(n)
  real(wp)::edge,fac
  integer::i,s
  fac=2/(4*pi);rho=fac*orb**2;tau=0
  do s=1,2
   tau(:,s)=fac*real(ell(s)*(ell(s)+1),wp)*h*orb(:,s)**2/w
   do i=1,n-1
    edge=fac*(i*h)**2*(orb(i+1,s)-orb(i,s))**2/h
    tau(i,s)=tau(i,s)+edge/(2*w(i));tau(i+1,s)=tau(i+1,s)+edge/(2*w(i+1))
   enddo
   if(mod(ell(s),2)==0)tau(n,s)=tau(n,s)+fac*2*radius**2*orb(n,s)**2/(h*w(n))
  enddo
  jp=real(ell(2),wp)*rho(:,2)/r
 end subroutine
end program
