program nv_hf_numerics
 use compose_eos_types,only:wp
 use compose_eos_nv_hf,only:nv_hf_radial_spectrum,nv_hf_gradient,nv_hf_gradient_adjoint,nv_hf_electrons
 implicit none
 real(wp),parameter::pi=acos(-1._wp),radius=20._wp,mass=20.721248538623252_wp
 real(wp),parameter::roots(3,0:3)=reshape([ &
  pi,2*pi,3*pi, &
  2.0815759778181006105_wp,5.9403699905727124972_wp,9.2058401429366649154_wp, &
  5.7634591968945497914_wp,9.0950113304763551563_wp,12.322940970566582052_wp, &
  4.5140996470322816772_wp,8.5837549563657667415_wp,11.972730032192526309_wp],[3,4])
 real(wp),allocatable::r(:),w(:),b(:),u(:),ss(:),rv(:,:),x(:),q(:),gx(:),aq(:)
 real(wp)::ev(3),exact(3),err(2),h,norm,maxerr,adjerror,ee,pe,mu,ne,eplus,eminus,dmu
 integer::mesh,n,l,i,k,ncheck
 character(1024)::status
 ncheck=0
 do mesh=1,2
  n=120*2**mesh
  allocate(r(n),w(n),b(n),u(n),ss(n),rv(n,3),x(n),q(n),gx(n),aq(n))
  h=radius/n
  do i=1,n
   r(i)=(i-.5_wp)*h;w(i)=((i*h)**3-((i-1)*h)**3)/3
  enddo
  b=mass;u=0;ss=0;maxerr=0
  do l=0,3
   call nv_hf_radial_spectrum(r,w,h,b,u,l,0._wp,ss,3,ev,rv,status)
   if(status/='ok')error stop 'radial eigensolver failed'
   exact=mass*(roots(:,l)/radius)**2
   maxerr=max(maxerr,maxval(abs(ev/exact-1)))
   do k=1,3
    norm=sum(w*rv(:,k)**2)
    if(abs(norm-1)>2e-12_wp)error stop 'radial normalization failure'
    ncheck=ncheck+1
   enddo
   write(*,'(a,i0,a,i0,3(1x,es18.10))')'FREE mesh=',n,' l=',l,ev
  enddo
  err(mesh)=maxerr
  x=exp(-(r/7)**2);q=r*exp(-(r/9)**2)
  call nv_hf_gradient(x,h,gx);call nv_hf_gradient_adjoint(q,w,h,aq)
  adjerror=abs(sum(w*gx*q)-sum(w*x*aq))
  if(adjerror>2e-12_wp)error stop 'density derivative adjoint failure'
  ncheck=ncheck+1
  deallocate(r,w,b,u,ss,rv,x,q,gx,aq)
 enddo
 if(err(2)>1e-4_wp.or.err(2)/err(1)>.27_wp)error stop 'radial second-order convergence failure'
 ncheck=ncheck+1
 do i=1,8
  ne=10._wp**(-12+i)
  call nv_hf_electrons(ne,ee,pe,mu)
  call nv_hf_electrons(ne*(1+1e-5_wp),eplus,dmu,norm)
  call nv_hf_electrons(ne*(1-1e-5_wp),eminus,dmu,norm)
  if(abs((eplus-eminus)/(2e-5_wp*ne)-mu)>1e-7_wp*max(mu,1._wp))error stop 'electron chemical potential derivative'
  if(abs(ne*mu-ee-pe)>1e-14_wp)error stop 'electron Euler identity'
  ncheck=ncheck+2
 enddo
 write(*,'(a,i0,3(a,es18.10))')'PASS checks=',ncheck,' free_error_coarse=',err(1), &
 ' free_error_fine=',err(2),' error_ratio=',err(2)/err(1)
end program
