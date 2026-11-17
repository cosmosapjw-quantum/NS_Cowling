! Independent explicit CLD energy checks for the conditional generated SPG cells.
! Nuclear bulk and lepton routines are shared; derivatives, surface, Coulomb and
! constrained-coordinate mapping are evaluated independently of spg_native_cell.
program spg_finish_validation
 use compose_eos_types,only:wp
 use compose_eos_gdfm_native
 use compose_eos_spg_native
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 real(wp),parameter::pi=acos(-1._wp),e2=1.43996448_wp,me=.51099895_wp
 real(wp),parameter::steps(5)=[1e-3_wp,3e-4_wp,1e-4_wp,3e-5_wp,1e-5_wp],scales(4)=[1._wp,.1_wp,.1_wp,1._wp]
 type(gdfm_native_model_t)::model
 type(spg_native_result_t)::state,near(4)
 character(32)::model_id,mode
 character(160)::status
 integer::row,critical,edge,ios,j,k,l
 real(wp)::n,sigma0,b,z(4),energy_ref,p_ref,mu_ref,e,g(4),ind,der(4),trial(4),h,ep,em,epp,emm
 real(wp)::pfd,mufd,cs2,pd,ed,offset(4),values(4),pg(4),maxk,ebar0
 call get_command_argument(1,mode)
 if(trim(mode)=='thermo')then
  write(*,'(a)')'model,row,n_b,edge,relative_step,mu_reoptimized_fd,mu_native,mu_error,pressure_reoptimized_fd,pressure_native,pressure_error,cs2_over_c2,max_neighbor_KKT'
 else
  write(*,'(a)')'model,row,n_b,relative_step,epsilon_independent_error,pressure_fd_error,mu_fd_error,gradient_t_error,gradient_ni_error,gradient_y_error,gradient_v_error,projected_gradient,gas_proton_margin,muon_margin'
 endif
 do
  read(*,*,iostat=ios)model_id,row,critical,edge,n,sigma0,b,z,energy_ref,p_ref,mu_ref
  if(ios<0)exit
  if(ios/=0)stop 2
  call load_gdfm_native(trim(model_id),model,status);call require_ok('load')
  call spg_native_cell(model,sigma0,b,n,z,e,g,state,status);call require_ok('cell')
  ind=explicit_energy(n,z)
  pg=g
  if(z(1)>=1-1e-8_wp)pg(1)=max(g(1),0._wp)
  if(z(4)<=1e-8_wp)pg(4)=min(g(4),0._wp)
  if(trim(mode)=='thermo')then
   if(critical==0)cycle
   do j=1,size(steps)
    h=n*steps(j)
    if(edge==1)then
     offset=[1._wp,2._wp,3._wp,4._wp]
    elseif(edge==-1)then
     offset=[-1._wp,-2._wp,-3._wp,-4._wp]
    else
     offset=[-2._wp,-1._wp,1._wp,2._wp]
    endif
    maxk=0
    do l=1,4
     call spg_native_solve(model,sigma0,b,n+offset(l)*h,near(l),status,start=z)
     call require_ok('thermo neighbor')
     if(.not.near(l)%stationary)stop 5
     maxk=max(maxk,near(l)%max_projected_gradient)
    enddo
    values=near%epsilon
    ed=derivative(energy_ref,values,h,edge)
    values=near%pressure
    pd=derivative(p_ref,values,h,edge)
    values=near%energy_per_baryon
    pfd=n*n*derivative(state%energy_per_baryon,values,h,edge)
    mufd=ed;cs2=pd/ed
    write(*,'(a,",",i0,",",es25.16e3,",",i0,",",*(es25.16e3,:,","))') &
     trim(model_id),row,n,edge,steps(j),mufd,mu_ref,mufd-mu_ref,pfd,p_ref,pfd-p_ref,cs2,maxk
   enddo
   cycle
  endif
  do j=1,size(steps)
   h=n*steps(j)
   ep=explicit_energy(n+h,z);em=explicit_energy(n-h,z)
   epp=explicit_energy(n+2*h,z);emm=explicit_energy(n-2*h,z)
   mufd=(-epp+8*ep-8*em+emm)/(12*h)
   pfd=n*n*(-epp/(n+2*h)+8*ep/(n+h)-8*em/(n-h)+emm/(n-2*h))/(12*h)
   der=0
   if(critical/=0)then
    do k=1,4
     h=steps(j)*scales(k)
     trial=z;trial(k)=z(k)+h
     if(k==1.and.z(k)+2*h>1)then
      trial(k)=z(k)-h;em=explicit_energy(n,trial)
      trial(k)=z(k)-2*h;emm=explicit_energy(n,trial)
      der(k)=(3*ind-4*em+emm)/(2*h*n)
     elseif(k==4.and.z(k)<2*h)then
      ep=explicit_energy(n,trial);trial(k)=z(k)+2*h;epp=explicit_energy(n,trial)
      der(k)=(-3*ind+4*ep-epp)/(2*h*n)
     else
      ep=explicit_energy(n,trial);trial(k)=z(k)-h;em=explicit_energy(n,trial)
      trial(k)=z(k)+2*h;epp=explicit_energy(n,trial)
      trial(k)=z(k)-2*h;emm=explicit_energy(n,trial)
      der(k)=(-epp+8*ep-8*em+emm)/(12*h*n)
     endif
    enddo
    der=der-g
   endif
   write(*,'(a,",",i0,",",*(es25.16e3,:,","))')trim(model_id),row,n,steps(j), &
    ind-energy_ref,pfd-p_ref,mufd-mu_ref,der,maxval(abs(pg)),state%gas_proton_margin,state%muon_margin
  enddo
 enddo
contains
 subroutine require_ok(context)
  character(*),intent(in)::context
  if(status/='ok')then
   write(0,*)trim(model_id),row,trim(context),trim(status)
   stop 4
  endif
 end subroutine
 function explicit_energy(nb,q)result(energy)
  real(wp),intent(in)::nb,q(4)
  real(wp)::energy,t,ni,yi,yg,ng,f,ne,ei,eg,ee,mn,mp,p,mue,sigma,phi,r,contrast
  t=q(1);ni=q(2);yi=q(3);yg=q(4)*yi
  f=t*nb/ni;ng=(nb-f*ni)/(1-f);ne=f*ni*yi+(1-f)*ng*yg
  ! Preserve the exact vacuum boundary despite algebraic roundoff.
  if(t==1._wp)ng=0
  call gdfm_native_baryons(model,ni,yi,ei,mn,mp,p,status);call require_ok('independent dense')
  call gdfm_native_baryons(model,ng,yg,eg,mn,mp,p,status);call require_ok('independent gas')
  call native_fermion(ne,me,ee,p,mue)
  sigma=sigma0*(b+16)/(b+yi**(-3)+(1-yi)**(-3))
  phi=(2-3*f**(1._wp/3)+f)/5
  contrast=ni*yi-ng*yg
  r=(3*sigma/(4*pi*e2*phi*contrast**2))**(1._wp/3)
  energy=f*ei+(1-f)*eg+ee+2*pi*e2*f*phi*r*r*contrast**2+3*f*sigma/r
  if(.not.ieee_is_finite(energy))stop 6
 end function
 function derivative(center,vals,h,edge)result(d)
  real(wp),intent(in)::center,vals(4),h
  integer,intent(in)::edge
  real(wp)::d
  if(edge/=0)then
   d=real(edge,wp)*(-25*center+48*vals(1)-36*vals(2)+16*vals(3)-3*vals(4))/(12*h)
  else
   d=(vals(1)-8*vals(2)+8*vals(3)-vals(4))/(12*h)
  endif
 end function
end program
