program spg_finish_drip_curvature
 use compose_eos_types,only:wp
 use compose_eos_gdfm_native
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 real(wp),parameter::pi=acos(-1._wp),e2=1.43996448_wp,me=.51099895_wp
 type(gdfm_native_model_t)::model
 character(32)::model_id
 character(160)::status
 integer::row,critical,edge,ios,j,k,l
 real(wp)::n,sigma0,b,z(4),er,pr,mr,base,h,q(4),a,c,d,mineig,ep,em
 real(wp),parameter::steps(3)=[1e-4_wp,3e-5_wp,1e-5_wp]
 write(*,'(a)')'model,row,h_density_and_fraction,Hni_ni_scaled,Hni_y_scaled,Hy_y_scaled,min_scaled_curvature'
 do
  read(*,*,iostat=ios)model_id,row,critical,edge,n,sigma0,b,z,er,pr,mr
  if(ios<0)exit
  if(ios/=0)stop 2
  if(row/=1)cycle
  call load_gdfm_native(trim(model_id),model,status);call require_ok('load')
  base=explicit_energy(n,z)
  do j=1,3
   h=steps(j);q=z;q(2)=z(2)+h;ep=explicit_energy(n,q)
   q(2)=z(2)-h;em=explicit_energy(n,q);a=(ep-2*base+em)/(n*h*h)*.01_wp
   q=z;q(3)=z(3)+h;ep=explicit_energy(n,q)
   q(3)=z(3)-h;em=explicit_energy(n,q);d=(ep-2*base+em)/(n*h*h)*.01_wp
   q=z;q(2:3)=z(2:3)+h;c=explicit_energy(n,q)
   q=z;q(2)=z(2)+h;q(3)=z(3)-h;c=c-explicit_energy(n,q)
   q=z;q(2)=z(2)-h;q(3)=z(3)+h;c=c-explicit_energy(n,q)
   q=z;q(2:3)=z(2:3)-h;c=(c+explicit_energy(n,q))/(4*n*h*h)*.01_wp
   mineig=(a+d-sqrt((a-d)**2+4*c*c))/2
   write(*,'(a,",",i0,",",*(es25.16e3,:,","))')trim(model_id),row,h,a,c,d,mineig
   if(mineig<=0)stop 3
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
end program
