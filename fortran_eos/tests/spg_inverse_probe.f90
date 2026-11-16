! Diagnostic inverse recovery, never author-input attribution.
! Inputs: model, row, training flag, nb [fm^-3], Ye, A, Z.
! Conditional convention: A=4*pi*R^3*ni/3, y=Z/A, proton-free gas.
module spg_inverse_diagnostic
 use compose_eos_types,only:wp
 use compose_eos_gdfm_native
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_is_finite
 implicit none
 real(wp),parameter::pi=acos(-1._wp),e2=1.43996448_wp,me=.51099895_wp
 logical::excess=.false.
 type inverse_t
  real(wp)::nb=0,ye=0,aa=0,zz=0,y=0,t=0,ni=0,ng=0,f=0,r=0,sigma=0,ds=0,b=0,sigma0=0
  real(wp)::mech=0,chem=0,beta=0,d=0,dp=0,eps=0,pressure=0,mue=0,mnjump=0,cbeta=0
  integer::roots=0,status=0
 end type
contains
 subroutine evaluate(model,a,ni,h)
  type(gdfm_native_model_t),intent(in)::model
  type(inverse_t),intent(inout)::a
  real(wp),intent(in)::ni
  real(wp),intent(out)::h
  real(wp)::phi,phip,c,ei,eg,mni,mpi,mng,mpg,pii,pg,ee,pe,mue,x
  character(160)::ss
  a%ni=ni
  if(excess)then
   a%ng=a%nb*(1-a%ye*a%aa/a%zz)
   a%f=a%nb*a%ye*a%aa/a%zz/(ni-a%ng)
   a%r=(3*a%aa/(4*pi*(ni-a%ng)))**(1._wp/3)
   a%y=a%nb*a%ye/(a%f*ni);a%t=a%f*ni/a%nb
  else
   a%f=a%t*a%nb/ni;a%ng=a%nb*(1-a%t)/(1-a%f)
   a%r=(3*a%aa/(4*pi*ni))**(1._wp/3)
  endif
  x=a%f**(1._wp/3);phi=((1-a%f)/(1+x+x*x))**2*(x+2)/5
  phip=(1-a%f**(-2._wp/3))/5
  c=2*pi*e2*phi*a%r**2*(ni*a%y)**2 ! Coulomb density divided by f
  call gdfm_native_baryons(model,ni,a%y,ei,mni,mpi,pii,ss)
  call gdfm_native_baryons(model,a%ng,0._wp,eg,mng,mpg,pg,ss)
  call native_fermion(a%nb*a%ye,me,ee,pe,mue)
  h=pii-pg-c*(1+a%f*phip/phi)
  a%mue=mue;a%mnjump=mni-mng;a%cbeta=2*c/(ni*a%y)
  a%mech=h;a%chem=mni+a%y*(mpi-mni+mue)-mng+2*c/ni
  a%sigma=2*c*a%r/3
  a%ds=-a%r*ni/3*(mpi-mni+mue+2*c/(ni*a%y))
  a%d=a%y**(-3)+(1-a%y)**(-3)
  a%dp=-3*a%y**(-4)+3*(1-a%y)**(-4)
  a%b=-a%d-a%dp*a%sigma/a%ds
  a%sigma0=a%sigma*(a%b+a%d)/(a%b+16)
  a%eps=a%f*ei+(1-a%f)*eg+ee+3*a%f*c
  ! At a stationary interior cell muB=mng. The unfitted chemical residual
  ! diagnoses whether this pressure formula is an equilibrium interpretation.
  a%pressure=a%nb*mng-a%eps
 end subroutine
 subroutine invert(model,a)
  type(gdfm_native_model_t),intent(in)::model
  type(inverse_t),intent(inout)::a
  real(wp)::lo,hi,xl,xr,hl,hr,xm,hm,root
  integer::k,j
  a%status=1
  if(.not.all(ieee_is_finite([a%nb,a%ye,a%aa,a%zz])))return
  if(min(a%nb,a%ye,a%aa,a%zz)<=0.or.a%zz>=a%aa/2)return
  a%y=a%zz/a%aa;a%t=a%ye/a%y
  if(a%t>=1.or.a%t<=0)return
  lo=max(.03_wp,a%nb*(1+1e-6_wp));hi=.25_wp
  root=0;xl=lo;call evaluate(model,a,xl,hl)
  do k=1,160
   xr=lo+(hi-lo)*k/160
   call evaluate(model,a,xr,hr)
   if(hl*hr<0)then
    block
     real(wp)::bl,br,bhl
     bl=xl;br=xr;bhl=hl
     do j=1,55
      xm=(bl+br)/2;call evaluate(model,a,xm,hm)
      if(bhl*hm<=0)then
       br=xm
      else
       bl=xm;bhl=hm
      endif
     enddo
     a%roots=a%roots+1;root=(bl+br)/2
    end block
   endif
   xl=xr;hl=hr
  enddo
  if(a%roots==0)then
   a%status=2
   root=ieee_value(0._wp,ieee_quiet_nan)
   a%ni=root;a%ng=root;a%f=root;a%r=root;a%sigma=root;a%ds=root;a%b=root;a%sigma0=root
   a%mech=root;a%chem=root;a%d=root;a%dp=root;a%eps=root;a%pressure=root;a%mue=root;a%mnjump=root;a%cbeta=root
   return
  endif
  call evaluate(model,a,root,hm)
  a%status=0
 end subroutine
 subroutine vacuum(model,aa,zz,s0,b,ni,en,der)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::aa,zz,s0,b,ni
  real(wp),intent(out)::en,der
  real(wp)::eps,mn,mp,p,sigma,y,r,s,c
  character(160)::ss
  y=zz/aa;r=(3*aa/(4*pi*ni))**(1._wp/3)
  sigma=s0*(b+16)/(b+y**(-3)+(1-y)**(-3))
  call gdfm_native_baryons(model,ni,y,eps,mn,mp,p,ss)
  s=4*pi*r*r*sigma/aa;c=3*e2*zz*zz/(5*r*aa)
  en=eps/ni-gdfm_mass+s+c;der=p/ni**2+(-2*s+c)/(3*ni)
 end subroutine
 subroutine vacuum_min(model,aa,zz,s0,b,ni,en,der,edge)
  type(gdfm_native_model_t),intent(in)::model
  real(wp),intent(in)::aa,zz,s0,b
  real(wp),intent(out)::ni,en,der
  integer,intent(out)::edge
  real(wp)::xx(41),ee(41),lo,hi,x1,x2,e1,e2r,gr,tmp
  integer::i,k
  do i=1,41
   xx(i)=.04_wp+.2_wp*(i-1)/40
   call vacuum(model,aa,zz,s0,b,xx(i),ee(i),tmp)
  enddo
  k=minloc(ee,dim=1);edge=0
  if(k==1.or.k==41)then
   edge=1;ni=xx(k);en=ee(k);call vacuum(model,aa,zz,s0,b,ni,en,der);return
  endif
  lo=xx(k-1);hi=xx(k+1)
  call vacuum(model,aa,zz,s0,b,lo,e1,gr)
  call vacuum(model,aa,zz,s0,b,hi,e2r,tmp)
  if(gr>=0.or.tmp<=0)then
   edge=2;ni=xx(k);call vacuum(model,aa,zz,s0,b,ni,en,der);return
  endif
  do i=1,55
   ni=(lo+hi)/2;call vacuum(model,aa,zz,s0,b,ni,en,der)
   if(der>0)then
    hi=ni
   else
    lo=ni
   endif
  enddo
  ni=(lo+hi)/2;call vacuum(model,aa,zz,s0,b,ni,en,der)
 end subroutine

end module
program spg_inverse_probe
 use spg_inverse_diagnostic
 use compose_eos_spg_native
 implicit none
 integer,parameter::maxn=15000
 type(inverse_t)::a(maxn),syn
 type(gdfm_native_model_t)::model
 type(spg_native_result_t)::cell
 character(16)::ids(maxn),cmd,id
 character(160)::ss
 integer::rows(maxn),train(maxn),n,i,m,j,ios,cnt,mi,edge
 real(wp)::aa,zz,obs,ni,en,der
 real(wp)::sd,su,sdd,sdu,det,alpha,intercept,bfit,s0fit,sfit,dsfit,br,rr,den,vals(2)
 call get_command_argument(1,cmd)
 select case(trim(cmd))
 case('','full','excess','synthetic','mass','forward')
 case default
  write(0,'(a)')'unknown diagnostic mode';stop 3
 end select
 excess=trim(cmd)=="excess"
 if(trim(cmd)=='forward')then
  print '(a)','model,row,nb,epsilon,pressure,Ye,A,Z,ni,f,projected_gradient,min_curvature,proton_insertion_margin,attempts,stationary,status'
  do
   read(*,*,iostat=ios)id,j,ni,s0fit,bfit
   if(ios>0)stop 4
   if(ios<0)exit
   call load_gdfm_native(trim(id),model,ss)
   call spg_native_solve(model,s0fit,bfit,ni,cell,ss)
   write(*,'(a,a,i0,11(a,es24.15e3),a,i0,a,l1,a,a)')trim(id),',',j,',',ni,',',cell%epsilon,',',cell%pressure, &
    ',',cell%ye,',',cell%aav,',',cell%zav,',',cell%z(2),',',cell%volume_fraction, &
    ',',cell%max_projected_gradient,',',cell%min_scaled_curvature,',',cell%gas_proton_margin, &
    ',',cell%attempt_count,',',cell%stationary,',',trim(ss)
  enddo
  stop
 endif
 if(trim(cmd)=='mass')then
  print '(a)','model,Z,A,observed_E_over_A,predicted_E_over_A,n_opt,dE_dn,binding_residual_MeV,boundary'
  do
   read(*,*,iostat=ios)id,s0fit,bfit,zz,aa,obs
   if(ios>0)stop 4
   if(ios<0)exit
   call load_gdfm_native(trim(id),model,ss)
   call vacuum_min(model,aa,zz,s0fit,bfit,ni,en,der,edge)
   write(*,'(a,7(a,es24.15e3),a,i0)')trim(id),',',zz,',',aa,',',obs,',',en,',',ni,',',der,',',(obs-en)*aa,',',edge
  enddo
  stop
 endif
 if(trim(cmd)=='synthetic')then
  print '(a)','model index status nb sigma0_local b_local ni_error chemical_residual mechanical_residual forward_gradient full_beta_residual factor1minusf_beta_residual'
  do m=1,5
   write(id,'(a,i1)')'SPG_M',m
   call load_gdfm_native(trim(id),model,ss)
   vals=[.008_wp,.03_wp]
   do j=1,2
    call spg_native_solve(model,1.1_wp,25._wp,vals(j),cell,ss)
    if(ss/='ok')then
     print *,trim(id),vals(j),'FORWARD_FAIL ',trim(ss);cycle
    endif
    syn=inverse_t();syn%nb=cell%nb;syn%ye=cell%ye;syn%aa=cell%aav;syn%zz=cell%zav
    call invert(model,syn)
    sfit=1.1_wp*41._wp/(25._wp+syn%d);dsfit=-sfit*syn%dp/(25._wp+syn%d)
    br=3*(dsfit-syn%ds)/(syn%r*syn%ni);rr=br-3*syn%f*dsfit/(syn%r*syn%ni)
    print '(a,1x,2i4,1x,*(es24.15e3,1x))',trim(id),j,syn%status,vals(j),syn%sigma0,syn%b, &
      syn%ni-cell%z(2),syn%chem,syn%mech,cell%max_projected_gradient,br,rr
   enddo
  enddo
  stop
 endif
 n=0
 do
  i=n+1
  if(i>maxn)stop 5
  read(*,*,iostat=ios)ids(i),rows(i),train(i),a(i)%nb,a(i)%ye,a(i)%aa,a(i)%zz
  if(ios>0)stop 4
  if(ios<0)exit
  n=i
  call load_gdfm_native(trim(ids(i)),model,ss)
  if(ss/='ok')stop 2
  call invert(model,a(i))
 enddo
 print '(a)','model,row,train,status,roots,nb,Ye,A,Z,y,t,ni,ng,f,R,sigma,dsigma,b_local,sigma0_local,mech_MeV_fm3,chemical_MeV,sigma_fit,dsigma_fit,beta_fit_MeV,sigma_relres,epsilon_MeV_fm3,pressure_MeV_fm3,mue_MeV,bulk_neutron_jump_MeV,coulomb_beta_MeV,beta_factor1minusf_MeV,neutron_equilibrium_MeV'
 do m=1,5
  write(id,'(a,i1)')'SPG_M',m
  cnt=0;sd=0;su=0;sdd=0;sdu=0
  do i=1,n
   if(ids(i)/=id.or.train(i)/=1.or.a(i)%status/=0)cycle
   cnt=cnt+1;sd=sd+a(i)%d;su=su+1/a(i)%sigma;sdd=sdd+a(i)%d**2;sdu=sdu+a(i)%d/a(i)%sigma
  enddo
  if(cnt<2)cycle
  det=cnt*sdd-sd*sd;alpha=(cnt*sdu-sd*su)/det;intercept=(sdd*su-sd*sdu)/det
  bfit=intercept/alpha;s0fit=1/(alpha*(bfit+16))
  write(0,'(a,1x,i4,1x,*(es24.15e3,1x))')trim(id),cnt,s0fit,bfit,det/(cnt*sdd)
  do i=1,n
   if(ids(i)/=id)cycle
   sfit=s0fit*(bfit+16)/(bfit+a(i)%d);dsfit=-sfit*a(i)%dp/(bfit+a(i)%d)
   br=3*(dsfit-a(i)%ds)/(a(i)%r*a(i)%ni);rr=(sfit-a(i)%sigma)/a(i)%sigma
   write(*,'(a,4(a,i0),27(a,es24.15e3))')trim(id),',',rows(i),',',train(i),',',a(i)%status,',',a(i)%roots, &
    ',',a(i)%nb,',',a(i)%ye,',',a(i)%aa,',',a(i)%zz,',',a(i)%y,',',a(i)%t,',',a(i)%ni,',',a(i)%ng,',',a(i)%f, &
    ',',a(i)%r,',',a(i)%sigma,',',a(i)%ds,',',a(i)%b,',',a(i)%sigma0,',',a(i)%mech,',',a(i)%chem, &
    ',',sfit,',',dsfit,',',br,',',rr,',',a(i)%eps,',',a(i)%pressure,',',a(i)%mue,',',a(i)%mnjump,',',a(i)%cbeta, &
    ',',br-3*a(i)%f*dsfit/(a(i)%r*a(i)%ni),',',a(i)%mnjump-3*a(i)%y*dsfit/(a(i)%r*a(i)%ni)
  enddo
 enddo
end program
