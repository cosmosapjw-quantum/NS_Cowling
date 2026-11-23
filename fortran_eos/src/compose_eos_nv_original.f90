! Original NV source binding and reusable Table-1 bulk functional.
! Negele & Vautherin, Nucl.Phys.A207 (1973) 298-320, DOI 10.1016/0375-9474(73)90349-7.
! Ground-state provider: original Table4/p317 author fit, numerically identical
! to NV1974_author_fit. The Table1 uniform reconstruction below is a diagnostic
! component, not a self-consistent inhomogeneous HF crust provider.
module compose_eos_nv_original
 use compose_eos_types,only:wp,eos_data_t,allocate_eos
 use compose_eos_nv_primary,only:nv_primary_fit_t,nv_primary_state_t,nv_primary_state
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 character(*),parameter::pdf_sha='5611ae66700d983e7f66881839802c2204e16f91cafbf19c670fbd178042bb5b'
 real(wp),parameter::mn=939.56542052_wp,mp=938.27208816_wp,me=.51099895_wp
 real(wp),parameter::pi=acos(-1._wp),hc=197.3269804_wp,alpha=1._wp/137.035999084_wp
 real(wp),parameter::clo=2.79e-4_wp,chi=7.89e-2_wp,cref=1e-4_wp
 real(wp),parameter::cg(0:7)=[-4._wp,.28822899_wp,.59150523_wp,.090185940_wp, &
  -.11025614_wp,.029377479_wp,-.0032618465_wp,.00013543555_wp]
 real(wp),parameter::cu(0:7)=[-4._wp,1.4821424_wp,-.040373482_wp,.060455728_wp, &
  -.015307639_wp,.0034774416_wp,-.00043627154_wp,.000023383473_wp]
 real(wp),parameter::table3(9,11)=reshape([ &
  2.79e35_wp,140._wp,40._wp,.2_wp,-26.8_wp,4e-5_wp,.53_wp,-1.425_wp,.436_wp, &
  4e35_wp,160._wp,40._wp,.3_wp,-29.4_wp,9.7e-5_wp,.53_wp,-.962_wp,.543_wp, &
  6e35_wp,210._wp,40._wp,.6_wp,-29.5_wp,2.6e-4_wp,.53_wp,-.462_wp,.692_wp, &
  8.79e35_wp,280._wp,40._wp,1._wp,-28.5_wp,4.8e-4_wp,.53_wp,-.050_wp,.865_wp, &
  1.59e36_wp,460._wp,40._wp,1.4_wp,-29.4_wp,.0012_wp,.52_wp,.541_wp,1.214_wp, &
  3.73e36_wp,900._wp,50._wp,2.6_wp,-33.6_wp,.003_wp,.46_wp,1.465_wp,1.926_wp, &
  5.77e36_wp,1050._wp,50._wp,3.3_wp,-34.5_wp,.0047_wp,.45_wp,1.996_wp,2.408_wp, &
  8.91e36_wp,1300._wp,50._wp,4.2_wp,-35.8_wp,.0078_wp,.44_wp,2.610_wp,2.981_wp, &
  2.04e37_wp,1750._wp,50._wp,6.5_wp,-43.6_wp,.0184_wp,.35_wp,4.097_wp,4.422_wp, &
  4.75e37_wp,1460._wp,40._wp,10.9_wp,-54._wp,.0436_wp,.28_wp,6.428_wp,6.660_wp, &
  7.89e37_wp,950._wp,32._wp,15._wp,-68.3_wp,.0737_wp,.16_wp,8.611_wp,8.657_wp],[9,11])
 real(wp),parameter::vt(28)=[0._wp,0._wp,-1.29058e3_wp,-2.60473e2_wp,3.48227e3_wp,5.41286e4_wp, &
  2.70730e3_wp,2.08122e4_wp,-9.01970e5_wp,1.68268e7_wp, &
  -2.97146e4_wp,-9.24196e4_wp,-1.69465e6_wp,-2.96778e7_wp,3.35473e8_wp, &
  1.83110e5_wp,8.17815e5_wp,2.40883e7_wp,-2.72504e8_wp,-2.39553e9_wp,4.55011e10_wp, &
  -4.30956e5_wp,-4.77095e6_wp,-2.89827e7_wp,2.15640e9_wp,-1.62585e10_wp,1.29625e9_wp,-1.88806e11_wp]
 type,public::nv_original_source_t
  type(nv_primary_fit_t)::fit
  real(wp)::uniform_coeff(0:7)=0,vnm(0:6,0:6)=0
  logical::valid=.false.
 end type
 type,public::nv_uniform_state_t
  real(wp)::nb=0,nn=0,np=0,e_offset=0,energy=0,pressure=0,mu_n=0,mu_p=0,mu_e=0
  real(wp)::beta_residual=0,beta_curvature=0,bulk_v=0,exchange=0
  logical::beta_curvature_defined=.false.
 end type
 public::load_nv_original_source,generate_nv_original_native,nv1973_bulk_potential
 public::nv1973_uniform_composition,nv1973_uniform_beta_state
contains
 subroutine tokens_read(u,t,nt,eof,status)
  integer,intent(in)::u
  character(*),intent(out)::t(:),status
  integer,intent(out)::nt
  logical,intent(out)::eof
  character(4096)::line
  integer::ios,i,j,ll
  status='ok';nt=0;t='';eof=.false.
  do
   read(u,'(a)',iostat=ios)line
   if(ios<0)then;eof=.true.;return;endif
   if(ios/=0)then;status='NV1973 source read failed';return;endif
   if(len_trim(line)==len(line))then;status='NV1973 source line too long';return;endif
   line=adjustl(line);ll=len_trim(line)
   if(ll==0)cycle
   if(line(1:1)=='#'.or.line(1:1)=='!')cycle
   exit
  enddo
  i=1
  do while(i<=ll)
   if(index(' '//achar(9)//achar(13),line(i:i))/=0)then;i=i+1;cycle;endif
   j=i
   do while(j<=ll)
    if(index(' '//achar(9)//achar(13),line(j:j))/=0)exit
    j=j+1
   enddo
   nt=nt+1
   if(nt>size(t))then;status='NV1973 source excess columns';return;endif
   if(j-i>len(t(1)))then;status='NV1973 source token too long';return;endif
   t(nt)=line(i:j-1);i=j
  enddo
 end subroutine
 subroutine number(text,x,status)
  character(*),intent(in)::text
  real(wp),intent(out)::x
  character(*),intent(out)::status
  integer::i,ie,ne,ex,ios
  status='NV1973 invalid numeric token';x=0
  if(len_trim(text)==0)return
  do i=1,len_trim(text)
   if(index('0123456789+-.eEdD',text(i:i))==0)return
  enddo
  ie=scan(text,'eEdD')
  if(ie>0)then
   ne=len_trim(text)-ie
   if(ne<1.or.ne>6)return
   read(text(ie+1:),*,iostat=ios)ex
   if(ios/=0)return
   if(abs(ex)>100)return
  endif
  read(text,*,iostat=ios)x
  if(ios/=0)return
  if(.not.ieee_is_finite(x))return
  status='ok'
 end subroutine
 subroutine load_nv_original_source(path,src,status)
  character(*),intent(in)::path
  type(nv_original_source_t),intent(out)::src
  character(*),intent(out)::status
  character(128)::t(16)
  real(wp)::v(9)
  logical::eof
  integer::u,ios,nt,r,k,i,j
  src=nv_original_source_t();status='NV1973 source cannot be opened'
  open(newunit=u,file=path,status='old',action='read',iostat=ios)
  if(ios/=0)return
  i=0;j=0
  do r=1,49
   call tokens_read(u,t,nt,eof,status)
   if(status/='ok')exit
   if(eof)then;status='NV1973 source truncated';exit;endif
   select case(r)
   case(1)
    if(nt/=1.or.t(1)/='NV1973_ORIGINAL_SOURCE_V1')then;status='NV1973 source schema mismatch';exit;endif
   case(2)
    if(nt/=2.or.t(1)/='DOCUMENT_SHA256'.or.t(2)/=pdf_sha)then;status='NV1973 document identity mismatch';exit;endif
   case(3,4)
    if(nt/=3)then;status='NV1973 masses/domain require two values';exit;endif
    if(r==3.and.t(1)/='MASS_MEV')then;status='NV1973 mass row mismatch';exit;endif
    if(r==4.and.t(1)/='DENSITY_DOMAIN_FM3')then;status='NV1973 domain row mismatch';exit;endif
    do k=1,2
     call number(t(k+1),v(k),status)
     if(status/='ok')exit
    enddo
    if(status/='ok')exit
    if(r==3)then
     src%fit%mn=v(1);src%fit%mp=v(2)
    else
     src%fit%nlo=v(1);src%fit%nhi=v(2)
    endif
   case(5)
    if(nt/=2.or.t(1)/='LOG_REFERENCE_FM3')then;status='NV1973 log reference row mismatch';exit;endif
    call number(t(2),src%fit%nref,status)
    if(status/='ok')exit
   case(6,7)
    if(nt/=9)then;status='NV1973 fit coefficients require eight values';exit;endif
    if(r==6.and.t(1)/='GROUND_COEFFICIENTS')then;status='NV1973 ground coefficient marker mismatch';exit;endif
    if(r==7.and.t(1)/='UNIFORM_COEFFICIENTS')then;status='NV1973 uniform coefficient marker mismatch';exit;endif
    do k=1,8
     call number(t(k+1),v(k),status)
     if(status/='ok')exit
    enddo
    if(status/='ok')exit
    if(r==6)then;src%fit%coeff=v(1:8);else;src%uniform_coeff=v(1:8);endif
   case(8)
    if(nt/=1.or.t(1)/='TABLE3_11X9')then;status='NV1973 Table3 marker mismatch';exit;endif
   case(9:19)
    if(nt/=9)then;status='NV1973 Table3 requires nine columns';exit;endif
    do k=1,9
     call number(t(k),v(k),status)
     if(status/='ok')exit
    enddo
    if(status/='ok')exit
    src%fit%table(:,r-8)=v
   case(20)
    if(nt/=1.or.t(1)/='TABLE1_LOWER_28')then;status='NV1973 Table1 marker mismatch';exit;endif
   case(21:48)
    if(nt/=3)then;status='NV1973 Table1 requires I,J,value';exit;endif
    do k=1,3
     call number(t(k),v(k),status)
     if(status/='ok')exit
    enddo
    if(status/='ok')exit
    if(v(1)/=real(i,wp).or.v(2)/=real(j,wp))then;status='NV1973 Table1 index order mismatch';exit;endif
    if(v(3)/=vt(r-20))then;status='NV1973 Table1 coefficient differs from source';exit;endif
    src%vnm(i,j)=v(3);src%vnm(j,i)=v(3)
    j=j+1
    if(j>i)then;i=i+1;j=0;endif
   case(49)
    if(nt/=1.or.t(1)/='END_NV1973_SOURCE')then;status='NV1973 source end marker mismatch';exit;endif
   end select
  enddo
  if(status=='ok')then
   call tokens_read(u,t,nt,eof,status)
   if(status=='ok'.and..not.eof)status='NV1973 source trailing data'
  endif
  close(u,iostat=ios)
  if(status/='ok')return
  if(ios/=0)then;status='NV1973 source close failed';return;endif
  if(src%fit%mn/=mn.or.src%fit%mp/=mp)then;status='NV1973 declared project mass convention mismatch';return;endif
  if(src%fit%nlo/=clo.or.src%fit%nhi/=chi.or.src%fit%nref/=cref)then;status='NV1973 domain/log convention mismatch';return;endif
  if(any(src%fit%coeff/=cg).or.any(src%uniform_coeff/=cu))then;status='NV1973 Table4 coefficient identity mismatch';return;endif
  if(any(src%fit%table/=table3))then;status='NV1973 Table3 identity mismatch; row7 energy is 1.996';return;endif
  src%fit%valid=.true.;src%valid=.true.;status='ok'
 end subroutine

 subroutine nv1973_bulk_potential(src,nn,np,v,gradient,hessian,status)
  type(nv_original_source_t),intent(in)::src
  real(wp),intent(in)::nn,np
  real(wp),intent(out)::v,gradient(2),hessian(2,2)
  character(*),intent(out)::status
  real(wp)::a,b,c,oldv,oldp
  integer::i,j
  v=0;gradient=0;hessian=0;status='NV1973 bulk requires validated source'
  if(.not.src%valid)return
  status='NV1973 bulk polynomial diagnostic domain requires 0<=nn,np<=0.15 fm^-3'
  if(.not.all(ieee_is_finite([nn,np])))return
  if(nn<0.or.np<0.or.nn>.15_wp.or.np>.15_wp)return
  ! Nested differentiated Horner. Evaluate only the published Table1 entries
  ! I,J=0..6; missing index-7 terms in the older prose notation are not invented.
  do i=6,0,-1
   a=src%vnm(i,6);b=0;c=0
   do j=5,0,-1
    c=c*np+2*b;b=b*np+a;a=a*np+src%vnm(i,j)
   enddo
   oldv=v;oldp=gradient(2)
   hessian(1,1)=hessian(1,1)*nn+2*gradient(1)
   gradient(1)=gradient(1)*nn+oldv
   hessian(1,2)=hessian(1,2)*nn+oldp
   gradient(2)=gradient(2)*nn+b
   hessian(2,2)=hessian(2,2)*nn+c
   v=v*nn+a
  enddo
  hessian(2,1)=hessian(1,2)
  if(.not.all(ieee_is_finite([v,gradient,reshape(hessian,[4])])) )then
   status='NV1973 bulk polynomial nonfinite';return
  endif
  status='ok'
 end subroutine

 pure subroutine electron_gas(ne,energy,pressure,mu)
  real(wp),intent(in)::ne
  real(wp),intent(out)::energy,pressure,mu
  real(wp)::x,x2,root,fac
  if(ne==0)then;energy=0;pressure=0;mu=me;return;endif
  x=hc*(3*pi*pi*ne)**(1._wp/3)/me;x2=x*x;root=sqrt(1+x2);mu=me*root
  if(x<.01_wp)then
   energy=ne*me*(1+3*x2/10-3*x2*x2/56+x2**3/48-15*x2**4/1408)
   pressure=ne*me*(x2/5-x2*x2/14+x2**3/24-5*x2**4/176)
  else
   fac=me**4/(8*pi*pi*hc**3)
   energy=fac*(x*(1+2*x2)*root-asinh(x))
   pressure=fac/3*(x*(2*x2-3)*root+3*asinh(x))
  endif
 end subroutine

 subroutine nv1973_uniform_composition(src,nb,np,state,status)
  type(nv_original_source_t),intent(in)::src
  real(wp),intent(in)::nb,np
  type(nv_uniform_state_t),intent(out)::state
  character(*),intent(out)::status
  real(wp)::nn,v,dv(2),hh(2,2),tcoef,kn,kp,tn,tp,ee,pe,mue,ec,upc,ecpp,xcoef,np13
  state=nv_uniform_state_t();status='NV1973 uniform diagnostic domain/composition invalid'
  if(.not.all(ieee_is_finite([nb,np])))return
  if(nb<clo.or.nb>chi.or.np<0.or.np>=nb)return
  nn=nb-np
  call nv1973_bulk_potential(src,nn,np,v,dv,hh,status)
  if(status/='ok')return
  ! Uniform limit of p304: A+B_n*tau_n+B_p*tau_p=V^NM.
  ! Common nucleon kinetic mass M=m_n is explicit; species rest masses differ.
  tcoef=hc*hc/(2*mn)*(3*pi*pi)**(2._wp/3)
  kn=tcoef*nn**(2._wp/3);kp=tcoef*np**(2._wp/3)
  tn=3._wp/5*nn*kn;tp=3._wp/5*np*kp
  xcoef=.75_wp*(3/pi)**(1._wp/3)*alpha*hc
  np13=np**(1._wp/3);ec=-xcoef*np*np13
  upc=0;ecpp=0
  if(np>0)then
   ! Differentiate powers before division: np**2 underflows for tiny inputs.
   upc=-4._wp/3*xcoef*np13;ecpp=-4._wp/9*xcoef/(np13*np13)
  endif
  call electron_gas(np,ee,pe,mue)
  state%nb=nb;state%nn=nn;state%np=np;state%bulk_v=v;state%exchange=ec
  state%energy=mn*nn+mp*np+tn+tp+v+ec+ee
  ! Store the offset without subtracting two ~939 MeV numbers.
  state%e_offset=((mp-mn)*np+tn+tp+v+ec+ee)/nb
  state%pressure=2._wp/3*(tn+tp)+nn*dv(1)+np*dv(2)-v+ec/3+pe
  state%mu_n=mn+kn+dv(1);state%mu_p=mp+kp+dv(2)+upc;state%mu_e=mue
  state%beta_residual=state%mu_p+state%mu_e-state%mu_n
  if(np>0)then
   state%beta_curvature_defined=.true.
   ! Use the electron Fermi momentum, avoiding cancellation in mu_e**2-me**2.
   state%beta_curvature=hh(1,1)+hh(2,2)-2*hh(1,2)+ &
    2._wp/3*(kn/nn+kp/np)+hc*hc*(3*pi*pi)**(2._wp/3)/(3*mue*np13)+ecpp
  endif
  if(.not.all(ieee_is_finite([state%energy,state%e_offset,state%pressure, &
   state%mu_n,state%mu_p,state%mu_e,state%beta_residual,state%beta_curvature])))then
   status='NV1973 uniform diagnostic is nonfinite';return
  endif
  status='ok'
 end subroutine

 subroutine nv1973_uniform_beta_state(src,nb,state,status)
  type(nv_original_source_t),intent(in)::src
  real(wp),intent(in)::nb
  type(nv_uniform_state_t),intent(out)::state
  character(*),intent(out)::status
  type(nv_uniform_state_t)::lo_state,hi_state
  real(wp)::lo,hi,mid
  integer::k
  state=nv_uniform_state_t()
  call nv1973_uniform_composition(src,nb,0._wp,lo_state,status)
  if(status/='ok')return
  call nv1973_uniform_composition(src,nb,nb/2,hi_state,status)
  if(status/='ok')return
  if(lo_state%beta_residual>=0.or.hi_state%beta_residual<=0)then
   status='NV1973 uniform beta root not bracketed on source neutron-rich branch';return
  endif
  lo=0;hi=nb/2
  do k=1,90
   mid=(lo+hi)/2
   call nv1973_uniform_composition(src,nb,mid,state,status)
   if(status/='ok')return
   if(abs(state%beta_residual)<2e-12_wp)exit
   if(mid==lo.or.mid==hi)exit
   if(state%beta_residual<0)then;lo=mid;else;hi=mid;endif
  enddo
  if(abs(state%beta_residual)>2e-10_wp.or.state%beta_curvature<=0)then
   status='NV1973 uniform beta stationary/minimum residual failed';return
  endif
  status='ok'
 end subroutine

 subroutine finish_write(u,ios,status,label)
  integer,intent(in)::u
  integer,intent(inout)::ios
  character(*),intent(out)::status
  character(*),intent(in)::label
  integer::cs
  close(u,iostat=cs)
  if(ios==0)ios=cs
  status='ok'
  if(ios/=0)status='NV1973 cannot write '//label
 end subroutine

 subroutine generate_nv_original_native(source_path,n_points,eos,output_dir,status)
  character(*),intent(in)::source_path,output_dir
  integer,intent(in)::n_points
  type(eos_data_t),intent(out)::eos
  character(*),intent(out)::status
  type(nv_original_source_t)::src
  type(nv_primary_fit_t)::uniform_fit
  type(nv_primary_state_t)::s,sg
  type(nv_primary_state_t),allocatable::states(:)
  type(nv_uniform_state_t)::uni
  character(64),parameter::files(4)=[character(64)::'nv1973_states.csv','nv1973_table3_comparison.csv', &
   'nv1973_uniform_reconstruction.csv','NV_ORIGINAL_GENERATION.json']
  real(wp)::n,p_euler,e74,max_e,max_gas_e,max_uniform,max_beta,min_mu,max_euler
  real(wp)::v,dv(2),hh(2,2),max_poly_symmetry
  integer::u,ios,i
  logical::exists
  eos=eos_data_t();status='NV1973 grid requires 3..100000 points'
  if(n_points<3.or.n_points>100000)return
  call load_nv_original_source(source_path,src,status)
  if(status/='ok')return
  do i=1,size(files)
   inquire(file=trim(output_dir)//'/'//trim(files(i)),exist=exists)
   if(exists)then;status='NV1973 diagnostic file already exists';return;endif
  enddo
  allocate(states(n_points));call allocate_eos(eos,n_points)
  eos%mn=mn;eos%mp=mp;eos%il=1;eos%beta_equilibrated=.true.
  eos%composition_available=.false.;eos%composition_file_present=.false.
  do i=1,n_points
   n=exp(log(clo)+real(i-1,wp)/real(n_points-1,wp)*log(chi/clo))
   if(i==1)n=clo
   if(i==n_points)n=chi
   call nv_primary_state(src%fit,n,states(i),status)
   if(status/='ok')return
   eos%nb(i)=n;eos%q(1,i)=states(i)%pressure/n
   eos%q(3,i)=(states(i)%e_offset+eos%q(1,i))/mn
   eos%q(6,i)=states(i)%e_offset/mn;eos%q(7,i)=eos%q(6,i)
  enddo
  open(newunit=u,file=trim(output_dir)//'/'//trim(files(1)),status='new',action='write',iostat=ios)
  if(ios/=0)then;status='NV1973 cannot create state diagnostics';return;endif
  write(u,'(a)',iostat=ios)'n_b_fm3,e_offset_mev,energy_mev_fm3,pressure_mev_fm3,mu_b_mev,cs2_over_c2,gamma'
  do i=1,n_points
   if(ios/=0)exit
   s=states(i)
   write(u,'(*(es25.16e3,:,","))',iostat=ios)s%nb,s%e_offset,s%epsilon,s%pressure,s%mu_b,s%sound2,s%gamma
  enddo
  call finish_write(u,ios,status,'state diagnostics')
  if(status/='ok')return
  uniform_fit=src%fit;uniform_fit%coeff=src%uniform_coeff
  max_e=0;max_gas_e=0;max_uniform=0;max_beta=0;max_poly_symmetry=0
  open(newunit=u,file=trim(output_dir)//'/'//trim(files(2)),status='new',action='write',iostat=ios)
  if(ios/=0)then;status='NV1973 cannot create Table3 comparison';return;endif
  write(u,'(a)',iostat=ios)'n_b_fm3,e_1973_table3_mev,e_1974_tableII_mev,e_ground_fit_mev,'// &
   'ground_fit_minus_1973_mev,mu_n_HF_table_mev,mu_fit_offset_mev,pressure_fit_mev_fm3,'// &
   'pressure_table_euler_diagnostic_mev_fm3,e_uniform_table3_mev,e_uniform_table4_fit_mev'
  do i=1,11
   if(ios/=0)exit
   n=src%fit%table(1,i)*1e-39_wp
   if(i==1)n=clo
   if(i==11)n=chi
   call nv_primary_state(src%fit,n,s,status)
   if(status/='ok')exit
   call nv_primary_state(uniform_fit,n,sg,status)
   if(status/='ok')exit
   p_euler=n*(src%fit%table(4,i)-src%fit%table(8,i))
   e74=src%fit%table(8,i)
   if(i==7)e74=1.966_wp
   write(u,'(*(es25.16e3,:,","))',iostat=ios)n,src%fit%table(8,i),e74,s%e_offset, &
    s%e_offset-src%fit%table(8,i),src%fit%table(4,i),s%e_offset+s%pressure/n,s%pressure, &
    p_euler,src%fit%table(9,i),sg%e_offset
   max_e=max(max_e,abs(s%e_offset-src%fit%table(8,i)))
   max_gas_e=max(max_gas_e,abs(sg%e_offset-src%fit%table(9,i)))
  enddo
  if(status/='ok')then;close(u,iostat=ios);return;endif
  call finish_write(u,ios,status,'Table3 comparison')
  if(status/='ok')return
  open(newunit=u,file=trim(output_dir)//'/'//trim(files(3)),status='new',action='write',iostat=ios)
  if(ios/=0)then;status='NV1973 cannot create uniform reconstruction';return;endif
  write(u,'(a)',iostat=ios)'n_b_fm3,n_n_fm3,n_p_fm3,e_uniform_reconstructed_mev,e_uniform_table3_mev,'// &
   'reconstructed_minus_table3_mev,energy_mev_fm3,pressure_mev_fm3,mu_n_mev,mu_p_mev,mu_e_mev,'// &
   'beta_residual_mev,beta_curvature_mev_fm3,V_NM_mev_fm3,proton_slater_exchange_mev_fm3'
  do i=1,11
   if(ios/=0)exit
   n=src%fit%table(1,i)*1e-39_wp
   if(i==1)n=clo
   if(i==11)n=chi
   call nv1973_uniform_beta_state(src,n,uni,status)
   if(status/='ok')exit
   call nv1973_bulk_potential(src,uni%np,uni%nn,v,dv,hh,status)
   if(status/='ok')exit
   max_poly_symmetry=max(max_poly_symmetry,abs(v-uni%bulk_v))
   write(u,'(*(es25.16e3,:,","))',iostat=ios)n,uni%nn,uni%np,uni%e_offset,src%fit%table(9,i), &
    uni%e_offset-src%fit%table(9,i),uni%energy,uni%pressure,uni%mu_n,uni%mu_p,uni%mu_e, &
    uni%beta_residual,uni%beta_curvature,uni%bulk_v,uni%exchange
   max_uniform=max(max_uniform,abs(uni%e_offset-src%fit%table(9,i)))
   max_beta=max(max_beta,abs(uni%beta_residual))
  enddo
  if(status/='ok')then;close(u,iostat=ios);return;endif
  call finish_write(u,ios,status,'uniform reconstruction')
  if(status/='ok')return
  min_mu=minval(states%mu_b-mn)
  max_euler=maxval(abs((states%epsilon+states%pressure)/states%nb-states%mu_b))
  open(newunit=u,file=trim(output_dir)//'/'//trim(files(4)),status='new',action='write',iostat=ios)
  if(ios/=0)then;status='NV1973 cannot create receipt';return;endif
  write(u,'(a)',iostat=ios)'{"schema":"ns_cowling.nv_original_generation.v1",'
  if(ios==0)write(u,'(a)',iostat=ios)'"model":"NV1973_author_fit","doi":"10.1016/0375-9474(73)90349-7",'
  if(ios==0)write(u,'(a)',iostat=ios)'"status":"ORIGINAL_PUBLISHED_FIT_GENERATED_WITH_NATIVE_BULK_RECONSTRUCTION",'
  if(ios==0)write(u,'(a)',iostat=ios)'"original_pdf_sha256":"'//pdf_sha//'",'
  if(ios==0)write(u,'(a)',iostat=ios)'"source_pages":[304,305,306,307,308,309,310,311,317],'
  if(ios==0)write(u,'(a)',iostat=ios)'"published_fit_semantically_equal_to":"NV1974_author_fit","semantic_equivalent_model":"NV1974_author_fit",'
  if(ios==0)write(u,'(a)',iostat=ios)'"distinct_new_physical_eos":false,'
  if(ios==0)write(u,'(a)',iostat=ios)'"kernel":"compose_eos_nv_primary:nv_primary_state","source_coefficients_modified":false,'
  if(ios==0)write(u,'(a)',iostat=ios)'"original_table3_row7_energy_mev":1.996,"1974_tableII_row7_energy_mev":1.966,'
  if(ios==0)write(u,'(a)',iostat=ios)'"composition_available":false,"electron_energy_included":true,"extrapolation":false,'
  if(ios==0)write(u,'(a)',iostat=ios)'"mass_convention":"existing_project_mn939.56542052_mp938.27208816_MeV",'
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)'"rows":',n_points,','
  if(ios==0)write(u,'(a,es25.16e3,a,es25.16e3,a)',iostat=ios)'"density_domain_fm3":[',clo,',',chi,'],'
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_ground_fit_minus_table3_mev":',max_e,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_uniform_fit_minus_table3_mev":',max_gas_e,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_uniform_reconstruction_minus_table3_mev":',max_uniform,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_uniform_beta_residual_mev":',max_beta,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_bulk_exchange_symmetry_roundoff_mev_fm3":',max_poly_symmetry,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"max_ground_euler_roundoff_mev":',max_euler,','
  if(ios==0)write(u,'(a,es25.16e3,a)',iostat=ios)'"minimum_fit_mu_minus_neutron_rest_mev":',min_mu,','
  if(ios==0)write(u,'(a,i0,a)',iostat=ios)'"sampled_fit_mu_below_neutron_rest_rows":',count(states%mu_b<mn),','
  if(ios==0)write(u,'(a)',iostat=ios)'"table1_polynomial":"only_printed_indices_0_to_6_old_prose_upper_limit7_unresolved",'
  if(ios==0)write(u,'(a)',iostat=ios)'"uniform_reconstruction":"Table1_bulk_plus_nonrelativistic_common_Mn_kinetic_relativistic_electron_proton_Slater",'
  if(ios==0)write(u,'(a)',iostat=ios)'"uniform_reconstruction_role":"diagnostic_not_exported_as_inhomogeneous_crust",'
  if(ios==0)write(u,'(a)',iostat=ios)'"full_hartree_fock_recomputed":false,'
  if(ios==0)write(u,'(a)',iostat=ios)'"full_hf_missing_inputs":["B_IJ","C_IJ","D_IJ","spin_orbit_W"],'
  if(ios==0)write(u,'(a)',iostat=ios)'"missing_inputs_reference_doi":"10.1103/PhysRevC.5.1472",'
  if(ios==0)write(u,'(a)',iostat=ios)'"physical_matching_acceptance_asserted":false}'
  call finish_write(u,ios,status,'generation receipt')
 end subroutine
end module compose_eos_nv_original
