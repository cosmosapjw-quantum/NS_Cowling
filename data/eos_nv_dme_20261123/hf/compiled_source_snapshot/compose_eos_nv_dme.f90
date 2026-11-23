! NV1972 Table I gradient/kinetic coefficients + NV1973 Table 1 bulk.
! Source: Phys.Rev.C5 (1972) 1472, DOI 10.1103/PhysRevC.5.1472,
! Eq(3.10)-(3.11), Table I p1487, Eq(6.1)-(6.4) p1485.
! Neutron spin-orbit is omitted variationally as prescribed on NV1973 p308:
! H_so = W0/2 * grad(rho_n+2*rho_p) dot J_p. No rest masses or Coulomb here.
! All operations are native Fortran; immutable coefficients permit OpenMP calls.
module compose_eos_nv_dme
 use compose_eos_types,only:wp
 use compose_eos_nv_original,only:nv_original_source_t,nv1973_bulk_potential
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 character(*),parameter,public::nv_dme_pdf_sha256= &
  'fd981fe5c9d7d914329c533fb57a2ac9375ed3d641f76ecc01b307bf82cd500b'
 real(wp),parameter,public::nv_dme_w0=120._wp
 integer,parameter,public::nv_dme_printed=0,nv_dme_d03_diagnostic=1
 character(*),parameter,public::nv_dme_d03_model='NV1973_DME_D03x100_diagnostic'
 ! This is a derived diagnostic repair, not an authenticated published erratum.
 ! The complete printed tables remain unchanged in nv_dme_table1_coefficients.
 real(wp),parameter::d03_repair=-6.62699e5_wp,d03_printed=-6.62699e3_wp
 ! Common kinetic mass convention matches the existing NV1973 bulk diagnostic.
 real(wp),parameter::mn=939.56542052_wp,hc=197.3269804_wp,pi=acos(-1._wp)
 real(wp),parameter,public::nv_dme_h2m=hc*hc/(2*mn)
 real(wp),parameter::tfactor=.6_wp*(3*pi*pi)**(2._wp/3)
 real(wp),parameter::btab(0:7,0:7)=transpose(reshape([ &
  0.0_wp, 3.18256e2_wp, -5.80513e3_wp, 7.23959e4_wp, &
  -5.12780e5_wp, 2.94816e6_wp, -1.57374e7_wp, 3.97592e7_wp, &
  1.21735e2_wp, -2.18424e3_wp, 6.92849e4_wp, -2.12527e5_wp, &
  -6.03181e6_wp, 2.00051e7_wp, 1.97485e8_wp, -6.60797e8_wp, &
  -4.13899e3_wp, 1.51287e4_wp, -8.00144e5_wp, -2.09223e6_wp, &
  1.03533e8_wp, -3.42514e8_wp, 1.03833e9_wp, -1.27850e10_wp, &
  7.62547e4_wp, -3.01156e5_wp, 1.05656e7_wp, -4.97705e6_wp, &
  6.84039e6_wp, -2.12738e9_wp, -1.65859e10_wp, 1.74936e11_wp, &
  -3.77066e5_wp, 5.00700e6_wp, -7.49536e7_wp, -3.82683e8_wp, &
  -9.27117e8_wp, -1.44377e9_wp, 1.18544e11_wp, 1.33994e11_wp, &
  -4.49404e6_wp, -3.56394e7_wp, 2.00503e8_wp, 4.53299e9_wp, &
  1.43914e10_wp, 1.69732e10_wp, -4.04503e11_wp, -7.43017e12_wp, &
  5.49780e7_wp, 8.87814e7_wp, -7.95987e7_wp, -9.21386e9_wp, &
  -1.67292e11_wp, 1.01556e12_wp, -2.21940e12_wp, 2.45505e13_wp, &
  -1.59317e8_wp, 0._wp, 0._wp, 0._wp, &
  0._wp, 0._wp, 0._wp, 0._wp ],[8,8]))

 real(wp),parameter::ctab(0:6,0:7)=transpose(reshape([ &
  1.30934e2_wp, -1.54662e3_wp, 3.05687e4_wp, -8.51161e4_wp, &
  -1.98427e6_wp, 8.08339e6_wp, 4.42696e7_wp, -1.65199e8_wp, &
  -4.74089e3_wp, 2.08118e4_wp, -4.64041e5_wp, -1.58197e6_wp, &
  5.66979e7_wp, -1.81885e8_wp, 5.19167e8_wp, -6.39250e9_wp, &
  1.04648e5_wp, -2.57851e5_wp, 6.49523e6_wp, 3.66425e6_wp, &
  -9.16113e5_wp, -1.59554e9_wp, -1.24394e10_wp, 1.31202e11_wp, &
  -5.66775e5_wp, 4.53068e6_wp, -6.26252e7_wp, -4.03940e8_wp, &
  -9.27117e8_wp, -1.44377e9_wp, 1.18544e11_wp, 1.33994e11_wp, &
  -1.07815e7_wp, -4.14671e7_wp, 2.25121e8_wp, 5.66624e9_wp, &
  1.79892e10_wp, 2.12165e10_wp, -5.05629e11_wp, -9.28771e12_wp, &
  1.44728e8_wp, 1.28070e8_wp, -1.19398e8_wp, -1.38208e10_wp, &
  -2.50938e11_wp, 1.52334e12_wp, -3.32910e12_wp, 3.68258e13_wp, &
  -4.75056e8_wp, 0._wp, 0._wp, 0._wp, &
  0._wp, 0._wp, 0._wp, 0._wp ],[8,7]))

 real(wp),parameter::dtab(0:6,0:6)=transpose(reshape([ &
  5.61303e2_wp, -8.45142e3_wp, 9.78215e4_wp, -6.62699e3_wp, &
  3.98432e6_wp, -2.69682e7_wp, 8.32712e7_wp, &
  -8.45142e3_wp, 1.48770e5_wp, -7.83358e5_wp, -4.55955e6_wp, &
  1.52685e7_wp, 3.45463e8_wp, -1.19619e9_wp, &
  9.78215e4_wp, -7.83358e5_wp, -5.28176e6_wp, 1.21991e8_wp, &
  -7.68296e8_wp, 4.95725e9_wp, -2.92841e10_wp, &
  -6.62699e3_wp, -4.55955e6_wp, 1.21991e8_wp, -6.47391e6_wp, &
  -3.58635e9_wp, -1.04875e10_wp, 1.38846e11_wp, &
  3.98432e6_wp, 1.52685e7_wp, -7.68296e8_wp, -3.58635e9_wp, &
  -3.60942e9_wp, 1.99033e11_wp, 1.50394e12_wp, &
  -2.69682e7_wp, 3.45463e8_wp, 4.95725e9_wp, -1.04875e10_wp, &
  1.99033e11_wp, -1.21351e12_wp, -1.63319e13_wp, &
  8.32712e7_wp, -1.19619e9_wp, -2.92841e10_wp, 1.38846e11_wp, &
  1.50394e12_wp, -1.63319e13_wp, 8.59269e13_wp ],[7,7]))
 type,public::nv_dme_coeff_t
  ! Species: 1=neutron, 2=proton. A derivative axis always indexes density.
  real(wp)::vnm=0,dvnm(2)=0,hvnm(2,2)=0
  real(wp)::a=0,da(2)=0,ha(2,2)=0
  real(wp)::b(2)=0,db(2,2)=0,hb(2,2,2)=0
  real(wp)::c(2)=0,dc(2,2)=0,hc(2,2,2)=0
  real(wp)::d=0,dd(2)=0,hd(2,2)=0
  real(wp)::tau_uniform(2)=0,dtau_uniform(2)=0
  ! A has a rho^(5/3) term: its second derivative can diverge on a density axis.
  ! ha is usable only when this flag is true. First derivatives remain defined.
  logical::a_hessian_defined=.false.
 end type
 public::nv_dme_coefficients,nv_dme_local,nv_dme_table1_coefficients
contains
 pure subroutine nv_dme_table1_coefficients(b,c,d)
  real(wp),intent(out)::b(0:7,0:7),c(0:6,0:7),d(0:6,0:6)
  b=btab;c=ctab;d=dtab
 end subroutine

 pure subroutine polynomial2(coef,x,y,value,gradient,hessian,need_hessian)
  real(wp),intent(in)::coef(0:,0:),x,y
  logical,intent(in)::need_hessian
  real(wp),intent(out)::value,gradient(2),hessian(2,2)
  real(wp)::a,b,c,oldv,oldy
  integer::i,j
  value=0;gradient=0;hessian=0
  ! Analytically differentiated nested Horner; no finite-difference derivatives.
  do i=ubound(coef,1),0,-1
   a=coef(i,ubound(coef,2));b=0;c=0
   do j=ubound(coef,2)-1,0,-1
    if(need_hessian)c=c*y+2*b
    b=b*y+a;a=a*y+coef(i,j)
   enddo
   oldv=value;oldy=gradient(2)
   if(need_hessian)hessian(1,1)=hessian(1,1)*x+2*gradient(1)
   gradient(1)=gradient(1)*x+oldv
   if(need_hessian)hessian(1,2)=hessian(1,2)*x+oldy
   gradient(2)=gradient(2)*x+b
   if(need_hessian)hessian(2,2)=hessian(2,2)*x+c
   value=value*x+a
  enddo
  hessian(2,1)=hessian(1,2)
 end subroutine

 subroutine nv_dme_coefficients(src,nn,np,co,status,need_hessian,coefficient_mode)
  type(nv_original_source_t),intent(in)::src
  real(wp),intent(in)::nn,np
  type(nv_dme_coeff_t),intent(out)::co
  logical,intent(in),optional::need_hessian
  integer,intent(in),optional::coefficient_mode
  logical::use_hessian
  character(*),intent(out)::status
  real(wp)::g(2),h(2,2),rho(2),root,ddtau(2)
  integer::q,i,j,mode
  co=nv_dme_coeff_t();use_hessian=.true.
  if(present(need_hessian))use_hessian=need_hessian
  mode=nv_dme_printed
  if(present(coefficient_mode))mode=coefficient_mode
  status='NV DME unknown coefficient mode'
  if(mode/=nv_dme_printed.and.mode/=nv_dme_d03_diagnostic)return
  ! Existing source loader pins the NV1973 source identity and domain.
  call nv1973_bulk_potential(src,nn,np,co%vnm,co%dvnm,co%hvnm,status)
  if(status/='ok')return
  call polynomial2(btab,nn,np,co%b(1),co%db(:,1),co%hb(:,:,1),use_hessian)
  call polynomial2(btab,np,nn,co%b(2),g,h,use_hessian)
  co%db(:,2)=[g(2),g(1)];co%hb(:,:,2)=h(2:1:-1,2:1:-1)
  call polynomial2(ctab,nn,np,co%c(1),co%dc(:,1),co%hc(:,:,1),use_hessian)
  call polynomial2(ctab,np,nn,co%c(2),g,h,use_hessian)
  co%dc(:,2)=[g(2),g(1)];co%hc(:,:,2)=h(2:1:-1,2:1:-1)
  call polynomial2(dtab,nn,np,co%d,co%dd,co%hd,use_hessian)
  if(mode==nv_dme_d03_diagnostic)then
   co%d=co%d+(d03_repair-d03_printed)*(nn**3+np**3)
   co%dd=co%dd+3*(d03_repair-d03_printed)*[nn*nn,np*np]
   if(use_hessian)then
    co%hd(1,1)=co%hd(1,1)+6*(d03_repair-d03_printed)*nn
    co%hd(2,2)=co%hd(2,2)+6*(d03_repair-d03_printed)*np
   endif
  endif
  rho=[nn,np];ddtau=0;co%a_hessian_defined=use_hessian
  do q=1,2
   root=rho(q)**(1._wp/3)
   co%tau_uniform(q)=tfactor*rho(q)*root*root
   co%dtau_uniform(q)=(5._wp/3)*tfactor*root*root
   if(rho(q)>0.and.use_hessian)then
    ddtau(q)=(10._wp/9)*tfactor/root
   else if(rho(q)==0.and.co%b(q)/=0)then
    ! rho=0 with the other species present: A_qq generally has no finite limit.
    ! Do not silently assert that zero-filled ha is the actual Hessian.
    co%a_hessian_defined=.false.
   endif
  enddo
  co%a=co%vnm-dot_product(co%b,co%tau_uniform)
  co%da=co%dvnm-matmul(co%db,co%tau_uniform)-co%b*co%dtau_uniform
  if(co%a_hessian_defined)then
   co%ha=co%hvnm
   do q=1,2
    co%ha=co%ha-co%hb(:,:,q)*co%tau_uniform(q)
    do i=1,2
     co%ha(i,q)=co%ha(i,q)-co%db(i,q)*co%dtau_uniform(q)
     co%ha(q,i)=co%ha(q,i)-co%db(i,q)*co%dtau_uniform(q)
    enddo
    co%ha(q,q)=co%ha(q,q)-co%b(q)*ddtau(q)
   enddo
  endif
  ! No coefficient clipping: finite domain errors are reported to the caller.
  do q=1,2
   if(.not.all(ieee_is_finite([co%b(q),co%c(q),co%db(:,q),co%dc(:,q)])))then
    status='NV DME nonfinite B/C polynomial';return
   endif
   do i=1,2
    do j=1,2
     if(.not.ieee_is_finite(co%hb(i,j,q)).or..not.ieee_is_finite(co%hc(i,j,q)))then
      status='NV DME nonfinite B/C Hessian';return
     endif
    enddo
   enddo
  enddo
  if(.not.all(ieee_is_finite([co%a,co%da,co%d,co%dd,reshape(co%hd,[4])])) )then
   status='NV DME nonfinite local functional';return
  endif
  if(co%a_hessian_defined)then
   if(.not.all(ieee_is_finite(co%ha)))then
    status='NV DME nonfinite A Hessian';return
   endif
  endif
  status='ok'
 end subroutine

 subroutine nv_dme_local(src,nn,np,tau,grad,jp,energy,un_local,bmass,qgrad,wspin,status,coefficient_mode)
  type(nv_original_source_t),intent(in)::src
  real(wp),intent(in)::nn,np,tau(2),grad(2),jp
  real(wp),intent(out)::energy,un_local(2),bmass(2),qgrad(2),wspin
  character(*),intent(out)::status
  integer,intent(in),optional::coefficient_mode
  type(nv_dme_coeff_t)::co
  real(wp)::delta_tau(2),cross
  energy=0;un_local=0;bmass=0;qgrad=0;wspin=0
  status='NV DME requires finite tau,grad,Jp and nonnegative tau'
  if(.not.all(ieee_is_finite([tau,grad,jp])))return
  if(any(tau<0))return
  call nv_dme_coefficients(src,nn,np,co,status,need_hessian=.false.,coefficient_mode=coefficient_mode)
  if(status/='ok')return
  bmass=nv_dme_h2m+co%b
  if(any(bmass<=0))then
   status='NV DME nonpositive inverse effective mass';return
  endif
  ! A+B*tau is evaluated as V+B*(tau-tau_TF) to preserve the uniform limit.
  delta_tau=tau-co%tau_uniform;cross=grad(1)*grad(2)
  wspin=.5_wp*nv_dme_w0*(grad(1)+2*grad(2))
  energy=nv_dme_h2m*sum(tau)+co%vnm+dot_product(co%b,delta_tau)+ &
   dot_product(co%c,grad*grad)+co%d*cross+wspin*jp
  ! Density derivatives hold tau,grad,Jp independent, as required in HF variation.
  un_local=co%dvnm+matmul(co%db,delta_tau)-co%b*co%dtau_uniform+ &
   matmul(co%dc,grad*grad)+co%dd*cross
  qgrad(1)=2*co%c(1)*grad(1)+co%d*grad(2)+.5_wp*nv_dme_w0*jp
  qgrad(2)=2*co%c(2)*grad(2)+co%d*grad(1)+nv_dme_w0*jp
  if(.not.all(ieee_is_finite([energy,un_local,bmass,qgrad,wspin])))then
   status='NV DME nonfinite local evaluation';return
  endif
  status='ok'
 end subroutine
end module compose_eos_nv_dme
