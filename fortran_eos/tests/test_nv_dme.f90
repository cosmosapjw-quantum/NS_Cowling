! Independent native-Fortran checks: decimal Table I, real128 monomial oracle,
! finite differences of the whole energy, symmetry, uniform chemical potentials.
program test_nv_dme
 use compose_eos_types,only:wp
 use compose_eos_nv_original,only:nv_original_source_t,load_nv_original_source
 use compose_eos_nv_dme
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 integer,parameter::qp=selected_real_kind(30,4000)
 type(nv_original_source_t)::src,badsrc
 type(nv_dme_coeff_t)::co,sw
 real(wp)::b(0:7,0:7),c(0:6,0:7),d(0:6,0:6)
 real(qp)::bq(0:7,0:7),cq(0:7,0:7),dq(0:7,0:7),vq(0:7,0:7)
 real(qp)::xq(7),xp(7),xm(7),xpp(7),xmm(7),qvalue,step,oracle,oref(6),xx,yy,hxy(2)
 integer,parameter::shifts(4)=[-2,-1,1,2],weights(4)=[1,-8,8,-1]
 real(wp)::energy,ul(2),bm(2),qg(2),ws,nn,np,tau(2),gr(2),jp,measure(4),jac(7)
 real(wp)::max_poly,max_gradient,max_hessian,max_local,max_uniform,e_expected,mu_expected(2)
 real(wp)::primal(6),twopi,finite_diff_a,diag_a,hh_step
 character(1024)::srcpath,tablepath,line,status,token
 character::kind
 integer::u,ios,i,j,k,l,nread,sample,ix,iy,ii,jj,checks,mode
 logical::seen_b(0:7,0:7),seen_c(0:6,0:7),seen_d(0:6,0:6)
 bq=0;cq=0;dq=0;vq=0;seen_b=.false.;seen_c=.false.;seen_d=.false.
 max_poly=0;max_gradient=0;max_hessian=0;max_local=0;max_uniform=0;checks=0
 call get_command_argument(1,srcpath);call get_command_argument(2,tablepath)
 mode=0
 call get_command_argument(3,token)
 if(len_trim(token)>0)then
  read(token,*,iostat=ios)mode
  if(ios/=0.or.mode<0.or.mode>1)error stop 'invalid test coefficient mode'
 endif
 if(len_trim(srcpath)==0.or.len_trim(tablepath)==0)error stop 'requires NV1973 source and NV1972 table TSV'
 call load_nv_original_source(trim(srcpath),src,status);call ok()
 call nv_dme_table1_coefficients(b,c,d)
 open(newunit=u,file=trim(tablepath),status='old',action='read',iostat=ios)
 if(ios/=0)error stop 'cannot read independent decimal table'
 nread=0
 do
  read(u,'(a)',iostat=ios)line
  if(ios<0)exit
  if(ios/=0)error stop 'table read failure'
  if(len_trim(line)==0)cycle
  if(line(1:1)=='#')cycle
  read(line,*,iostat=ios)kind,i,j,token
  if(ios/=0)error stop 'malformed table row'
  read(token,*,iostat=ios)qvalue
  if(ios/=0)error stop 'malformed decimal coefficient'
  select case(kind)
  case('B')
   if(i<0.or.i>7.or.j<0.or.j>7)error stop 'B index'
   if(seen_b(i,j))error stop 'B duplicate'
   seen_b(i,j)=.true.;bq(i,j)=qvalue
   if(b(i,j)/=real(qvalue,wp))error stop 'B source transcription mismatch'
  case('C')
   if(i<0.or.i>6.or.j<0.or.j>7)error stop 'C index'
   if(seen_c(i,j))error stop 'C duplicate'
   seen_c(i,j)=.true.;cq(i,j)=qvalue
   if(c(i,j)/=real(qvalue,wp))error stop 'C source transcription mismatch'
  case('D')
   if(i<0.or.i>6.or.j<0.or.j>6)error stop 'D index'
   if(seen_d(i,j))error stop 'D duplicate'
   seen_d(i,j)=.true.;dq(i,j)=qvalue
   if(d(i,j)/=real(qvalue,wp))error stop 'D source transcription mismatch'
  case default
   error stop 'unknown table coefficient'
  end select
  nread=nread+1
 enddo
 close(u)
 if(nread/=169.or..not.all(seen_b).or..not.all(seen_c).or..not.all(seen_d))error stop 'incomplete table'
 if(any(d/=transpose(d)))error stop 'D symmetry'
 vq(0:6,0:6)=real(src%vnm,qp)
 if(mode==1)then
  dq(0,3)=-6.62699e5_qp;dq(3,0)=dq(0,3)
 endif
 ! Exercise both dilute crust densities and the local nuclear interior domain.
 do ix=1,10
  nn=.0001_wp+real(ix-1,wp)/9*.1399_wp
  do iy=1,10
   np=.00005_wp+real(iy-1,wp)/9*.13995_wp
   call nv_dme_coefficients(src,nn,np,co,status,coefficient_mode=mode);call ok()
   call nv_dme_coefficients(src,np,nn,sw,status,coefficient_mode=mode);call ok()
   if(maxval(abs(co%b-sw%b(2:1:-1)))>1e-11_wp)error stop 'B species exchange'
   if(maxval(abs(co%c-sw%c(2:1:-1)))>1e-11_wp)error stop 'C species exchange'
   if(abs(co%d-sw%d)>2e-10_wp)error stop 'D exchange'
   do k=1,2
    if(k==1)then;xx=real(nn,qp);yy=real(np,qp)
    else;xx=real(np,qp);yy=real(nn,qp);endif
    call direct_poly(bq,xx,yy,oref)
    if(k==1)then
     primal=[co%b(k),co%db(:,k),co%hb(1,1,k),co%hb(1,2,k),co%hb(2,2,k)]
    else
     primal=[co%b(k),co%db(2,k),co%db(1,k),co%hb(2,2,k),co%hb(1,2,k),co%hb(1,1,k)]
    endif
    call compare_poly(primal,oref)
    call direct_poly(cq,xx,yy,oref)
    if(k==1)then
     primal=[co%c(k),co%dc(:,k),co%hc(1,1,k),co%hc(1,2,k),co%hc(2,2,k)]
    else
     primal=[co%c(k),co%dc(2,k),co%dc(1,k),co%hc(2,2,k),co%hc(1,2,k),co%hc(1,1,k)]
    endif
    call compare_poly(primal,oref)
   enddo
   call direct_poly(dq,real(nn,qp),real(np,qp),oref)
   primal=[co%d,co%dd,co%hd(1,1),co%hd(1,2),co%hd(2,2)]
   call compare_poly(primal,oref)
   tau=co%tau_uniform*[1.13_wp,.89_wp]
   gr=[-.003_wp,.0017_wp];jp=.00012_wp
   call nv_dme_local(src,nn,np,tau,gr,jp,energy,ul,bm,qg,ws,status,coefficient_mode=mode);call ok()
   xq=real([nn,np,tau,gr,jp],qp);jac=[ul,bm,qg,ws]
   oracle=energy_oracle(xq)
   max_local=max(max_local,real(abs(real(energy,qp)-oracle)/max(1._qp,abs(oracle)),wp))
   ! Independent real128 fourth-order finite differences of the complete EDF.
   do k=1,7
    step=1e-5_qp*max(abs(xq(k)),1e-5_qp)
    xp=xq;xm=xq;xpp=xq;xmm=xq
    xp(k)=xp(k)+step;xm(k)=xm(k)-step
    xpp(k)=xpp(k)+2*step;xmm(k)=xmm(k)-2*step
    oracle=(-energy_oracle(xpp)+8*energy_oracle(xp)-8*energy_oracle(xm)+energy_oracle(xmm))/(12*step)
    max_local=max(max_local,real(abs(real(jac(k),qp)-oracle)/max(1._qp,abs(oracle)),wp))
   enddo
   ! An independent numerical Hessian of A verifies the fractional-power terms.
   if(.not.co%a_hessian_defined)error stop 'interior A Hessian missing'
   do k=1,2
    step=1e-4_qp*xq(k)
    xp=xq;xm=xq;xpp=xq;xmm=xq
    xp(k)=xp(k)+step;xm(k)=xm(k)-step
    xpp(k)=xpp(k)+2*step;xmm(k)=xmm(k)-2*step
    oracle=(-a_oracle(xpp(1:2))+16*a_oracle(xp(1:2))-30*a_oracle(xq(1:2))+ &
      16*a_oracle(xm(1:2))-a_oracle(xmm(1:2)))/(12*step*step)
    max_hessian=max(max_hessian,real(abs(real(co%ha(k,k),qp)-oracle)/max(1._qp,abs(oracle)),wp))
   enddo
   oracle=0;hxy=1e-4_qp*xq(1:2)
   do ii=1,4
    do jj=1,4
     xp=xq
     xp(1:2)=xp(1:2)+hxy*[shifts(ii),shifts(jj)]
     oracle=oracle+weights(ii)*weights(jj)*a_oracle(xp(1:2))
    enddo
   enddo
   oracle=oracle/(144*hxy(1)*hxy(2))
   max_hessian=max(max_hessian,real(abs(real(co%ha(1,2),qp)-oracle)/max(1._qp,abs(oracle)),wp))
   if(co%ha(1,2)/=co%ha(2,1))error stop 'A Hessian not symmetric'
   ! Uniform cancellation is an independent physical identity, including mu.
   call nv_dme_local(src,nn,np,co%tau_uniform,[0._wp,0._wp],0._wp,energy,ul,bm,qg,ws,status,coefficient_mode=mode);call ok()
   e_expected=nv_dme_h2m*sum(co%tau_uniform)+co%vnm
   mu_expected=nv_dme_h2m*co%dtau_uniform+co%dvnm
   max_uniform=max(max_uniform,abs(energy-e_expected),maxval(abs(ul+bm*co%dtau_uniform-mu_expected)))
   checks=checks+1
  enddo
 enddo
 ! Vacuum and species-axis first derivatives remain finite, even where A'' does not.
 call nv_dme_local(src,0._wp,0._wp,[0._wp,0._wp],[0._wp,0._wp],0._wp,energy,ul,bm,qg,ws,status,coefficient_mode=mode);call ok()
 if(energy/=0.or.any(ul/=0))error stop 'vacuum limit'
 call nv_dme_coefficients(src,.01_wp,0._wp,co,status,coefficient_mode=mode);call ok()
 if(co%a_hessian_defined)error stop 'singular A Hessian incorrectly declared finite'
 call nv_dme_local(src,.01_wp,0._wp,co%tau_uniform,[0._wp,0._wp],0._wp,energy,ul,bm,qg,ws,status,coefficient_mode=mode);call ok()
 call nv_dme_local(src,-1e-8_wp,0._wp,[0._wp,0._wp],[0._wp,0._wp],0._wp,energy,ul,bm,qg,ws,status)
 if(status=='ok')error stop 'negative density accepted'
 call nv_dme_local(src,.150001_wp,0._wp,[0._wp,0._wp],[0._wp,0._wp],0._wp,energy,ul,bm,qg,ws,status)
 if(status=='ok')error stop 'out-of-source domain accepted'
 call nv_dme_local(src,.01_wp,.001_wp,[-1._wp,0._wp],[0._wp,0._wp],0._wp,energy,ul,bm,qg,ws,status)
 if(status=='ok')error stop 'negative kinetic density accepted'
 call nv_dme_local(badsrc,.01_wp,.001_wp,[0._wp,0._wp],[0._wp,0._wp],0._wp,energy,ul,bm,qg,ws,status)
 if(status=='ok')error stop 'unvalidated source accepted'
 call nv_dme_local(src,.01_wp,.001_wp,[0._wp,0._wp],[ieee_value(0._wp,ieee_quiet_nan),0._wp], &
  0._wp,energy,ul,bm,qg,ws,status)
 if(status=='ok')error stop 'nonfinite gradient accepted'
 call nv_dme_local(src,.01_wp,.001_wp,[0._wp,0._wp],[0._wp,0._wp], &
  0._wp,energy,ul,bm,qg,ws,status,coefficient_mode=2)
 if(status=='ok')error stop 'unknown coefficient mode accepted'
 write(*,'(a,i0)')'coefficient_mode=',mode
 write(*,'(a,i0)')'table_coefficients_checked=',nread
 write(*,'(a,i0)')'interior_density_points=',checks
 write(*,'(a,es24.15)')'max_polynomial_scaled_error=',max_poly
 write(*,'(a,es24.15)')'max_gradient_scaled_error=',max_gradient
 write(*,'(a,es24.15)')'max_hessian_scaled_error=',max_hessian
 write(*,'(a,es24.15)')'max_full_local_variation_scaled_error=',max_local
 write(*,'(a,es24.15)')'max_uniform_identity_absolute_error=',max_uniform
 if(max_poly>2e-10_wp.or.max_gradient>2e-9_wp.or.max_hessian>2e-8_wp)error stop 'polynomial oracle mismatch'
 if(max_local>2e-9_wp)error stop 'local variational derivative mismatch'
 if(max_uniform>2e-11_wp)error stop 'uniform energy or chemical potential mismatch'
 write(*,'(a)')'NV_DME_NATIVE_TESTS_PASS'
contains
 subroutine ok()
  if(status/='ok')then
   write(*,'(a)')trim(status)
   error stop 'unexpected source/functional error'
  endif
 end subroutine
 subroutine direct_poly(coeff,x,y,f)
  real(qp),intent(in)::coeff(0:7,0:7),x,y
  real(qp),intent(out)::f(6)
  real(qp)::a
  integer::i,j
  f=0
  do i=0,7
   do j=0,7
    a=coeff(i,j)
    f(1)=f(1)+a*x**i*y**j
    if(i>=1)f(2)=f(2)+i*a*x**(i-1)*y**j
    if(j>=1)f(3)=f(3)+j*a*x**i*y**(j-1)
    if(i>=2)f(4)=f(4)+i*(i-1)*a*x**(i-2)*y**j
    if(i>=1.and.j>=1)f(5)=f(5)+i*j*a*x**(i-1)*y**(j-1)
    if(j>=2)f(6)=f(6)+j*(j-1)*a*x**i*y**(j-2)
   enddo
  enddo
 end subroutine
 subroutine compare_poly(actual,expected)
  real(wp),intent(in)::actual(6)
  real(qp),intent(in)::expected(6)
  real(wp)::error(6)
  error=real(abs(real(actual,qp)-expected)/max(1._qp,abs(expected)),wp)
  max_poly=max(max_poly,error(1));max_gradient=max(max_gradient,maxval(error(2:3)))
  max_hessian=max(max_hessian,maxval(error(4:6)))
 end subroutine
 function value_poly(coeff,x,y) result(value)
  real(qp),intent(in)::coeff(0:7,0:7),x,y
  real(qp)::value
  integer::i,j
  value=0
  do i=0,7
   do j=0,7
    value=value+coeff(i,j)*x**i*y**j
   enddo
  enddo
 end function
 function a_oracle(rho) result(value)
  real(qp),intent(in)::rho(2)
  real(qp)::value,t(2),bb(2),qpi
  qpi=acos(-1._qp)
  t=.6_qp*(3*qpi*qpi)**(2._qp/3)*rho**(5._qp/3)
  bb=[value_poly(bq,rho(1),rho(2)),value_poly(bq,rho(2),rho(1))]
  value=value_poly(vq,rho(1),rho(2))-dot_product(bb,t)
 end function
 function energy_oracle(x) result(value)
  real(qp),intent(in)::x(7)
  real(qp)::value,bb(2),cc(2),dd
  bb=[value_poly(bq,x(1),x(2)),value_poly(bq,x(2),x(1))]
  cc=[value_poly(cq,x(1),x(2)),value_poly(cq,x(2),x(1))]
  dd=value_poly(dq,x(1),x(2))
  value=real(nv_dme_h2m,qp)*sum(x(3:4))+a_oracle(x(1:2))+dot_product(bb,x(3:4))+ &
   dot_product(cc,x(5:6)**2)+dd*x(5)*x(6)+60._qp*(x(5)+2*x(6))*x(7)
 end function
end program test_nv_dme
