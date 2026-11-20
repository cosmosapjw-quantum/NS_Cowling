! Independent all-Fortran C10 Brent eigenfrequency verification.
! Smooth Schwarzschild interior control + negative counterexample.
program cowling_c10_brent_root_test
  use iso_c_binding, only: c_int,c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_dvode, only: cowling_native_setup, &
      cowling_native_scan,cowling_native_root,cowling_native_shot
  implicit none
  integer,parameter :: n=4097,nf=121
  real(c_double),parameter :: pi=3.1415926535897932384626433832795_c_double
  real(c_double),parameter :: ckm=299792.458_c_double
  real(c_double),parameter :: compactness=1.e-4_c_double, radius=10._c_double
  real(c_double) :: s(n),v(5,n),dv(5,n),f(nf),scan(nf)
  real(c_double) :: root,det,roots(2),exact,match_s,surface_s
  real(c_double) :: hstart,hc,h,z,zsurface,z0,r,ep,max_ds
  real(c_double) :: si(2),so(2),inner(2,2),outer(2,2),residual
  integer(c_int) :: setup_status,flags(nf-1),nfev,status,niter,shot_nfev
  integer :: i,idx,method,mf

  ! Negative control from the earlier C8 synthetic coefficients: an
  ! artificial determinant sign discontinuity near 1728 Hz is not a root.
  do i=1,n
    s(i)=real(i-1,c_double)/real(n-1,c_double)
    v(:,i)=[.1_c_double+12._c_double*s(i), &
      -.3_c_double+.08_c_double*s(i), &
       .05_c_double+.12_c_double*s(i), &
       .011_c_double+.004_c_double*s(i), &
       .12_c_double+.035_c_double*s(i)]
    dv(:,i)=[12._c_double,.08_c_double,.12_c_double, &
      .004_c_double,.035_c_double]
  enddo
  call cowling_native_setup(n,s,v,dv,setup_status)
  if(setup_status/=0)stop 10
  call cowling_native_root(1700._c_double,1800._c_double,2_c_int, &
    1.e-9_c_double,1.e-12_c_double,.005_c_double,10_c_int, &
    .5_c_double,1._c_double,1.e-6_c_double,root,det, &
    niter,nfev,status)
  print *, 'SYNTHETIC_FALSE_SIGN_FLIP_STATUS',status
  if(status/=33.or.root/=0._c_double.or.det/=0._c_double)stop 11

  ! Positive physical control: exact Schwarzschild constant-density
  ! interior. The reference weak-gravity f-mode is
  ! f0 = (c/2pi)*sqrt(2 compactness/R^2) for ell=2 (Kelvin limit).
  zsurface=sqrt(1._c_double-2._c_double*compactness)
  z0=sqrt(1._c_double-2._c_double*compactness*(.01_c_double/radius)**2)
  hc=log(2._c_double*zsurface/(3._c_double*zsurface-1._c_double))
  hstart=log(2._c_double*zsurface/(3._c_double*zsurface-z0))
  surface_s=sqrt(hstart-1.e-10_c_double)
  match_s=sqrt(hstart-.15_c_double*hc)
  max_ds=min(surface_s/200._c_double,1.e-6_c_double/(2._c_double*surface_s))
  do i=1,n
    s(i)=real(i-1,c_double)/real(n-1,c_double)*surface_s
    h=hstart-s(i)*s(i)
    z=zsurface*(3._c_double-2._c_double*exp(-h))
    r=radius*sqrt(max(0._c_double,(1._c_double-z*z)/(2._c_double*compactness)))
    ep=(3._c_double*zsurface-z)/2._c_double
    v(:,i)=[r,log(ep),-log(z), &
       compactness*r/(radius*radius*z*ep),1.e100_c_double]
  enddo
  do i=2,n-1
    dv(:,i)=(v(:,i+1)-v(:,i-1))/(s(i+1)-s(i-1))
  enddo
  dv(:,1)=(v(:,2)-v(:,1))/(s(2)-s(1))
  dv(:,n)=(v(:,n)-v(:,n-1))/(s(n)-s(n-1))
  call cowling_native_setup(n,s,v,dv,setup_status)
  if(setup_status/=0)then
    print *, 'UNIFORM_DENSITY_SETUP_FAIL',setup_status
    stop 12
  endif
  exact=ckm*sqrt(2._c_double*compactness/(radius*radius))/(2._c_double*pi)
  si=[0._c_double,match_s]
  so=[surface_s,match_s]
  do method=1,2
    mf=10
    if(method==2)mf=22
    call cowling_native_root(.99_c_double*exact,1.01_c_double*exact, &
      2_c_int,1.e-9_c_double,1.e-12_c_double,max_ds,int(mf,c_int), &
      match_s,surface_s,1.e-6_c_double,root,det,niter,nfev,status)
    if(status/=0)then
      print *, 'PHYSICAL_ROOT_STATUS',method,status
      stop 13
    endif
    if(.not.ieee_is_finite(root).or.niter<1.or.nfev<1)stop 14
    if(abs(root/exact-1._c_double)>5.e-4_c_double)stop 15
    if(abs(det)>1.e-5_c_double)stop 16
    call cowling_native_shot(2._c_double*pi*root/ckm, &
      2_c_int,1.e-9_c_double,1.e-12_c_double,max_ds, &
      int(mf,c_int),2_c_int,si,2_c_int,so, &
      inner,outer,residual,shot_nfev,status)
    if(status/=0.or.abs(det-residual)>1.e-12_c_double)stop 17
    roots(method)=root
    write(*,'("UNIFORM_STAR METHOD ",I1," F_HZ ",F18.9, &
      " RESID ",ES15.7," ITER ",I4)') method,root,det,niter
  enddo
  if(abs(roots(1)-roots(2))>0.05_c_double)stop 18
  print *, 'C10_NATIVE_ROOT_UNIFORM_PHYSICAL_AND_FALSE_FLIP_PASS'
end program cowling_c10_brent_root_test
