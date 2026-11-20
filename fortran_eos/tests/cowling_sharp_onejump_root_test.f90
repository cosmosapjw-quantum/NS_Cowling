! Independent Fortran executable: source-locked sharp one-Maxwell interface.
! Python prepares an immutable one-sided PCHIP input file; this executable
! loads it and performs all native ODE segments, slow/rapid junction transfers
! and eigenfrequency polishing, without any Python numerical callback.
program cowling_sharp_onejump_brent_test
  use iso_c_binding, only: c_double,c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_sharp_root, only: cowling_native_sharp_onejump_brent
  implicit none
  character(len=1024)::input
  character(len=40)::selected_arg
  integer::unit,ios,nh,nl,i,ncase,j,mf,conversion,selected_case
  real(c_double)::smatch,smax,rtol,atol,max_ds,xtol
  real(c_double)::iface(7),bad(7),flo,fhi,fref,root,resid
  real(c_double),allocatable::sh(:),sl(:),vh(:,:),vl(:,:),dh(:,:),dl(:,:)
  integer(c_int)::steps,nfev,status

  call get_command_argument(1,input)
  if(len_trim(input)==0)stop 10
  selected_case=0
  call get_command_argument(2,selected_arg)
  if(len_trim(selected_arg)>0)then
    read(selected_arg,*,iostat=ios)selected_case
    if(ios/=0.or.selected_case<1.or.selected_case>4)stop 9
  endif
  open(newunit=unit,file=trim(input),status='old',action='read',iostat=ios)
  if(ios/=0)stop 11
  read(unit,*,iostat=ios)nh,nl
  if(ios/=0.or.nh<4.or.nl<4)stop 12
  allocate(sh(nh),sl(nl),vh(5,nh),vl(5,nl),dh(5,nh),dl(5,nl))
  read(unit,*,iostat=ios)smatch,smax,rtol,atol,max_ds,xtol
  if(ios/=0)stop 13
  read(unit,*,iostat=ios)iface
  if(ios/=0)stop 14
  do i=1,nh
    read(unit,*,iostat=ios)sh(i),vh(:,i),dh(:,i)
    if(ios/=0)stop 15
  enddo
  do i=1,nl
    read(unit,*,iostat=ios)sl(i),vl(:,i),dl(:,i)
    if(ios/=0)stop 16
  enddo
  read(unit,*,iostat=ios)ncase
  if(ios/=0.or.ncase/=4)stop 17
  if(abs(sl(nl)-smax)>1.e-12_c_double)stop 18
  do j=1,ncase
    read(unit,*,iostat=ios)conversion,flo,fhi,fref
    if(ios/=0)stop 19
    if(selected_case/=0.and.selected_case/=j)cycle
    do mf=10,22,12
      call cowling_native_sharp_onejump_brent(flo,fhi, &
         2_c_int,rtol,atol,max_ds,int(mf,c_int), &
         int(nh,c_int),sh,vh,dh,int(nl,c_int),sl,vl,dl, &
         smatch,iface,int(conversion,c_int),xtol, &
         root,resid,steps,nfev,status)
      if(status/=0)then
        print *, 'FORTRAN_SHARP_ROOT_FAILURE',j,conversion,mf,status
        stop 20
      endif
      if(.not.ieee_is_finite(root).or.root<=flo.or.root>=fhi)stop 21
      if(abs(root-fref)>.005_c_double)then
        print *, 'FORTRAN_FREQ_OUTSIDE_PCHIP_PARITY',root,fref
        stop 22
      endif
      if(abs(resid)>1.e-5_c_double)stop 23
      if(steps<1.or.nfev<1)stop 24
      write(*,'("FORTRAN_SHARP_CASE ",I2," CONV ",I1," MF ",I2, &
        " FREQ_HZ ",F18.9," REF_DIFF_HZ ",ES15.7, &
        " RESID ",ES15.7)') j,conversion,mf,root,root-fref,resid
    enddo
  enddo
  close(unit)
  if(selected_case/=0)then
    print *, 'SHARP_ONEJUMP_NATIVE_SELECTED_CASE_PASS',selected_case
  else
  ! Missing sign-change must not be promoted to a certified eigenfrequency.
  call cowling_native_sharp_onejump_brent(5600._c_double,5700._c_double, &
    2_c_int,rtol,atol,max_ds,10_c_int, &
    int(nh,c_int),sh,vh,dh,int(nl,c_int),sl,vl,dl, &
    smatch,iface,2_c_int,xtol,root,resid,steps,nfev,status)
  if(status/=3.or.root/=0._c_double.or.resid/=0._c_double)stop 25
  ! A nonphysical Maxwell latent heat must be rejected before DVODE starts.
  bad=iface
  bad(7)=bad(6)*0.9_c_double
  call cowling_native_sharp_onejump_brent(5110._c_double,5170._c_double, &
    2_c_int,rtol,atol,max_ds,10_c_int, &
    int(nh,c_int),sh,vh,dh,int(nl,c_int),sl,vl,dl, &
    smatch,bad,2_c_int,xtol,root,resid,steps,nfev,status)
  if(status/=1.or.root/=0._c_double.or.resid/=0._c_double)stop 26
  print *, 'SHARP_ONEJUMP_NATIVE_ADAMS_BDF_SLOW_RAPID_AND_BAD_INPUT_PASS'
  endif
end program cowling_sharp_onejump_brent_test
