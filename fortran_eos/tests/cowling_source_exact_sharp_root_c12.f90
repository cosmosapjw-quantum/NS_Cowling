! Source-locked independent Fortran scientific regression, no Python runtime.
! Reads exactly the original one-sided EOS laws and PCHIP geometry tables.
! Tests native source c_s² + native strong knot restarts + native DVODE,
! slow/rapid Maxwell transfer and complete Fortran bracketed polishing.
program cowling_c12_original_eos_sharp_root
  use iso_c_binding, only: c_double,c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_source_exact, only: cowling_native_source_eos_setup
  use compose_eos_cowling_sharp_source_root, only: cowling_native_sharp_source_brent
  implicit none
  character(len=1024)::input
  character(len=32)::repeats_arg
  integer::u,ios,no,nu,jump,nh,nl,i,j,mf,conversion,expectedcuts,nrepeats,trial
  real(c_double)::first_roots(2)
  real(c_double),allocatable::oe(:,:),ue(:,:),sh(:),sl(:)
  real(c_double),allocatable::vh(:,:),dh(:,:),vl(:,:),dl(:,:)
  real(c_double)::surface(4),bridge(7),geom(7),hstart,smatch
  real(c_double)::rtol,atol,max_ds,xtol,relative_jump
  real(c_double)::flo,fhi,fref,root,det,badgeom(7)
  integer(c_int)::st,iterations,calls,ncuts
  call get_command_argument(1,input)
  if(len_trim(input)==0)stop 10
  nrepeats=1
  call get_command_argument(2,repeats_arg)
  if(len_trim(repeats_arg)>0)then
    read(repeats_arg,*,iostat=ios)nrepeats
    if(ios/=0.or.nrepeats<1.or.nrepeats>1000)stop 9
  endif
  open(newunit=u,file=trim(input),action='read',status='old',iostat=ios)
  if(ios/=0)stop 11
  read(u,*,iostat=ios)no,nu,jump
  if(ios/=0.or.no<2.or.nu<4)stop 12
  allocate(oe(4,no),ue(3,nu))
  read(u,*,iostat=ios)surface
  if(ios/=0)stop 13
  read(u,*,iostat=ios)bridge
  if(ios/=0)stop 14
  do i=1,no
    read(u,*,iostat=ios)oe(:,i)
    if(ios/=0)stop 15
  enddo
  do i=1,nu
    read(u,*,iostat=ios)ue(:,i)
    if(ios/=0)stop 16
  enddo
  read(u,*,iostat=ios)nh,nl
  if(ios/=0.or.nh<4.or.nl<4)stop 17
  read(u,*,iostat=ios)hstart,smatch,rtol,atol,max_ds,xtol,relative_jump
  if(ios/=0)stop 18
  read(u,*,iostat=ios)geom
  if(ios/=0)stop 19
  allocate(sh(nh),sl(nl),vh(5,nh),vl(5,nl),dh(5,nh),dl(5,nl))
  do i=1,nh
    read(u,*,iostat=ios)sh(i),vh(:,i),dh(:,i)
    if(ios/=0)stop 20
  enddo
  do i=1,nl
    read(u,*,iostat=ios)sl(i),vl(:,i),dl(:,i)
    if(ios/=0)stop 21
  enddo
  read(u,*,iostat=ios)conversion,flo,fhi,fref,expectedcuts
  if(ios/=0)stop 22
  close(u)
  call cowling_native_source_eos_setup(int(no,c_int),oe,int(nu,c_int),ue, &
       surface,bridge,int(jump,c_int),st)
  if(st/=0)then
    print *, 'SOURCE_CONTRACT_FAIL',st
    stop 23
  endif
  do trial=1,nrepeats
  do j=1,2
    mf=10
    if(j==2)mf=22
    call cowling_native_sharp_source_brent(flo,fhi,2_c_int, &
       rtol,atol,max_ds,int(mf,c_int), &
       int(nh,c_int),sh,vh,dh,int(nl,c_int),sl,vl,dl, &
       hstart,smatch,geom,int(conversion,c_int),xtol,relative_jump, &
       root,det,iterations,calls,ncuts,st)
    if(st/=0)then
      print *, 'SHARP_SOURCE_ROOT_FAILED',mf,st
      stop 24
    endif
    if(.not.ieee_is_finite(root).or.root<=flo.or.root>=fhi)stop 25
    if(abs(root-fref)>.05_c_double)then
      print *, 'SHARP_SOURCE_ROOT_OUT_OF_REFERENCE',mf,root,fref
      stop 26
    endif
    if(abs(det)>1.e-5_c_double)stop 27
    if(iterations<1.or.calls<1)stop 28
    if(ncuts/=expectedcuts)stop 29
    if(trial==1)then
      first_roots(j)=root
      write(*,'("C12_NATIVE_SHARP_ROOT MF ",I2," CONVERSION ",I1, &
        " F_HZ ",F18.9," REF_DIFF_HZ ",ES13.5, &
        " DET ",ES13.5," KNOTS ",I3)') mf,conversion,root,root-fref,det,ncuts
    else
      if(root/=first_roots(j))stop 32
    endif
  enddo
  enddo
  ! Fail closed for degenerate bracket, invalid density jump and conversion.
  call cowling_native_sharp_source_brent(flo,flo,2_c_int, &
       rtol,atol,max_ds,10_c_int, &
       int(nh,c_int),sh,vh,dh,int(nl,c_int),sl,vl,dl, &
       hstart,smatch,geom,int(conversion,c_int),xtol,relative_jump, &
       root,det,iterations,calls,ncuts,st)
  if(st==0.or.root/=0._c_double.or.det/=0._c_double)stop 30
  badgeom=geom
  badgeom(7)=badgeom(6)*.9_c_double
  call cowling_native_sharp_source_brent(flo,fhi,2_c_int, &
       rtol,atol,max_ds,10_c_int, &
       int(nh,c_int),sh,vh,dh,int(nl,c_int),sl,vl,dl, &
       hstart,smatch,badgeom,int(conversion,c_int),xtol,relative_jump, &
       root,det,iterations,calls,ncuts,st)
  if(st==0.or.root/=0._c_double.or.det/=0._c_double)stop 31
  print *, 'C12_FORTRAN_SOURCE_EXACT_SHARP_ROOT_ADAMS_BDF_PASS'
  print *, 'C12_ROOT_REPEAT_DETERMINISTIC',nrepeats
end program cowling_c12_original_eos_sharp_root
