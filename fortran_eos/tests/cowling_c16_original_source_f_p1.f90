! C16 standalone native original-EOS f and p1 detector (NO PYTHON).
! Original high/low Maxwell stellar profiles and EOS contract come from
! source-locked C12 files. All scan points and root polishing are Fortran.
program cowling_c16_original_source_f_p1
  use iso_c_binding,only:c_double,c_int
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use compose_eos_cowling_source_exact,only:cowling_native_source_eos_setup
  use compose_eos_cowling_source_scan,only:cowling_native_sharp_source_scan
  use compose_eos_cowling_sharp_source_root,only:cowling_native_sharp_source_brent
  use compose_eos_cowling_mode_nodes,only:cowling_native_phase_nodes, &
       cowling_native_select_mode_branch
  use compose_eos_cowling_mode_overlap,only:cowling_native_mode_shape_overlap
  implicit none
  character(len=1024)::input
  character(len=64)::arg
  integer::u,ios,no,nu,jump,nh,nl,i,j,mf,conversion,expectedcuts
  integer::repeat_count,trial
  real(c_double),allocatable::oe(:,:),ue(:,:),sh(:),sl(:),vh(:,:),vl(:,:),dh(:,:),dl(:,:)
  real(c_double)::surface(4),bridge(7),geom(7),hstart,smatch,rtol,atol,dsmax
  real(c_double)::xtol,reljump,flo,fhi,p1ref,fref,freqs(4),dets(4)
  real(c_double)::root,det,first_roots(2),reference
  real(c_double)::rad(8),wr(8),v(8),total,inner,outer
  integer(c_int)::st,statuses(4),calls(4),source_count,iters,neval,ncuts
  integer(c_int)::phase(8),nodes_w,nodes_v,selected,nodestatus
  integer(c_int)::nodecountw(2),nodecountv(2),good(2)
  real(c_double)::cand(2)
  call get_command_argument(1,input)
  call get_command_argument(2,arg)
  if(len_trim(input)==0.or.len_trim(arg)==0)stop 10
  read(arg,*,iostat=ios)fref
  if(ios/=0.or.fref<=0._c_double)stop 11
  repeat_count=1
  call get_command_argument(3,arg)
  if(len_trim(arg)>0)then
    read(arg,*,iostat=ios)repeat_count
    if(ios/=0.or.repeat_count<1.or.repeat_count>100)stop 12
  endif
  open(newunit=u,file=trim(input),status='old',action='read',iostat=ios)
  if(ios/=0)stop 13
  read(u,*,iostat=ios)no,nu,jump
  if(ios/=0.or.no<2.or.nu<4)stop 14
  allocate(oe(4,no),ue(3,nu))
  read(u,*,iostat=ios)surface
  if(ios/=0)stop 15
  read(u,*,iostat=ios)bridge
  if(ios/=0)stop 16
  do i=1,no
    read(u,*,iostat=ios)oe(:,i)
    if(ios/=0)stop 17
  enddo
  do i=1,nu
    read(u,*,iostat=ios)ue(:,i)
    if(ios/=0)stop 18
  enddo
  read(u,*,iostat=ios)nh,nl
  if(ios/=0.or.nh<4.or.nl<4)stop 19
  read(u,*,iostat=ios)hstart,smatch,rtol,atol,dsmax,xtol,reljump
  if(ios/=0)stop 20
  read(u,*,iostat=ios)geom
  if(ios/=0)stop 21
  allocate(sh(nh),sl(nl),vh(5,nh),vl(5,nl),dh(5,nh),dl(5,nl))
  do i=1,nh
    read(u,*,iostat=ios)sh(i),vh(:,i),dh(:,i)
    if(ios/=0)stop 22
  enddo
  do i=1,nl
    read(u,*,iostat=ios)sl(i),vl(:,i),dl(:,i)
    if(ios/=0)stop 23
  enddo
  read(u,*,iostat=ios)conversion,flo,fhi,p1ref,expectedcuts
  if(ios/=0)stop 24
  close(u)
  call cowling_native_source_eos_setup(int(no,c_int),oe,int(nu,c_int),ue, &
    surface,bridge,int(jump,c_int),st)
  if(st/=0)stop 25
  freqs=[fref-40._c_double,fref+40._c_double, &
          p1ref-60._c_double,p1ref+60._c_double]
  if(any(freqs(2:4)<=freqs(1:3)))stop 26
  do trial=1,repeat_count
    call cowling_native_sharp_source_scan(4_c_int,freqs,2_c_int, &
       rtol,atol,dsmax,10_c_int, &
       int(nh,c_int),sh,vh,dh,int(nl,c_int),sl,vl,dl, &
       hstart,smatch,geom,int(conversion,c_int),reljump, &
       dets,calls,statuses,source_count,st)
    if(st/=0.or.any(statuses/=0).or.any(calls<=0))stop 27
    if(source_count/=expectedcuts)stop 28
    if(dets(1)*dets(2)>=0._c_double.or.dets(3)*dets(4)>=0._c_double)then
      print *,'C16_SCAN_NO_SIGN_CROSSING',dets
      stop 29
    endif
    do j=1,2
      if(j==1)then
        i=1;reference=fref
      else
        i=3;reference=p1ref
      endif
      call cowling_native_sharp_source_brent(freqs(i),freqs(i+1), &
        2_c_int,rtol,atol,dsmax,10_c_int, &
        int(nh,c_int),sh,vh,dh,int(nl,c_int),sl,vl,dl, &
        hstart,smatch,geom,int(conversion,c_int),xtol,reljump, &
        root,det,iters,neval,ncuts,st)
      if(st/=0)then
        print *,'C16_ROOT_FAIL',st
        stop 30
      endif
      if(abs(root-reference)>.02_c_double.or.abs(det)>1.e-5_c_double)stop 31
      if(iters<1.or.neval<1.or.ncuts/=expectedcuts)stop 32
      if(trial==1)then
        first_roots(j)=root
        write(*,'("C16_FORTRAN_F_P1 LABEL ",I1," HZ ",F16.9, &
          " REFERENCE_DIFF_HZ ",ES14.5," KNOTS ",I3)') j,root,root-reference,ncuts
      else
        if(root/=first_roots(j))stop 33
      endif
    enddo
  enddo
  ! Shape invariant tests are fully independent of Python input.
  phase=[0,0,0,0,1,1,1,1]
  rad=[0._c_double,.2_c_double,.4_c_double,.6_c_double, &
       .6_c_double,.8_c_double,1._c_double,1.2_c_double]
  wr=[0._c_double,.7_c_double,1._c_double,.4_c_double, &
      .4_c_double,-.1_c_double,-.4_c_double,-.2_c_double]
  v=wr
  call cowling_native_phase_nodes(8_c_int,phase,wr,wr,1.e-9_c_double, &
     nodes_w,nodes_v,nodestatus)
  if(nodestatus/=0.or.nodes_w<1.or.nodes_v<1)stop 34
  call cowling_native_mode_shape_overlap(8_c_int,phase,rad,wr, &
       8_c_int,phase,rad,wr,total,inner,outer,st)
  if(st/=0.or.abs(total-1._c_double)>1.e-12_c_double)stop 35
  cand=first_roots
  nodecountw=[0_c_int,1_c_int]
  nodecountv=[0_c_int,1_c_int]
  good=0
  call cowling_native_select_mode_branch(2_c_int,cand,nodecountw,nodecountv,good, &
      0_c_int,fref,.12_c_double,selected,st)
  if(st/=0.or.selected/=0)stop 36
  call cowling_native_select_mode_branch(2_c_int,cand,nodecountw,nodecountv,good, &
      1_c_int,p1ref,.12_c_double,selected,st)
  if(st/=0.or.selected/=1)stop 37
  print *, 'C16_NATIVE_ORIGINAL_SOURCE_SCAN_ROOT_SELECTOR_OVERLAP_PASS'
end program cowling_c16_original_source_f_p1
