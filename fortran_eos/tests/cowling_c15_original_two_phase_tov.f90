! Independent source-locked original Maxwell TOV batch verification.
! Runs with NO Python interpreter: two physical EOS phases with a sharp
! epsilon jump and continuous pressure, mass and radius. The original
! piecewise EOS source sample arrays are already archived with SHA-256.
program cowling_c15_original_two_phase_tov
  use iso_c_binding, only:c_double,c_int
  use, intrinsic::ieee_arithmetic,only:ieee_is_finite
  use compose_eos_cowling_tov_two_phase_sweep,only:tov_native_two_phase_sweep
  implicit none
  character(len=1024)::path
  character(len=32)::repeats_arg
  integer::u,ios,nh,nl,nctrl,i,reps,trial
  real(c_double)::r0,m0,refm,refr,relative_r,relative_m
  real(c_double),allocatable::sh(:),ph(:),eh(:),sl(:),pl(:),el(:),bad_e(:)
  real(c_double),allocatable::rtols(:),atols(:),dsmax(:),result_h(:,:,:),result_l(:,:,:)
  real(c_double),allocatable::fixed_m(:),fixed_r(:)
  integer(c_int),allocatable::mfs(:),calls_h(:),calls_l(:),statuses(:)
  call get_command_argument(1,path)
  if(len_trim(path)==0)stop 10
  reps=1
  call get_command_argument(2,repeats_arg)
  if(len_trim(repeats_arg)>0)then
    read(repeats_arg,*,iostat=ios)reps
    if(ios/=0.or.reps<1.or.reps>500)stop 11
  endif
  open(newunit=u,file=trim(path),status="old",action="read",iostat=ios)
  if(ios/=0)stop 12
  read(u,*,iostat=ios)nh,nl,nctrl
  if(ios/=0.or.nh<4.or.nl<4.or.nctrl<4)stop 13
  read(u,*,iostat=ios)r0,m0,refm,refr
  if(ios/=0)stop 14
  allocate(sh(nh),ph(nh),eh(nh),sl(nl),pl(nl),el(nl),bad_e(nl))
  allocate(rtols(nctrl),atols(nctrl),dsmax(nctrl),mfs(nctrl))
  allocate(result_h(2,nh,nctrl),result_l(2,nl,nctrl))
  allocate(calls_h(nctrl),calls_l(nctrl),statuses(nctrl))
  allocate(fixed_m(nctrl),fixed_r(nctrl))
  do i=1,nh
    read(u,*,iostat=ios)sh(i),ph(i),eh(i)
    if(ios/=0)stop 15
  enddo
  do i=1,nl
    read(u,*,iostat=ios)sl(i),pl(i),el(i)
    if(ios/=0)stop 16
  enddo
  do i=1,nctrl
    read(u,*,iostat=ios)rtols(i),atols(i),dsmax(i),mfs(i)
    if(ios/=0)stop 17
  enddo
  close(u)
  if(eh(nh)<=el(1))stop 18
  if(abs(ph(nh)-pl(1))>1.e-8_c_double*max(ph(nh),pl(1)))stop 19
  do trial=1,reps
    call tov_native_two_phase_sweep(int(nctrl,c_int),int(nh,c_int),sh,ph,eh, &
      int(nl,c_int),sl,pl,el,r0,m0,rtols,atols,dsmax,mfs, &
      result_h,result_l,calls_h,calls_l,statuses)
    if(any(statuses/=[0,1,0,0]))then
      print *,'C15_TOV_STATUS_INVALID',statuses
      stop 20
    endif
    if(any(result_h(:,:,2)/=0._c_double).or.any(result_l(:,:,2)/=0._c_double))stop 21
    if(calls_h(2)/=0.or.calls_l(2)/=0)stop 22
    do i=1,nctrl
      if(statuses(i)/=0)cycle
      if(calls_h(i)<=0.or.calls_l(i)<=0)stop 23
      if(.not.ieee_is_finite(result_l(1,nl,i)))stop 24
      if(result_l(1,nl,i)<=2._c_double*result_l(2,nl,i))stop 25
      if(any(result_h(:,nh,i)/=result_l(:,1,i)))stop 26
      if(result_l(1,nl,i)<result_h(1,nh,i))stop 27
      if(result_l(2,nl,i)<result_h(2,nh,i))stop 28
      relative_r=abs(result_l(1,nl,i)/refr-1._c_double)
      relative_m=abs(result_l(2,nl,i)/refm-1._c_double)
      if(max(relative_r,relative_m)>2.e-4_c_double)then
        print *, 'C15_ORIGINAL_STAR_MASS_RADIUS_DRIFT',relative_m,relative_r
        stop 29
      endif
      if(trial==1)then
        fixed_r(i)=result_l(1,nl,i)
        fixed_m(i)=result_l(2,nl,i)
        write(*,'("C15_TWO_PHASE_TOV_RUN METHOD ",I2," RTOL ",ES12.4, &
         " RADIUS_KM ",F16.10," MASS_GEOM_KM ",F16.10, &
         " NFEV ",I8," SOURCE_CONTINUITY 1")') &
            mfs(i),rtols(i),fixed_r(i),fixed_m(i),calls_h(i)+calls_l(i)
      else
        if(result_l(1,nl,i)/=fixed_r(i).or.result_l(2,nl,i)/=fixed_m(i))stop 30
      endif
    enddo
  enddo
  bad_e=el
  bad_e(1)=2._c_double*eh(nh)
  call tov_native_two_phase_sweep(int(nctrl,c_int),int(nh,c_int),sh,ph,eh, &
      int(nl,c_int),sl,pl,bad_e,r0,m0,rtols,atols,dsmax,mfs, &
      result_h,result_l,calls_h,calls_l,statuses)
  if(any(statuses==0))stop 31
  if(any(result_h/=0._c_double).or.any(result_l/=0._c_double))stop 32
  print '("C15_STANDALONE_ORIGINAL_SOURCE_FORTRAN_TWO_PHASE_TOV_PASS REPEATS ",I4)',reps
end program cowling_c15_original_two_phase_tov
