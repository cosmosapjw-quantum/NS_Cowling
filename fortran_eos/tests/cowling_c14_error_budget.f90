! C14 original-source controlled error budget, standalone Fortran.
! Three independent axes: retained stellar geometry knots (factor 1/2/4/8),
! DVODE relative tolerance and DVODE maximum step (units: s=sqrt(h0-h)).
! All source c_s², PCHIP slopes, native roots and Maxwell traction transfer
! are calculated in Fortran. Strong original EOS knots are pinned in
! geometry coarsening; no boundary interpolation across a Maxwell jump.
program cowling_c14_error_budget
  use iso_c_binding, only: c_double,c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use compose_eos_cowling_source_exact, only: cowling_native_source_eos_setup
  use compose_eos_cowling_pchip, only: cowling_native_pchip_slopes
  use compose_eos_cowling_error_budget, only: cowling_native_sharp_source_sweep
  implicit none
  integer :: u,ios,no,nu,jump,nh,nl,conversion,oldcuts,i,fac
  real(c_double),allocatable :: eos_outer(:,:),eos_upper(:,:)
  real(c_double),allocatable :: shi(:),slo(:),vhi(:,:),vlo(:,:),dhi(:,:),dlo(:,:)
  real(c_double),allocatable :: hs(:),ls(:),hv(:,:),lv(:,:),hd(:,:),ld(:,:)
  real(c_double),allocatable :: allknots(:)
  real(c_double)::surface(4),bridge(7),geom(7),hstart,smatch,rtol,atol
  real(c_double)::max_ds,xtol,relative_jump,flo,fhi,fref
  integer(c_int)::status
  character(len=1024)::path
  call get_command_argument(1,path)
  if(len_trim(path)==0)stop 10
  open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
  if(ios/=0)stop 11
  read(u,*,iostat=ios)no,nu,jump
  if(ios/=0.or.no<2.or.nu<4)stop 12
  allocate(eos_outer(4,no),eos_upper(3,nu))
  read(u,*,iostat=ios)surface
  if(ios/=0)stop 13
  read(u,*,iostat=ios)bridge
  if(ios/=0)stop 14
  do i=1,no
    read(u,*,iostat=ios)eos_outer(:,i)
    if(ios/=0)stop 15
  enddo
  do i=1,nu
    read(u,*,iostat=ios)eos_upper(:,i)
    if(ios/=0)stop 16
  enddo
  read(u,*,iostat=ios)nh,nl
  if(ios/=0.or.nh<4.or.nl<4)stop 17
  read(u,*,iostat=ios)hstart,smatch,rtol,atol,max_ds,xtol,relative_jump
  if(ios/=0)stop 18
  read(u,*,iostat=ios)geom
  if(ios/=0)stop 19
  allocate(shi(nh),slo(nl),vhi(5,nh),vlo(5,nl),dhi(5,nh),dlo(5,nl))
  do i=1,nh
    read(u,*,iostat=ios)shi(i),vhi(:,i),dhi(:,i)
    if(ios/=0)stop 20
  enddo
  do i=1,nl
    read(u,*,iostat=ios)slo(i),vlo(:,i),dlo(:,i)
    if(ios/=0)stop 21
  enddo
  read(u,*,iostat=ios)conversion,flo,fhi,fref,oldcuts
  if(ios/=0)stop 22
  close(u)
  call cowling_native_source_eos_setup(int(no,c_int),eos_outer, &
    int(nu,c_int),eos_upper,surface,bridge,int(jump,c_int),status)
  if(status/=0)then
    print *, 'C14_EOS_SETUP_FAILED',status
    stop 23
  endif
  allocate(allknots(no+nu+3))
  call build_source_knots()
  print *, 'C14_BEGIN',no,nu,nh,nl,conversion,flo,fhi,fref,oldcuts
  ! Geometry convergence: factors 1/2/4/8, numerical controls fixed.
  do fac=1,8
    if(fac/=1.and.fac/=2.and.fac/=4.and.fac/=8)cycle
    call geometry(fac)
    call run_controls("GEOMETRY",fac)
    deallocate(hs,ls,hv,lv,hd,ld)
  enddo
  ! At exactly one fixed full geometry, vary the ODE tolerances and maximum
  ! s-step *separately* while keeping all physical/source parameters fixed.
  call geometry(1)
  call run_controls("RTOL",1)
  call run_controls("MAXSTEP",1)
  print *, 'C14_FORTRAN_3AXIS_ERROR_BUDGET_COMPLETE'
contains
  subroutine build_source_knots()
    integer::j
    real(c_double)::rest,mu
    rest=surface(1)
    allknots(1)=sqrt(max(0._c_double,hstart - &
       log((rest+surface(4)/(surface(4)-1._c_double)*surface(3)/surface(2))/rest)))
    do j=1,no
      mu=eos_outer(3,j)+eos_outer(2,j)/eos_outer(1,j)
      allknots(j+1)=sqrt(max(0._c_double,hstart-log(mu/rest)))
    enddo
    do j=1,nu
      mu=(eos_upper(3,j)+eos_upper(2,j))/eos_upper(1,j)
      allknots(no+1+j)=sqrt(max(0._c_double,hstart-log(mu/rest)))
    enddo
    allknots(no+nu+2)=0._c_double
    allknots(no+nu+3)=slo(nl)
  end subroutine

  ! Retain one physical source geometry phase, one-sided endpoint, every
  ! integer stride and all original source chemical-potential boundary knots
  ! found on the immutable high-resolution grid. No new interpolated source
  ! EOS state is created by this coarsening experiment.
  subroutine phase_subsample(stride,s,y,orig_d,so,yo,deriv)
    integer,intent(in)::stride
    real(c_double),intent(in)::s(:),y(:,:),orig_d(:,:)
    real(c_double),allocatable,intent(out)::so(:),yo(:,:),deriv(:,:)
    logical,allocatable :: keep(:)
    real(c_double)::target,err
    integer::n,nkeep,j,k,p,where
    integer(c_int)::st
    n=size(s)
    allocate(keep(n))
    keep=.false.
    keep(1)=.true.;keep(n)=.true.
    do j=1,n,stride
      keep(j)=.true.
    enddo
    do j=1,size(allknots)
      target=allknots(j)
      if(target<=s(1).or.target>=s(n))cycle
      where=minloc(abs(s-target),dim=1)
      err=abs(s(where)-target)
      if(err<=1.e-9_c_double)keep(where)=.true.
    enddo
    nkeep=count(keep)
    allocate(so(nkeep),yo(5,nkeep),deriv(5,nkeep))
    k=0
    do j=1,n
      if(.not.keep(j))cycle
      k=k+1
      so(k)=s(j);yo(:,k)=y(:,j)
    enddo
    if(k/=nkeep)stop 30
    call cowling_native_pchip_slopes(int(nkeep,c_int),so,yo,deriv,st)
    if(st/=0)then
      print *, 'C14_PCHIP_COEFFICIENTS_REJECTED',stride,st
      stop 31
    endif
    if(stride==1)then
      if(any(abs(deriv-orig_d)>5.e-9_c_double+2.e-10_c_double*abs(orig_d)))then
        print *, 'C14_FULL_GRID_PCHIP_NOT_EQUIVALENT'
        stop 32
      endif
    endif
  end subroutine

  subroutine geometry(stride)
    integer,intent(in)::stride
    call phase_subsample(stride,shi,vhi,dhi,hs,hv,hd)
    call phase_subsample(stride,slo,vlo,dlo,ls,lv,ld)
  end subroutine

  subroutine run_controls(axis,factor)
    character(len=*),intent(in)::axis
    integer,intent(in)::factor
    real(c_double)::tols(8),steps(8),roots(8),resid(8)
    integer(c_int)::methods(8),iters(8),calls(8),cuts(8),stat(8)
    integer::n,j,m
    n=0
    if(axis=="GEOMETRY")then
      n=2
      tols(1:n)=rtol
      steps(1:n)=max_ds
      methods(1:2)=[10_c_int,22_c_int]
    else if(axis=="RTOL")then
      n=8
      do j=1,4
        select case(j)
          case(1);tols(2*j-1)=3.e-8_c_double
          case(2);tols(2*j-1)=rtol
          case(3);tols(2*j-1)=3.e-10_c_double
          case(4);tols(2*j-1)=1.e-10_c_double
        end select
        tols(2*j)=tols(2*j-1)
        steps(2*j-1:2*j)=max_ds
        methods(2*j-1:2*j)=[10_c_int,22_c_int]
      enddo
    else if(axis=="MAXSTEP")then
      n=8
      do j=1,4
        tols(2*j-1:2*j)=rtol
        select case(j)
          case(1);steps(2*j-1)=2._c_double*max_ds
          case(2);steps(2*j-1)=max_ds
          case(3);steps(2*j-1)=max_ds/2._c_double
          case(4);steps(2*j-1)=max_ds/4._c_double
        end select
        steps(2*j)=steps(2*j-1)
        methods(2*j-1:2*j)=[10_c_int,22_c_int]
      enddo
    else
      stop 40
    endif
    call cowling_native_sharp_source_sweep(int(n,c_int),flo,fhi,2_c_int, &
      atol,hstart,smatch,geom,int(conversion,c_int),xtol,relative_jump, &
      int(size(hs),c_int),hs,hv,hd,int(size(ls),c_int),ls,lv,ld, &
      tols(1:n),steps(1:n),methods(1:n),roots(1:n),resid(1:n), &
      iters(1:n),calls(1:n),cuts(1:n),stat(1:n))
    do j=1,n
      print *, 'C14_ROW',axis,factor,size(hs),size(ls),methods(j), &
        tols(j),steps(j),roots(j),resid(j),iters(j),calls(j),cuts(j),stat(j)
      if(stat(j)==0_c_int)then
        if(abs(resid(j))>1.e-5_c_double)stop 41
        if(cuts(j)/=oldcuts)stop 42
        if(.not.ieee_is_finite(roots(j)))stop 43
        if(abs(roots(j)-fref)>.1_c_double)then
          print *, 'C14_EXCESSIVE_FREQ_DEVIATION',axis,factor,roots(j),fref
          stop 44
        endif
      endif
    enddo
  end subroutine
end program cowling_c14_error_budget
