! Read successful campaign triples grouped by (cell,request), descending step.
! Input fields: cell request max_dlogp M_msun R_km Lambda Pc_MeV_fm3.
program summarize_layered_mesh
  use iso_fortran_env, only: real64, iostat_end
  implicit none
  integer, parameter :: wp=real64
  character(24), parameter :: requests(*)=[character(24)::'Pc25','Pc50','Pc100','M1.4']
  real(wp), parameter :: steps(*)=[.04_wp,.02_wp,.01_wp]
  character(24) :: cell,request,first_cell,first_request
  integer :: j,k,ios,count(4)
  character(24) :: report
  real(wp) :: step,obs(4,3),maximum(4,2,4),target_error
  call get_command_argument(1,report)
  if(report=='method-spread')then
    call method_spread()
    stop
  endif
  count=0
  maximum=0
  target_error=0
  do
    do j=1,3
      read(*,*,iostat=ios)cell,request,step,obs(:,j)
      if(j==1.and.ios==iostat_end)exit
      if(ios/=0)error stop 'incomplete or malformed mesh series'
      if(j==1)then
        first_cell=cell
        first_request=request
      endif
      if(cell/=first_cell.or.request/=first_request.or.step/=steps(j))error stop 'inconsistent mesh triple'
    enddo
    if(ios==iostat_end)exit
    k=findloc(requests,request,dim=1)
    if(k==0.or.any(obs<=0))error stop 'unknown request or invalid observables'
    do j=1,2
      maximum(:,j,k)=max(maximum(:,j,k),abs(obs(:,j)/obs(:,j+1)-1))
    enddo
    if(k==4)target_error=max(target_error,maxval(abs(obs(1,:)-1.4_wp)))
    count(k)=count(k)+1
  enddo
  write(*,'(a)') 'request,n,coarse_step,fine_step,max_rel_M,max_rel_R,max_rel_Lambda,max_rel_Pc'
  do k=1,size(requests)
    do j=1,2
      write(*,'(a,",",i0,6(",",es17.10e2))')trim(requests(k)),count(k),steps(j),steps(j+1),maximum(:,j,k)
    enddo
  enddo
  write(0,'(a,es12.5)')'maximum absolute target-mass residual (solar masses) = ',target_error
contains
  subroutine method_spread()
    character(80) :: core,group,input_core,input_group
    real(wp) :: r(5),lambda(5)
    integer :: i,ios
    write(*,'(a)') 'core_id,causal_interval_status,n_methods,min_R_km,max_R_km,spread_R_km,min_Lambda,max_Lambda,spread_Lambda'
    do
      do i=1,5
        read(*,*,iostat=ios)input_core,input_group,r(i),lambda(i)
        if(i==1.and.ios==iostat_end)return
        if(ios/=0)error stop 'missing matching method'
        if(i==1)then
          core=input_core
          group=input_group
        endif
        if(core/=input_core.or.group/=input_group)error stop 'inconsistent matching group'
      enddo
      write(*,'(a,",",a,",5",6(",",es20.12e2))')trim(core),trim(group),minval(r),maxval(r),maxval(r)-minval(r), &
          minval(lambda),maxval(lambda),maxval(lambda)-minval(lambda)
    enddo
  end subroutine
end program
