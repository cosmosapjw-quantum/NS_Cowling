program test_core_subsaturation
 use compose_eos_types, only: wp,eos_data_t
 use compose_eos_core_native, only: core_state_t,core_model_bounds,solve_core_state,generate_core_native
 implicit none
 character(len=12),parameter::models(3)=[character(len=12)::'GM_GM1','RG_KDE0v1','RG_SLY4']
 real(wp),parameter::grid(4)=[.02_wp,.04_wp,.06_wp,.08_wp]
 type(core_state_t)::s,stencil(0:4)
 type(core_state_t),allocatable::states(:)
 type(eos_data_t)::eos
 character(len=240)::status
 real(wp)::lo,hi,mn,mp,n,h,de,errors(2)
 integer::i,j,k,l
 write(*,'(a)')'model,n_fm3,epsilon_mev_fm3,pressure_mev_fm3,mu_mev,firstlaw_abs_h,firstlaw_abs_half_h'
 do i=1,3
  call core_model_bounds(trim(models(i)),lo,hi,mn,mp,status)
  call require(status=='ok','declared model')
  call require(lo==.02_wp,'computational lower bound is .02')
  call generate_core_native(trim(models(i)),grid,eos,states,status)
  call require(status=='ok','native subsaturation generation: '//trim(status))
  call require(eos%n==4.and.all(eos%nb==grid),'requested source grid retained')
  do j=1,3
   n=grid(j)
   call solve_core_state(trim(models(i)),n,s,status)
   call require(status=='ok','central source state')
   call require(abs(s%yp-s%ye-s%ymu)<1.e-10_wp,'charge neutrality')
   call require(abs(s%mu_n-s%mu_p-s%mu_e)<1.e-8_wp,'beta equilibrium')
   call require(abs(s%epsilon+s%pressure-n*s%mu_n)<1.e-8_wp,'Euler identity')
   do l=1,2
    h=n*1.e-4_wp/real(l,wp)
    if(j==1)then
     ! Fourth-order one-sided stencil remains within declared support at .02.
     do k=0,4
      call solve_core_state(trim(models(i)),n+k*h,stencil(k),status)
      call require(status=='ok','forward derivative state')
     enddo
     de=(-25*stencil(0)%epsilon+48*stencil(1)%epsilon-36*stencil(2)%epsilon+ &
         16*stencil(3)%epsilon-3*stencil(4)%epsilon)/(12*h)
    else
     do k=0,4
      call solve_core_state(trim(models(i)),n+(k-2)*h,stencil(k),status)
      call require(status=='ok','central derivative state')
     enddo
     de=(stencil(0)%epsilon-8*stencil(1)%epsilon+8*stencil(3)%epsilon-stencil(4)%epsilon)/(12*h)
    endif
    errors(l)=abs(de-s%mu_n)
    call require(errors(l)<5.e-5_wp,'independent subsaturation first law')
   enddo
   write(*,'(a,6(a,es24.16e3))')trim(models(i)),',',n,',',s%epsilon,',',s%pressure,',',s%mu_n,',',errors(1),',',errors(2)
  enddo
 enddo
 call core_model_bounds('APR_APR',lo,hi,mn,mp,status)
 call require(lo==.085_wp,'APR declared bounds unchanged')
 call core_model_bounds('BBB_BHF2',lo,hi,mn,mp,status)
 call require(lo==.08_wp,'BBB source bounds unchanged')
 call solve_core_state('BBB_BHF2',.06_wp,s,status)
 call require(status/='ok','BBB energy table is not extrapolated')
 write(*,'(a)')'PASS native subsaturation domain and independent first-law checks'
contains
 subroutine require(condition,label)
  logical,intent(in)::condition
  character(len=*),intent(in)::label
  if(.not.condition)then
   write(*,'(a)')'FAIL '//trim(label)
   stop 1
  endif
 end subroutine
end program
