program native_lower_crust_driver
 use compose_eos_types,only:wp,eos_data_t
 use compose_eos_lower_crust_native
 implicit none
 type(mass_catalogue_t)::cat
 type(lower_crust_state_t)::s
 type(lower_crust_state_t),allocatable::states(:),lower_states(:)
 type(bsk_anchor_t)::anchor
 type(eos_data_t)::eos,lower_eos
 type(lower_source_seam_t)::seam
 character(len=128)::action,arg,mode,status
 real(wp)::n,grid(17),res,jump,z,a,excess,radius
 integer::i,points,refinements,count
 call get_command_argument(1,action)
 call get_command_argument(2,arg)
 if(len_trim(arg)>0)then;read(arg,*)n;else;n=0;endif
 if(action=='outer'.or.action=='drip'.or.action=='sequence'.or.action=='source_segments'.or.action=='source_rows')then
  call load_mass_catalogue('data/outer_crust/hfb22_ame2016_cu.csv',cat,status);call check()
 endif
 select case(trim(action))
 case('isotope')
  call get_command_argument(3,arg);read(arg,*)z
  call get_command_argument(4,arg);read(arg,*)a
  call get_command_argument(5,arg);read(arg,*)excess
  call get_command_argument(6,arg);read(arg,*)radius
  call outer_isotope_state(n,z,a,excess,radius,s,status);call check()
 case('source_rows')
  call prepare_lower_source_native(cat,nint(n),17,eos,lower_eos,states,lower_states,seam,status);call check()
  do i=1,size(states)
   call emit(0,states(i))
  enddo
  do i=1,size(lower_states)
   call emit(1,lower_states(i))
  enddo
  stop
 case('source_segments')
  call prepare_lower_source_native(cat,nint(n),17,eos,lower_eos,states,lower_states,seam,status);call check()
  write(*,'(4(i0,1x),6(es25.16e3,1x))')eos%n,lower_eos%n,merge(1,0,seam%density_overlap),merge(1,0,seam%pressure_continuous.and.seam%mu_continuous),seam%outer_nb,seam%lower_nb,seam%delta_nb,seam%delta_pressure,seam%delta_mu_b,seam%delta_epsilon
  stop
 case('outer')
  call outer_state(cat,n,s,status);call check()
 case('drip')
  call outer_neutron_drip(cat,s,status);call check()
 case('bsk')
  call get_command_argument(3,mode)
  call bsk22_state(n,trim(mode),s,status);call check()
 case('anchor','carrier')
  call bsk22_make_anchor(.0002694320064010902_wp,.0004992059167689819_wp,939.5654205200001_wp,anchor,status,outer_epsilon=.2526497904790187_wp)
  call check()
  if(action=='anchor')then
   call bsk22_anchored_state(anchor,n,s,status,.true.);call check()
  else
   do i=1,17
    grid(i)=anchor%n0*(.01_wp/anchor%n0)**(real(i-1,wp)/16)
   enddo
   grid(17)=.01_wp
   call bsk22_generate(grid,'outer_anchored',eos,states,status,anchor,.true.);call check()
   write(*,'(4(i0,1x))')eos%n,merge(1,0,all(eos%yn_present)),merge(1,0,all(eos%yh_present)),merge(1,0,all(eos%row_status==2))
   stop
  endif
 case('sequence')
  points=nint(n)
  call outer_sequence(cat,points,eos,states,status,refinements=refinements);call check()
  res=0;jump=huge(1._wp)
  do i=2,size(states)
   if(states(i)%pressure==states(i-1)%pressure)then
    res=max(res,abs(states(i)%mu_b-states(i-1)%mu_b))
    jump=min(jump,states(i)%nb-states(i-1)%nb)
   endif
  enddo
  count=maxval(eos%iphase)-9
  write(*,'(2(i0,1x),2(es24.16e3,1x))')eos%n,count,res,jump
  stop
 case default
  error stop 'unknown native lower test action'
 end select
 write(*,'(15(es25.16e3,1x))')s%nb,s%epsilon,s%pressure,s%mu_b,s%mu_e,s%yn,s%yp,s%ye,s%ymu,s%yh,s%aav,s%zav,s%pressure_residual,s%quadrature_error,s%composition_mu_residual
contains
 subroutine emit(segment,row)
  integer,intent(in)::segment
  type(lower_crust_state_t),intent(in)::row
  write(*,'(i0,1x,15(es25.16e3,1x))')segment,row%nb,row%epsilon,row%pressure,row%mu_b,row%mu_e,row%yn,row%yp,row%ye,row%ymu,row%yh,row%aav,row%zav,row%pressure_residual,row%quadrature_error,row%composition_mu_residual
 end subroutine
 subroutine check()
  if(status/='ok')then;write(*,'(a)')trim(status);stop 2;endif
 end subroutine
end program
