program test_verify_domains_20261123
  use compose_eos_types
  use compose_eos_verify
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(eos_data_t)::e
  type(verify_report_t)::r
  integer::failures,j
  real(wp)::n,p,en,mu,b,c
  character(len=64)::mode
  failures=0;mode='all'
  if(command_argument_count()>0)call get_command_argument(1,mode)

  if(selected('valid_smooth'))then
    call source(e,32)
    call check_consistency(e,r)
    call require(r%overall,'valid_smooth')
  endif
  if(selected('negative_energy'))then
    call hidden_bad_phase(e,.true.)
    call check_consistency(e,r)
    call require(.not.r%overall,'negative_energy')
  endif
  if(selected('negative_internal_slope'))then
    call hidden_bad_phase(e,.false.)
    call check_consistency(e,r)
    call require(.not.r%overall.and.r%invalid_derivative_points==2,'negative_internal_slope')
  endif
  if(selected('nan_entropy'))then
    call source(e,32)
    e%q(2,17)=ieee_value(0._wp,ieee_quiet_nan)
    call check_consistency(e,r)
    call require(.not.r%overall,'nan_entropy')
  endif
  if(selected('zero_density'))then
    call source(e,32);e%nb(17)=0._wp
    call check_consistency(e,r)
    call require(.not.r%overall,'zero_density')
  endif
  if(selected('duplicate_density'))then
    call source(e,32);e%nb(17)=e%nb(16)
    call check_consistency(e,r)
    call require(.not.r%overall,'duplicate_density')
  endif
  if(selected('missing_q'))then
    call source(e,32);deallocate(e%q)
    call check_consistency(e,r)
    call require(.not.r%overall,'missing_q')
  endif
  if(selected('bad_q_shape'))then
    call source(e,32);deallocate(e%q);allocate(e%q(7,31));e%q=0._wp
    call check_consistency(e,r)
    call require(.not.r%overall,'bad_q_shape')
  endif
  if(selected('missing_phase'))then
    call source(e,32);deallocate(e%iphase)
    call check_consistency(e,r)
    call require(.not.r%overall,'missing_phase')
  endif
  if(selected('negative_pressure'))then
    call source(e,32);e%q(1,:)=-e%q(1,:)
    e%q(3,:)=e%q(6,:)+e%q(1,:)/e%mn
    call check_consistency(e,r)
    call require(.not.r%overall,'negative_pressure')
  endif
  if(selected('nan_metadata'))then
    call source(e,32);e%charge_fraction_grid=ieee_value(0._wp,ieee_quiet_nan)
    call check_consistency(e,r)
    call require(.not.r%overall,'nan_metadata')
  endif
  if(selected('invalid_tolerance'))then
    call source(e,32)
    call check_consistency(e,r,max_log_step=-1._wp)
    call require(.not.r%overall,'invalid_tolerance')
  endif
  if(selected('il0_effective_mu'))then
    call source(e,32);e%il=0;e%charge_fraction_grid=.2_wp
    e%q(4,:)=1000._wp*e%nb/e%mn
    e%q(3,:)=e%q(3,:)-e%charge_fraction_grid*e%q(4,:)
    call check_consistency(e,r)
    call require(r%overall.and.r%c7.and.r%c9,'il0_effective_mu')
  endif
  if(selected('valid_maxwell'))then
    call allocate_eos(e,6)
    e%mn=939._wp;e%mp=938._wp;e%il=1;e%beta_equilibrated=.true.
    e%temperature_mev=0._wp;e%charge_fraction_grid=0._wp
    e%nb=[.010_wp,.011_wp,.012_wp,.014_wp,.015_wp,.016_wp]
    e%iphase=0
    b=10._wp*(.012_wp**2-.014_wp**2)
    c=939._wp+20._wp*(.012_wp-.014_wp)
    do j=1,6
      n=e%nb(j)
      if(j<=3)then
        p=10._wp*n*n;en=939._wp*n+10._wp*n*n;mu=939._wp+20._wp*n
      else
        p=10._wp*n*n+b;en=c*n+10._wp*n*n-b;mu=c+20._wp*n
      endif
      e%q(1,j)=p/n;e%q(3,j)=mu/e%mn-1._wp
      e%q(6,j)=en/(n*e%mn)-1._wp;e%q(7,j)=e%q(6,j)
    enddo
    call check_consistency(e,r)
    call require(r%overall.and.r%thermodynamic_jump_boundaries==1.and. &
                 r%thermodynamic_jump_failures==0,'valid_maxwell')
    e%iphase(4:)=1
    call check_consistency(e,r)
    call require(r%overall,'valid_labelled_maxwell')
  endif
  if(selected('maxwell_near_endpoint'))then
    call allocate_eos(e,6)
    e%mn=939._wp;e%mp=938._wp;e%il=1;e%beta_equilibrated=.true.
    e%temperature_mev=0._wp;e%charge_fraction_grid=0._wp
    e%nb=[.010_wp,.011_wp,.014_wp,.015_wp,.016_wp,.017_wp]
    b=10._wp*(e%nb(2)**2-e%nb(3)**2)
    c=939._wp+20._wp*(e%nb(2)-e%nb(3))
    do j=1,6
      n=e%nb(j)
      if(j<=2)then
        p=10._wp*n*n;en=939._wp*n+10._wp*n*n;mu=939._wp+20._wp*n
      else
        p=10._wp*n*n+b;en=c*n+10._wp*n*n-b;mu=c+20._wp*n
      endif
      e%q(1,j)=p/n;e%q(3,j)=mu/e%mn-1._wp
      e%q(6,j)=en/(n*e%mn)-1._wp;e%q(7,j)=e%q(6,j)
    enddo
    call check_consistency(e,r)
    call require(r%overall.and.r%derivative_points==3.and.r%invalid_derivative_points==0, &
                 'maxwell_near_endpoint')
  endif
  print '(a,i0)','FAILURES=',failures
  if(failures>0)stop 1
contains
  logical function selected(name)
    character(len=*),intent(in)::name
    selected=trim(mode)=='all'.or.trim(mode)==name
  end function selected
  subroutine source(eos,k)
    type(eos_data_t),intent(out)::eos
    integer,intent(in)::k
    integer::i
    real(wp)::nb
    call allocate_eos(eos,k)
    eos%mn=939._wp;eos%mp=938._wp;eos%beta_equilibrated=.true.
    do i=1,k
      nb=.01_wp+real(i-1,wp)*.001_wp
      eos%nb(i)=nb;eos%q(1,i)=10._wp*nb
      eos%q(3,i)=20._wp*nb/eos%mn
      eos%q(6,i)=10._wp*nb/eos%mn;eos%q(7,i)=eos%q(6,i)
    enddo
  end subroutine source
  subroutine hidden_bad_phase(eos,negative)
    type(eos_data_t),intent(out)::eos
    logical,intent(in)::negative
    integer::i
    real(wp)::nb,pr,ed
    call source(eos,6)
    eos%nb=[.010_wp,.011_wp,.012_wp,.020_wp,.021_wp,.022_wp]
    eos%iphase(4:)=1
    do i=1,6
      nb=eos%nb(i)
      if(i<=3)then
        pr=10._wp*nb*nb;ed=939._wp*nb+pr
      else
        pr=1000._wp*nb*nb
        ed=.1_wp-nb
        if(negative)ed=-.1_wp-nb
      endif
      eos%q(1,i)=pr/nb;eos%q(3,i)=(ed+pr)/(nb*eos%mn)-1._wp
      eos%q(6,i)=ed/(nb*eos%mn)-1._wp;eos%q(7,i)=eos%q(6,i)
    enddo
  end subroutine hidden_bad_phase
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    print '(a,1x,l1,1x,a,1x,l1,1x,a,1x,a,1x,a,i0,1x,a,i0)', &
      trim(label),ok,'overall',r%overall,'failed',trim(r%failed_checks), &
      'checked=',r%derivative_points,'masked=',r%derivative_masked
    if(.not.ok)failures=failures+1
  end subroutine require
end program test_verify_domains_20261123
