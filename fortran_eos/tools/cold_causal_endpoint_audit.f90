! Audit a sampled one-parameter, cold EOS using only n, epsilon and P.
! CSV input: nb_fm3,epsilon_MeV_fm3,pressure_MeV_fm3,Iphase,region.
! Frozen-composition h=(epsilon+P)/n is along that same compression path.
program cold_causal_endpoint_audit
  use iso_fortran_env, only: real64, real128, iostat_end, output_unit
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  character(*),parameter :: numeric='0123456789EeDd.-+'
  integer, parameter :: q=real128
  ! Input CSV values are rounded binary64 data. Evaluate integrals in quad and
  ! permit conservative absolute uncertainty in each subtracted physical value.
  real(q), parameter :: roundoff=128._q*real(epsilon(1._real64),q)
  integer, parameter :: bad_order=1, bad_ratio=2, bad_pressure_area=3
  integer, parameter :: compatible=0, bad_secant=4, equal_support=5
  type :: point
    real(q) :: n,energy,p
    character(24) :: region
  end type
  character(4096) :: path,line,detail_path
  character(128) :: option
  type(point) :: left,right
  integer :: first_input,detail_unit,file_idx,io,unit,rows,code,bad(4),first_bad,cat,equal_pairs,secant_bad,counts(4),why(3)
  real(q) :: dv,maximum_violation,low,high,mbias,worst_relative,uq
  logical :: ordinary_secant_ok
  uq=roundoff
  detail_path=""
  detail_unit=output_unit
  call get_command_argument(1,path)
  if (path=='--self-test') then
    call self_test()
    stop
  end if
  first_input=1
  do while(first_input<=command_argument_count())
    call get_command_argument(first_input,option,status=io)
    if(index(option,'--')/=1) exit
    if(first_input==command_argument_count()) error stop 'option requires a value'
    select case(option)
    case('--relative-uncertainty')
      call get_command_argument(first_input+1,option,status=io)
      if(io/=0.or.len_trim(option)==0.or.verify(trim(option),numeric)/=0) error stop 'invalid relative uncertainty'
      read(option,*,iostat=io)uq
      if(io/=0) error stop 'invalid relative uncertainty'
      if(.not.ieee_is_finite(uq).or.uq<0.or.uq>=1) error stop 'need 0 <= relative uncertainty < 1'
      uq=max(uq,roundoff)
    case('--intervals')
      call get_command_argument(first_input+1,detail_path,status=io)
      if(io/=0.or.detail_path=='') error stop 'invalid interval report path'
    case default
      error stop 'unrecognized option'
    end select
    first_input=first_input+2
  end do
  if(first_input>command_argument_count()) error stop 'usage: cold_causal_endpoint_audit [--relative-uncertainty R] [--intervals FILE] <CSV> ...'
  if(detail_path/='') then
    open(newunit=detail_unit,file=trim(detail_path),status='new',action='write')
    write(detail_unit,'(a)')'input,left_row,region_left,region_right,reason,n_left,n_right,epsilon_left,epsilon_right,'// &
      'P_left,P_right,h_left,h_right,lower_P_area_bound,upper_P_area_bound,guaranteed_pressure_excess,relative_excess'
  end if
  print '(a)','input,rows,crust_pairs,bridge_pairs,core_pairs,other_pairs,bad_crust,bad_bridge,bad_core,bad_other,'// &
    'order_fail,ratio_fail,pressure_area_fail,secant_fail,near_equal_P_h,first_bad_left_row,max_pressure_area_excess_MeV_fm3,max_relative_excess'
  do file_idx=first_input,command_argument_count()
    call get_command_argument(file_idx,path,status=io)
    if(io/=0) error stop 'path length'
    open(newunit=unit,file=trim(path),status='old',action='read')
    read(unit,'(a)')line
    if(trim(line)/='nb_fm3,epsilon_MeV_fm3,pressure_MeV_fm3,Iphase,region') error stop 'expected physical CSV schema'
    rows=0; bad=0; first_bad=0; equal_pairs=0; secant_bad=0; counts=0; why=0; maximum_violation=0;worst_relative=0
    do
      read(unit,'(a)',iostat=io)line
      if(io==iostat_end) exit
      if(io/=0) error stop 'CSV stream read failed'
      call read_point(line,right)
      if(.not.all(ieee_is_finite([right%n,right%energy,right%p]))) error stop 'nonfinite EOS value'
      if(right%n<=0.or.right%energy<=0.or.right%p<0) error stop 'nonphysical EOS input'
      if(rows>0) then
        cat=4
        if(left%region==right%region) then
          select case(left%region)
          case('crust'); cat=1
          case('interface'); cat=2
          case('core'); cat=3
          end select
        end if
        counts(cat)=counts(cat)+1
        call check_pair(left,right,code,dv,ordinary_secant_ok,low,high)
        if(.not.ordinary_secant_ok) secant_bad=secant_bad+1
        if(code==equal_support) then
          equal_pairs=equal_pairs+1
        else if(code/=compatible) then
          if(first_bad==0) first_bad=rows
          bad(cat)=bad(cat)+1
          if(code<=size(why)) why(code)=why(code)+1
          maximum_violation=max(maximum_violation,dv)
          mbias=dv/max(right%p,left%p,tiny(1._q))
          worst_relative=max(worst_relative,mbias)
          if(detail_path/='') then
            write(detail_unit,'(a,",",i0,3(",",a),12(",",es26.17e3))')csv_text(trim(path)),rows,trim(left%region),trim(right%region), &
              trim(reason_name(code)),left%n,right%n,left%energy,right%energy,left%p,right%p, &
              (left%energy+left%p)/left%n,(right%energy+right%p)/right%n,low,high,dv,mbias
          end if
        end if
      end if
      rows=rows+1;left=right
    end do
    close(unit)
    if(rows<2) error stop 'at least two EOS rows required'
    write(*,'(a,15(",",i0),2(",",es18.10e3))')csv_text(trim(path)),rows,counts,bad,why,secant_bad,equal_pairs,first_bad,maximum_violation,worst_relative
  end do
  if(detail_path/="") close(detail_unit)
contains
  function csv_text(raw) result(quoted)
    character(*),intent(in) :: raw
    character(:),allocatable :: quoted
    integer :: k
    quoted='"'
    do k=1,len(raw)
      if(raw(k:k)=='"') quoted=quoted//'"'
      quoted=quoted//raw(k:k)
    end do
    quoted=quoted//'"'
  end function

  subroutine read_point(record,value)
    character(*),intent(in) :: record
    type(point),intent(out) :: value
    character(len=128) :: field
    integer :: pos,offset,k,stat,phase
    real(q) :: x(3)
    pos=1
    do k=1,4
      offset=index(record(pos:),',')
      if(offset<=1.or.offset>len(field)+1) error stop 'empty or oversized EOS field'
      field=adjustl(record(pos:pos+offset-2))
      if(len_trim(field)==0) error stop 'empty numeric field'
      if(k<4) then
        if(verify(trim(field),numeric)/=0) error stop 'invalid number token in EOS row'
        read(field,*,iostat=stat)x(k)
      else
        if(verify(trim(field),'0123456789-+')/=0) error stop 'invalid phase integer'
        read(field,*,iostat=stat)phase
      end if
      if(stat/=0) error stop 'invalid EOS row field'
      pos=pos+offset
    end do
    if(len_trim(record(pos:))>len(value%region)) error stop 'oversized EOS region'
    value%region=adjustl(record(pos:))
    if(len_trim(value%region)==0.or.verify(trim(value%region),'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_')/=0) &
      error stop 'invalid EOS region'
    value%n=x(1);value%energy=x(2);value%p=x(3)
  end subroutine

  function reason_name(code) result(name)
    integer,intent(in) :: code
    character(24) :: name
    select case(code)
    case(bad_order); name='order'
    case(bad_ratio); name='acausal_endpoint_ratio'
    case(bad_pressure_area); name='pressure_area'
    case(bad_secant); name='superluminal_secant'
    case default; name='unrecognized_reason'
    end select
  end function

  subroutine check_pair(l,r,result,excess,secant_ok,area_low,area_high)
    type(point),intent(in) :: l,r
    integer,intent(out) :: result
    real(q),intent(out) :: excess
    logical,intent(out) :: secant_ok
    real(q),intent(out),optional :: area_low,area_high
    real(q) :: hl,hr,hll,hlh,hrl,hrh,fmin,fmax,dpmin,dpmax,demax,lo,hi,n_tol
    fmin=1-uq;fmax=1+uq
    hl=(l%energy+l%p)/l%n;hr=(r%energy+r%p)/r%n
    ! Conservative endpoint boxes in n, epsilon and P; ignore correlations.
    hll=hl*fmin/fmax;hlh=hl*fmax/fmin
    hrl=hr*fmin/fmax;hrh=hr*fmax/fmin
    dpmin=r%p*fmin-l%p*fmax;dpmax=r%p*fmax-l%p*fmin
    demax=r%energy*fmax-l%energy*fmin
    secant_ok=r%n*fmax>l%n*fmin.and.demax>0.and.dpmax>=0.and.max(dpmin,0._q)<=demax
    ! dP=n dh and n(h)/h is nondecreasing for 0 <= c_s^2/c^2 <= 1.
    ! Monotonicity of the extremal areas in each input gives outward bounds
    ! over the supplied uncertainty boxes. Equal h allows a density jump.
    lo=l%n*fmin*max(0._q,hrl-hlh)*(hrl+hlh)/(2*hlh)
    hi=r%n*fmax*max(0._q,hrh-hll)*(hrh+hll)/(2*hrh)
    if(present(area_low)) area_low=lo
    if(present(area_high)) area_high=hi
    excess=0;result=bad_order
    if(r%n*fmax<=l%n*fmin.or.demax<=0.or.dpmax<0.or.hrh<hll) return
    if(.not.secant_ok) then
      result=bad_secant
      return
    end if
    if(dpmin<=0.and.hrl<=hlh) then
      ! A near-equality at the stated precision cannot establish coexistence.
      result=equal_support
      return
    end if
    ! The n/h test is necessary even before applying the pressure integral.
    n_tol=128*epsilon(1._q)*max(r%n/hr,l%n/hl)
    if(r%n*fmax/hrl-l%n*fmin/hlh < -n_tol) then
      result=bad_ratio
      return
    end if
    excess=max(0._q,lo-dpmax,dpmin-hi)
    if(excess>128*epsilon(1._q)*max(l%p,r%p,r%n*hr,l%n*hl)) then
      result=bad_pressure_area
    else
      result=compatible
    end if
  end subroutine

  subroutine self_test()
    type(point) :: a,b
    real(q) :: err
    integer :: verdict
    logical :: ok
    a=point(1._q,0.9_q,0.1_q,'core')
    b=point(2._q,1.9_q,0.1_q,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=equal_support) error stop 'plateau'
    b=point(2._q,2.5_q-0.36_q,0.36_q,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=bad_pressure_area.or..not.ok) error stop 'lower area bound missed'
    b=point(2._q,2.5_q-0.57_q,0.57_q,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=bad_pressure_area.or..not.ok) error stop 'upper area bound missed'
    b=point(2._q,4._q-1.6_q,1.6_q,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=compatible.or..not.ok) error stop 'luminal limiting line'
    b=point(4._q,8._q-0.1_q-7._q/3,0.1_q+7._q/3,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=compatible.or..not.ok) error stop 'c_s^2/c^2 = 1/2'
    b=point(1.1_q,2._q,0.2_q,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=bad_ratio.or..not.ok) error stop 'n/h test'
    b=point(1.5_q,1.1_q,0.11_q,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=bad_order.or..not.ok) error stop 'negative slope in h'
    ! A tiny tabulation error changes this near-luminal lower-area verdict.
    b=point(2._q,2.5_q-(0.1_q+0.28125_q-1e-8_q),0.1_q+0.28125_q-1e-8_q,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=bad_pressure_area) error stop 'near-bound nominal test'
    uq=1e-6_q
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=compatible) error stop 'uncertainty envelope does not contain witness'
    uq=roundoff
    b=point(1.1_q,1.0_q,0.4_q,'core')
    call check_pair(a,b,verdict,err,ok)
    if(verdict/=bad_secant.or.ok) error stop 'superluminal pressure secant'
    print '(a)','cold causal endpoint self-tests passed (10)'
  end subroutine
end program
