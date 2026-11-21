! C16 phase-local physical eigenfunction nodal classifier.
! Counts actual strict sign changes after ignoring amplitudes below a
! relative component/phase peak cutoff. Does not connect nodes across
! the sharp Maxwell density discontinuity.
module compose_eos_cowling_mode_nodes
  use iso_c_binding,only:c_int,c_double
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::cowling_native_phase_nodes,cowling_native_select_mode_branch
contains
  pure integer function count_nodes(x,tolerance) result(n)
    real(c_double),intent(in)::x(:),tolerance
    real(c_double)::threshold,peak
    integer::i,last
    n=0
    peak=maxval(abs(x))
    if(peak==0._c_double)return
    threshold=tolerance*peak
    last=0
    do i=1,size(x)
      if(abs(x(i))<=threshold)cycle
      if(last/=0)then
        if(sign(1._c_double,x(i))/=sign(1._c_double,x(last)))n=n+1
      endif
      last=i
    enddo
  end function

  subroutine cowling_native_phase_nodes(n,phase,w,v,tolerance,nw,nv,status) bind(c)
    integer(c_int),value::n
    integer(c_int),intent(in)::phase(n)
    real(c_double),intent(in)::w(n),v(n)
    real(c_double),value::tolerance
    integer(c_int),intent(out)::nw,nv,status
    integer::j,first,last
    nw=0;nv=0;status=1
    if(n<2)return
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0._c_double &
        .or.tolerance>=1._c_double)return
    if(.not.all(ieee_is_finite(w)).or..not.all(ieee_is_finite(v)))return
    if(phase(1)<0)return
    if(any(phase(2:n)<phase(1:n-1)))return
    if(any(phase(2:n)>phase(1:n-1)+1))return
    first=1
    do while(first<=n)
      last=first
      do while(last<n)
        if(phase(last+1)/=phase(first))exit
        last=last+1
      enddo
      nw=nw+count_nodes(w(first:last),tolerance)
      nv=nv+count_nodes(v(first:last),tolerance)
      first=last+1
    enddo
    status=0
  end subroutine
  ! A unique, physically accepted local branch successor with the
  ! same pair of radial node counts and a bounded frequency increment.
  ! We deliberately do NOT choose the closest of two eligible roots:
  ! ambiguous mode identity / an avoided-crossing case remains blocked.
  ! Selected index uses Python/C zero-based convention.
  subroutine cowling_native_select_mode_branch(n,freq,nw,nv,accept, &
     target_nodes,previous_hz,max_fraction,index,status) bind(c)
    integer(c_int),value::n,target_nodes
    real(c_double),intent(in)::freq(n)
    integer(c_int),intent(in)::nw(n),nv(n),accept(n)
    real(c_double),value::previous_hz,max_fraction
    integer(c_int),intent(out)::index,status
    integer::j,ngood,which
    index=-1;status=1
    if(n<1.or.n>128.or.target_nodes<0.or.target_nodes>16)return
    if(.not.ieee_is_finite(previous_hz).or.previous_hz<=0._c_double)return
    if(.not.ieee_is_finite(max_fraction).or.max_fraction<=0._c_double &
       .or.max_fraction>=1._c_double)return
    if(.not.all(ieee_is_finite(freq)))return
    if(any(freq<=0._c_double))return
    ngood=0;which=-1
    do j=1,n
      if(accept(j)/=0_c_int)cycle
      if(nw(j)/=target_nodes.or.nv(j)/=target_nodes)cycle
      if(abs(freq(j)-previous_hz)>max_fraction*previous_hz)cycle
      ngood=ngood+1
      which=j-1
    enddo
    if(ngood==0)then
      status=4
    else if(ngood>1)then
      status=3
    else
      index=which;status=0
    endif
  end subroutine
end module compose_eos_cowling_mode_nodes
