! C17 extended original-source Maxwell branch control and onset diagnostic.
!
! Local continuation is ADMITTED only for exactly one physical root with
! (i) the same W/V nodal labels, (ii) a bounded local frequency step,
! and (iii) separately acceptable high/low phase W-shape overlaps.
! Overlap is a descriptive areal-radius W profile correlation, not a
! conserved inner product. It must never be used to infer global branch
! identity or absence of an avoided crossing.
!
! The standard small-core Maxwell density jump test from
! Lindblom, Phys Rev D 58 024008 (1998), gr-qc/9802072:
! delta=(eps_plus-eps_minus)/(eps_minus+P_t)
! delta_crit=(eps_minus+3 P_t)/(2(eps_minus+P_t)).
! Applicable as a transition-onset diagnostic, NOT as a full Cowling
! mode stability certificate. Inputs P and epsilon share physical units.
module compose_eos_cowling_extended_branch
  use iso_c_binding, only: c_int,c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: cowling_native_select_mode_branch_overlap
  public :: cowling_native_maxwell_onset_threshold
contains
  subroutine cowling_native_select_mode_branch_overlap(n,frequency,nw,nv, &
       roots,high_overlap,low_overlap,target_nodes,previous_hz, &
       max_fraction,min_overlap,index,status) bind(c)
    integer(c_int),value :: n,target_nodes
    real(c_double),intent(in)::frequency(n),high_overlap(n),low_overlap(n)
    integer(c_int),intent(in)::nw(n),nv(n),roots(n)
    real(c_double),value :: previous_hz,max_fraction,min_overlap
    integer(c_int),intent(out)::index,status
    integer::i,neig,npass,last
    index=-1;status=1
    if(n<1.or.n>128.or.target_nodes<0.or.target_nodes>16)return
    if(.not.all(ieee_is_finite([previous_hz,max_fraction,min_overlap])))return
    if(previous_hz<=0._c_double.or.max_fraction<=0._c_double &
          .or.max_fraction>=1._c_double)return
    if(min_overlap<=0._c_double.or.min_overlap>1._c_double)return
    if(.not.all(ieee_is_finite(frequency)))return
    if(any(frequency<=0._c_double))return
    neig=0;npass=0;last=-1
    do i=1,n
      if(roots(i)/=0_c_int)cycle
      if(nw(i)/=target_nodes.or.nv(i)/=target_nodes)cycle
      if(abs(frequency(i)-previous_hz)>max_fraction*previous_hz)cycle
      neig=neig+1
      if(.not.ieee_is_finite(high_overlap(i)).or. &
         .not.ieee_is_finite(low_overlap(i)))return
      if(high_overlap(i)<min_overlap.or.low_overlap(i)<min_overlap)cycle
      if(high_overlap(i)>1._c_double.or.low_overlap(i)>1._c_double)return
      npass=npass+1;last=i-1
    enddo
    if(npass==1)then
      index=last;status=0
    else if(npass>1)then
      status=3  ! legitimate ambiguity: never silently choose closest
    else if(neig==0)then
      status=4  ! no root with correct node identity and freq step
    else
      status=5  ! frequency/node candidate exists, physical overlap fails
    endif
  end subroutine

  subroutine cowling_native_maxwell_onset_threshold( &
      eps_minus,eps_plus,p_transition,delta,delta_critical,margin,above,status) bind(c)
    real(c_double),value::eps_minus,eps_plus,p_transition
    real(c_double),intent(out)::delta,delta_critical,margin
    integer(c_int),intent(out)::above,status
    real(c_double)::den
    delta=0._c_double;delta_critical=0._c_double
    margin=0._c_double;above=-1;status=1
    if(.not.all(ieee_is_finite([eps_minus,eps_plus,p_transition])))return
    if(eps_minus<=0._c_double.or.eps_plus<=eps_minus.or.p_transition<0._c_double)return
    den=eps_minus+p_transition
    if(den<=0._c_double)return
    delta=(eps_plus-eps_minus)/den
    delta_critical=(eps_minus+3._c_double*p_transition)/(2._c_double*den)
    margin=delta-delta_critical
    if(.not.all(ieee_is_finite([delta,delta_critical,margin])))then
      delta=0._c_double;delta_critical=0._c_double
      margin=0._c_double;return
    endif
    if(margin>0._c_double)then
      above=1
    else
      above=0
    endif
    status=0
  end subroutine
end module compose_eos_cowling_extended_branch
