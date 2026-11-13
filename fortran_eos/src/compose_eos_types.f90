module compose_eos_types
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  integer, parameter, public :: wp=real64
  type, public :: eos_data_t
    integer :: n=0, il=1, iT=1, iYq=1, inb_start=1, nadd=0
    real(wp) :: mn=939._wp, mp=939._wp, temperature_mev=0._wp, charge_fraction_grid=0._wp
    logical :: composition_file_present=.false., composition_available=.false., beta_equilibrated=.false.
    !--- M10 generic CompOSE axes -------------------------------------------
    ! Full T and Yq axis storage: the grid values and their index bounds.
    ! nT and nYq are the axis lengths; t_first/t_last and yq_first/yq_last are
    ! the declared first/last index labels (CompOSE allows 0-based or 1-based).
    ! The cold slice selection is recorded in selected_iT and selected_iYq,
    ! which are the actual iT/iYq labels used in the thermo rows (may differ
    ! from t_first when the source uses a different index convention).
    integer :: nT=0, nYq=0, t_first=0, t_last=-1, yq_first=0, yq_last=-1
    integer :: selected_iT=0, selected_iYq=0
    !--- E390-10B phase 1: full multidim tensor view -------------------------
    ! nadd_max: the declared maximum extra-quantity count across all cells
    !   (0 when the table is Nadd=0; 1 for the current cold-slice cohort).
    ! q_full(nT,nYq,7,n): the full thermo tensor (7 quantities per cell),
    !   opt-in. Unallocated by default: a reader that only needs the cold slice
    !   (the current EOS390 production path) never allocates it, so M3-M6 and
    !   the app commands keep working unchanged. select_slice materializes the
    !   selected (iT,iYq) cell into the 1D q(7,n)/q_extra(n) surface.
    ! q_extra_full(nadd_max,nT,nYq,n): per-cell extra quantities; unallocated by
    !   default and only allocated when the table is genuinely multidim + a
    !   caller requests the full view.
    ! The existing q(7,n), q_extra(n), nadd fields remain and continue to hold
    !   the selected slice (the cold slice by default). See select_slice.
    integer :: nadd_max=0
    real(wp), allocatable :: q_full(:,:,:,:), q_extra_full(:,:,:,:)
    real(wp), allocatable :: t_axis(:), yq_axis(:)
    real(wp), allocatable :: nb(:), q(:,:), q_extra(:)
    integer, allocatable :: iphase(:)
    real(wp), allocatable :: yn(:), yp(:), ye(:), ymu(:), abar(:), zbar(:), yh(:)
    !--- M2A per-row presence masks (zero != missing) -----------------------
    ! row_status: 0=NOT_AVAILABLE (no compo row), 1=PARTIAL (some required field
    !   missing), 2=CHECKED (all required fields present and in physical range),
    !   3=FAIL (present but out of range or C3/C4 violated).
    ! *_present: per-row logical, .true. iff the particle-id/group field was
    !   explicitly written in the eos.compo row (so a written zero is .true.,
    !   a missing entry is .false.). This is the zero-vs-missing distinction.
    ! c3_applicable: per-row, .true. iff the complete baryon inventory needed
    !   for Yn+Yp+Yh=1 is present (Yn, Yp, and the group baryon fraction Yh).
    ! c4_applicable: per-row, .true. iff the complete charge inventory needed
    !   for Yp+(Z/A)Yh=Ye+Ymu is present (Yp, Abar>0 with Zbar, Yh, Ye, Ymu).
    integer, allocatable :: row_status(:)
    logical, allocatable :: yn_present(:), yp_present(:), ye_present(:), ymu_present(:)
    logical, allocatable :: abar_present(:), zbar_present(:), yh_present(:)
    logical, allocatable :: c3_applicable(:), c4_applicable(:)
  end type eos_data_t
  public :: allocate_eos, validate_eos, eos_pressure, eos_energy, eos_mu_b, eos_mu_q, eos_mu_l, eos_free_energy, select_slice
contains
  subroutine allocate_eos(eos,n)
    type(eos_data_t),intent(inout)::eos
    integer,intent(in)::n
    if(allocated(eos%nb)) deallocate(eos%nb)
    if(allocated(eos%q)) deallocate(eos%q)
    if(allocated(eos%q_extra)) deallocate(eos%q_extra)
    ! E390-10B phase 1: the full multidim tensors are opt-in and unallocated by
    ! default. Deallocate them if a caller previously allocated them so
    ! allocate_eos never leaves a stale tensor bound to the new 1D slice.
    if(allocated(eos%q_full)) deallocate(eos%q_full)
    if(allocated(eos%q_extra_full)) deallocate(eos%q_extra_full)
    if(allocated(eos%iphase)) deallocate(eos%iphase)
    if(allocated(eos%yn)) deallocate(eos%yn)
    if(allocated(eos%yp)) deallocate(eos%yp)
    if(allocated(eos%ye)) deallocate(eos%ye)
    if(allocated(eos%ymu)) deallocate(eos%ymu)
    if(allocated(eos%abar)) deallocate(eos%abar)
    if(allocated(eos%zbar)) deallocate(eos%zbar)
    if(allocated(eos%yh)) deallocate(eos%yh)
    if(allocated(eos%row_status)) deallocate(eos%row_status)
    if(allocated(eos%yn_present)) deallocate(eos%yn_present)
    if(allocated(eos%yp_present)) deallocate(eos%yp_present)
    if(allocated(eos%ye_present)) deallocate(eos%ye_present)
    if(allocated(eos%ymu_present)) deallocate(eos%ymu_present)
    if(allocated(eos%abar_present)) deallocate(eos%abar_present)
    if(allocated(eos%zbar_present)) deallocate(eos%zbar_present)
    if(allocated(eos%yh_present)) deallocate(eos%yh_present)
    if(allocated(eos%c3_applicable)) deallocate(eos%c3_applicable)
    if(allocated(eos%c4_applicable)) deallocate(eos%c4_applicable)
    ! M10 axes (t_axis/yq_axis) are intentionally NOT deallocated here: they
    ! are populated by read_compose BEFORE allocate_eos is called, and must
    ! survive so the cold-slice selection is recorded on the resulting eos.
    eos%n=n
    allocate(eos%nb(n),eos%q(7,n),eos%q_extra(n),eos%iphase(n))
    allocate(eos%yn(n),eos%yp(n),eos%ye(n),eos%ymu(n),eos%abar(n),eos%zbar(n),eos%yh(n))
    allocate(eos%row_status(n))
    allocate(eos%yn_present(n),eos%yp_present(n),eos%ye_present(n),eos%ymu_present(n))
    allocate(eos%abar_present(n),eos%zbar_present(n),eos%yh_present(n))
    allocate(eos%c3_applicable(n),eos%c4_applicable(n))
    eos%nb=0._wp; eos%q=0._wp; eos%q_extra=0._wp; eos%iphase=0
    eos%yn=0._wp; eos%yp=0._wp; eos%ye=0._wp; eos%ymu=0._wp
    eos%abar=0._wp; eos%zbar=0._wp; eos%yh=0._wp
    ! M2A defaults: no compo row yet (NOT_AVAILABLE), nothing present, C3/C4 not applicable.
    eos%row_status=0
    eos%yn_present=.false.; eos%yp_present=.false.; eos%ye_present=.false.; eos%ymu_present=.false.
    eos%abar_present=.false.; eos%zbar_present=.false.; eos%yh_present=.false.
    eos%c3_applicable=.false.; eos%c4_applicable=.false.
  end subroutine allocate_eos

  subroutine validate_eos(eos,ok,message)
    type(eos_data_t),intent(in)::eos
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    integer::i
    ok=.false.; message='EOS arrays are not allocated'
    if(eos%n<3) then
      message='EOS requires at least three density points';return
    endif
    if(.not.allocated(eos%nb).or..not.allocated(eos%q)) return
    if(size(eos%nb)/=eos%n.or.any(shape(eos%q)/=[7,eos%n])) then
      message='EOS array shapes are inconsistent';return
    endif
    if(.not.all(ieee_finite(eos%nb)).or..not.all(ieee_finite(eos%q))) then
      message='EOS primary arrays contain nonfinite values';return
    endif
    if(any(eos%nb<=0._wp).or.any(eos%nb(2:)<=eos%nb(:eos%n-1))) then
      message='EOS density grid must be positive and strictly increasing';return
    endif
    if(eos%mn<=0._wp.or.eos%mp<=0._wp.or..not.ieee_finite(eos%mn+eos%mp)) then
      message='invalid nucleon masses';return
    endif
    do i=1,eos%n
      if(eos%q(1,i)*eos%nb(i)<=0._wp) then
        message='nonpositive pressure encountered';return
      endif
    enddo
    ok=.true.;message='ok'
  end subroutine validate_eos

  pure function eos_pressure(eos) result(p)
    type(eos_data_t),intent(in)::eos
    real(wp)::p(eos%n)
    p=eos%q(1,:)*eos%nb
  end function
  pure function eos_mu_b(eos) result(mu)
    type(eos_data_t),intent(in)::eos
    real(wp)::mu(eos%n)
    mu=(eos%q(3,:)+1._wp)*eos%mn
  end function
  pure function eos_mu_q(eos) result(mu)
    type(eos_data_t),intent(in)::eos
    real(wp)::mu(eos%n)
    mu=eos%q(4,:)*eos%mn
  end function
  pure function eos_mu_l(eos) result(mu)
    type(eos_data_t),intent(in)::eos
    real(wp)::mu(eos%n)
    mu=eos%q(5,:)*eos%mn
  end function
  pure function eos_free_energy(eos) result(f)
    type(eos_data_t),intent(in)::eos
    real(wp)::f(eos%n)
    f=(eos%q(6,:)+1._wp)*eos%mn
  end function
  pure function eos_energy(eos) result(e)
    type(eos_data_t),intent(in)::eos
    real(wp)::e(eos%n)
    e=eos%nb*(eos%q(7,:)+1._wp)*eos%mn
  end function
  pure elemental logical function ieee_finite(x)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(wp),intent(in)::x
    ieee_finite=ieee_is_finite(x)
  end function

  !--- E390-10B phase 1: select_slice --------------------------------------
  ! Materialize the (iT,iYq) cell of the full multidim thermo tensor into the
  ! 1D q(7,n)/q_extra(n) surface that M3-M6 and the app commands consume. This
  ! bridges the tensor to the cold-slice surface API without changing their
  ! call sites. The 1D fields q/q_extra are the contract for downstream code;
  ! select_slice is the only sanctioned way to populate them from q_full.
  !
  ! Convention: iT/iYq are 1-based cell indices into q_full(nT,nYq,n). The
  ! caller is responsible for converting any 0-based CompOSE axis labels to
  ! 1-based before calling. selected_iT/selected_iYq are updated to record the
  ! chosen cell labels for provenance.
  !
  ! status: 'ok' on success; otherwise a human-readable reason:
  !   'q_full not allocated' (cold-slice reader or opt-out caller)
  !   'iT/iYq out of range [1,nT]/[1,nYq]'
  !   'q_full cell shape mismatch'
  subroutine select_slice(eos,iT,iYq,status)
    type(eos_data_t),intent(inout)::eos
    integer,intent(in)::iT,iYq
    character(len=*),intent(out)::status
    integer::n
    status='uninitialized'
    if(.not.allocated(eos%q_full)) then
      status='q_full not allocated';return
    endif
    n=eos%n
    if(n<1) then;status='eos%n must be >= 1';return;endif
    if(iT<1.or.iT>eos%nT.or.iYq<1.or.iYq>eos%nYq) then
      write(status,'(a,2(a,i0),a)')'iT/iYq out of range',' (iT=',iT,', nT=',eos%nT,')'
      return
    endif
    ! Composition/Iphase/presence masks currently belong only to the first
    ! materialized (T,Yq) cell. Do not copy a different thermo tensor cell into
    ! the legacy 1D view while leaving stale composition attached to it.
    if((eos%composition_file_present.or.eos%composition_available).and.(iT/=1.or.iYq/=1)) then
      status='composition is materialized only for slice (iT=1, iYq=1); full multidimensional composition support is required'
      return
    endif
    if(size(eos%q_full,1)/=eos%nT.or.size(eos%q_full,2)/=eos%nYq.or.size(eos%q_full,3)/=7 &
       .or.size(eos%q_full,4)/=n) then
      status='q_full cell shape mismatch';return
    endif
    if(.not.allocated(eos%q).or.size(eos%q,2)/=n) then
      status='q surface not allocated to eos%n';return
    endif
    eos%q(:,1:n)=eos%q_full(iT,iYq,:,1:n)
    ! extra quantities: only copy if the tensor carries them and the surface
    ! exists. nadd records how many extra columns this cell uses (<= nadd_max).
    if(allocated(eos%q_extra_full)) then
      if(size(eos%q_extra_full,1)/=eos%nadd_max.or.size(eos%q_extra_full,2)/=eos%nT &
         .or.size(eos%q_extra_full,3)/=eos%nYq.or.size(eos%q_extra_full,4)/=n) then
        status='q_extra_full cell shape mismatch';return
      endif
      if(allocated(eos%q_extra).and.size(eos%q_extra)/=n) then
        status='q_extra surface not allocated to eos%n';return
      endif
      eos%q_extra(1:n)=0._wp
      if(eos%nadd_max>=1.and.eos%nadd>=1) then
        eos%q_extra(1:n)=eos%q_extra_full(1,iT,iYq,1:n)
      endif
    else
      if(allocated(eos%q_extra)) eos%q_extra(1:n)=0._wp
    endif
    eos%selected_iT=iT;eos%selected_iYq=iYq
    status='ok'
  end subroutine select_slice
end module compose_eos_types
