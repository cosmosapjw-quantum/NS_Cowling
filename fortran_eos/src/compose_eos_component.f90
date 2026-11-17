module compose_eos_component
  use compose_eos_types
  implicit none
  private
  public :: extract_core_rows, extract_crust_rows, verify_core_roundtrip
contains
  subroutine extract_core_rows(source,first_row,last_row,expected_source_n,expected_phase,core,status)
    type(eos_data_t),intent(in)::source
    integer,intent(in)::first_row,last_row,expected_source_n,expected_phase
    type(eos_data_t),intent(out)::core
    character(len=*),intent(out)::status
    integer::nout
    logical::ok
    character(len=256)::message

    core=eos_data_t()
    status='uninitialized'
    if(source%n/=expected_source_n)then
      status='source point count does not match the admitted component manifest';return
    endif
    if(source%nT/=1.or.source%nYq/=1)then
      status='core-only constructor supports exactly one cold T/Yq slice; multidimensional input is fail-closed';return
    endif
    if(abs(source%temperature_mev)>1.e-12_wp)then
      status='core-only constructor requires the source T=0 slice';return
    endif
    if(first_row<1.or.last_row>source%n.or.last_row<first_row)then
      status='admitted core row range is outside the source EOS';return
    endif
    nout=last_row-first_row+1
    if(nout<3)then
      status='standalone core EOS requires at least three admitted density rows';return
    endif
    if(expected_phase>=0)then
      if(.not.source%composition_file_present)then
        status='model-specific homogeneous-core phase guard requires eos.compo';return
      endif
      if(any(source%iphase(first_row:last_row)/=expected_phase))then
        status='source rows do not match the admitted model-specific homogeneous-core Iphase code';return
      endif
    endif
    if(source%nadd>1.and..not.allocated(source%q_extra_full))then
      status='source Nadd>1 but the full Q-extra tensor is unavailable; refusing lossy core export';return
    endif
    if(.not.allocated(source%row_status).or..not.allocated(source%yn_present).or. &
       .not.allocated(source%yp_present).or..not.allocated(source%ye_present).or. &
       .not.allocated(source%ymu_present).or..not.allocated(source%abar_present).or. &
       .not.allocated(source%zbar_present).or..not.allocated(source%yh_present).or. &
       .not.allocated(source%c3_applicable).or..not.allocated(source%c4_applicable))then
      status='source row-level composition masks are incomplete';return
    endif

    ! Strict source-row projection only: no interpolation, extrapolation,
    ! smoothing, re-normalization or closure repair is permitted.
    core=source
    core%n=nout
    core%nb=source%nb(first_row:last_row)
    core%q=source%q(:,first_row:last_row)
    core%q_extra=source%q_extra(first_row:last_row)
    core%iphase=source%iphase(first_row:last_row)
    core%yn=source%yn(first_row:last_row)
    core%yp=source%yp(first_row:last_row)
    core%ye=source%ye(first_row:last_row)
    core%ymu=source%ymu(first_row:last_row)
    core%abar=source%abar(first_row:last_row)
    core%zbar=source%zbar(first_row:last_row)
    core%yh=source%yh(first_row:last_row)
    core%row_status=source%row_status(first_row:last_row)
    core%yn_present=source%yn_present(first_row:last_row)
    core%yp_present=source%yp_present(first_row:last_row)
    core%ye_present=source%ye_present(first_row:last_row)
    core%ymu_present=source%ymu_present(first_row:last_row)
    core%abar_present=source%abar_present(first_row:last_row)
    core%zbar_present=source%zbar_present(first_row:last_row)
    core%yh_present=source%yh_present(first_row:last_row)
    core%c3_applicable=source%c3_applicable(first_row:last_row)
    core%c4_applicable=source%c4_applicable(first_row:last_row)
    ! composition_available is a slice-local summary, not the presence of eos.compo.
    ! The parent unified table may have abundances only in its excluded crust.
    core%composition_available=any(core%yn_present.or.core%yp_present.or.core%ye_present.or. &
      core%ymu_present.or.core%abar_present.or.core%zbar_present.or.core%yh_present)
    ! beta_equilibrated is likewise a property of the selected barotrope,
    ! derived from Q5 on this component slice rather than copied globally.
    core%beta_equilibrated=maxval(abs(core%q(5,:)))<1.e-8_wp
    if(allocated(source%q_full))core%q_full=source%q_full(:,:,:,first_row:last_row)
    if(allocated(source%q_extra_full))core%q_extra_full=source%q_extra_full(:,:,:,first_row:last_row)

    ! Preserve original CompOSE row labels rather than renumbering the core.
    core%inb_start=source%inb_start+first_row-1
    core%iT=source%selected_iT
    core%iYq=source%selected_iYq
    core%selected_iT=source%selected_iT
    core%selected_iYq=source%selected_iYq

    call validate_eos(core,ok,message)
    if(.not.ok)then
      status='extracted core table failed native EOS validation: '//trim(message);return
    endif
    status='ok'
  end subroutine extract_core_rows

  subroutine verify_core_roundtrip(expected,actual,status)
    type(eos_data_t),intent(in)::expected,actual
    character(len=*),intent(out)::status
    status='uninitialized'
    if(expected%n/=actual%n.or.expected%inb_start/=actual%inb_start)then
      status='native CompOSE readback changed core row count or original density-row labels';return
    endif
    if(expected%mn/=actual%mn.or.expected%mp/=actual%mp.or.expected%il/=actual%il.or. &
       expected%temperature_mev/=actual%temperature_mev.or. &
       expected%charge_fraction_grid/=actual%charge_fraction_grid.or. &
       expected%selected_iT/=actual%selected_iT.or.expected%selected_iYq/=actual%selected_iYq.or. &
       expected%t_first/=actual%t_first.or.expected%t_last/=actual%t_last.or. &
       expected%yq_first/=actual%yq_first.or.expected%yq_last/=actual%yq_last.or. &
       expected%nadd/=actual%nadd.or.expected%nadd_max/=actual%nadd_max)then
      status='native CompOSE readback changed core thermodynamic metadata or axes';return
    endif
    if(expected%composition_file_present.neqv.actual%composition_file_present)then
      status='native CompOSE readback changed eos.compo sidecar presence';return
    endif
    if(expected%composition_available.neqv.actual%composition_available)then
      status='native CompOSE readback changed slice-local abundance availability';return
    endif
    if(expected%beta_equilibrated.neqv.actual%beta_equilibrated)then
      status='native CompOSE readback changed Q5-derived beta-equilibrium summary';return
    endif
    if(any(expected%nb/=actual%nb).or.any(expected%q/=actual%q).or. &
       any(expected%q_extra/=actual%q_extra).or.any(expected%iphase/=actual%iphase))then
      status='native CompOSE readback changed nB, Q1-Q7, Q-extra or Iphase';return
    endif
    if(any(expected%yn/=actual%yn).or.any(expected%yp/=actual%yp).or. &
       any(expected%ye/=actual%ye).or.any(expected%ymu/=actual%ymu).or. &
       any(expected%abar/=actual%abar).or.any(expected%zbar/=actual%zbar).or. &
       any(expected%yh/=actual%yh))then
      status='native CompOSE readback changed composition values';return
    endif
    if(any(expected%row_status/=actual%row_status).or. &
       any(expected%yn_present.neqv.actual%yn_present).or. &
       any(expected%yp_present.neqv.actual%yp_present).or. &
       any(expected%ye_present.neqv.actual%ye_present).or. &
       any(expected%ymu_present.neqv.actual%ymu_present).or. &
       any(expected%abar_present.neqv.actual%abar_present).or. &
       any(expected%zbar_present.neqv.actual%zbar_present).or. &
       any(expected%yh_present.neqv.actual%yh_present).or. &
       any(expected%c3_applicable.neqv.actual%c3_applicable).or. &
       any(expected%c4_applicable.neqv.actual%c4_applicable))then
      status='native CompOSE readback changed per-row composition presence or C3/C4 masks';return
    endif
    ! A cold one-point CompOSE axis is serialized by its mandatory grid
    ! file and is already represented by (first,last,value) scalar metadata.
    ! The allocatable t_axis/yq_axis arrays are a convenience in-memory view;
    ! their allocation state is not part of the one-point CompOSE contract.
    ! Multi-point axes, when ever admitted here, still require exact arrays.
    if(expected%nT>1.or.actual%nT>1)then
      if(.not.allocated(expected%t_axis).or..not.allocated(actual%t_axis))then
        status='native CompOSE readback lost a required multi-point T axis';return
      endif
      if(size(expected%t_axis)/=size(actual%t_axis))then
        status='native CompOSE readback changed T axis size';return
      endif
      if(any(expected%t_axis/=actual%t_axis))then
        status='native CompOSE readback changed T axis values';return
      endif
    else if(allocated(expected%t_axis).and.allocated(actual%t_axis))then
      if(size(expected%t_axis)/=size(actual%t_axis))then
        status='native CompOSE readback changed T axis size';return
      endif
      if(any(expected%t_axis/=actual%t_axis))then
        status='native CompOSE readback changed T axis values';return
      endif
    endif
    if(expected%nYq>1.or.actual%nYq>1)then
      if(.not.allocated(expected%yq_axis).or..not.allocated(actual%yq_axis))then
        status='native CompOSE readback lost a required multi-point Yq axis';return
      endif
      if(size(expected%yq_axis)/=size(actual%yq_axis))then
        status='native CompOSE readback changed Yq axis size';return
      endif
      if(any(expected%yq_axis/=actual%yq_axis))then
        status='native CompOSE readback changed Yq axis values';return
      endif
    else if(allocated(expected%yq_axis).and.allocated(actual%yq_axis))then
      if(size(expected%yq_axis)/=size(actual%yq_axis))then
        status='native CompOSE readback changed Yq axis size';return
      endif
      if(any(expected%yq_axis/=actual%yq_axis))then
        status='native CompOSE readback changed Yq axis values';return
      endif
    endif
    ! In a selected cold 1D slice, q(7,n) and q_extra(n) are the
    ! serialized authorities. q_extra_full is an auxiliary in-memory tensor
    ! view allocated whenever Nadd>=1; its allocation state is not itself a
    ! CompOSE data field. For Nadd>1, however, every extra slot must survive.
    if(expected%nYq>1.and.allocated(expected%q_full))then
      if(.not.allocated(actual%q_full))then
        status='native CompOSE readback lost required Q tensor';return
      endif
      if(any(shape(expected%q_full)/=shape(actual%q_full)))then
        status='native CompOSE readback changed Q tensor shape';return
      endif
      if(any(expected%q_full/=actual%q_full))then
        status='native CompOSE readback changed Q tensor values';return
      endif
    endif
    if(expected%nadd>1.and.allocated(expected%q_extra_full))then
      if(.not.allocated(actual%q_extra_full))then
        status='native CompOSE readback lost required Q-extra tensor';return
      endif
      if(any(shape(expected%q_extra_full)/=shape(actual%q_extra_full)))then
        status='native CompOSE readback changed Q-extra tensor shape';return
      endif
      if(any(expected%q_extra_full/=actual%q_extra_full))then
        status='native CompOSE readback changed Q-extra tensor values';return
      endif
    endif
    status='ok'
  end subroutine verify_core_roundtrip

  subroutine extract_crust_rows(source,expected_source_n,outer_last,expected_outer_phase, &
      expected_inner_phase,crust,status)
    type(eos_data_t),intent(in)::source
    integer,intent(in)::expected_source_n,outer_last,expected_outer_phase,expected_inner_phase
    type(eos_data_t),intent(out)::crust
    character(len=*),intent(out)::status

    crust=eos_data_t()
    status='uninitialized'
    if(source%n/=expected_source_n)then
      status='source point count does not match the admitted crust manifest';return
    endif
    if(outer_last<1.or.outer_last>=source%n)then
      status='adjudicated outer-crust row boundary is outside the source grid';return
    endif
    if(.not.source%composition_file_present.or..not.allocated(source%iphase))then
      status='crust phase coverage requires the source eos.compo Iphase field';return
    endif
    if(any(source%iphase(1:outer_last)/=expected_outer_phase))then
      status='source outer-crust rows do not match the adjudicated Iphase code';return
    endif
    if(any(source%iphase(outer_last+1:source%n)/=expected_inner_phase))then
      status='source inner-crust rows do not match the adjudicated Iphase code';return
    endif
    ! The BSk22 physical neutron-drip marker is the 217 -> 218 boundary.
    ! The separate provider splice is nB=2e-3 fm^-3; source grids contain no
    ! sampled rows in the interval between these markers. Never fill that gap.
    if(abs(source%nb(outer_last)-2.69021e-4_wp)>2.0e-9_wp)then
      status='source outer-crust endpoint disagrees with the BSk22 neutron-drip marker';return
    endif
    if(source%nb(outer_last+1)<2.0e-3_wp)then
      status='source has an unadjudicated sampled regime between neutron drip and the 2e-3 provider splice';return
    endif

    ! E390-03B is a source-faithful exact projection. The source thermo values
    ! were generated through the model's admitted CLD/Thomas-Fermi route;
    ! no new microphysical solver or numerical repair runs in this constructor.
    call extract_core_rows(source,1,source%n,expected_source_n,-1,crust,status)
    if(status/='ok')return
    if(crust%n/=source%n.or.crust%inb_start/=source%inb_start)then
      status='crust projection changed source density-grid support or row labels';return
    endif
    status='ok'
  end subroutine extract_crust_rows

end module compose_eos_component
