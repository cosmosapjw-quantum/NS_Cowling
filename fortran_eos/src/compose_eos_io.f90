module compose_eos_io
  use compose_eos_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: read_compose,write_compose,read_internal_csv,write_internal_csv,write_normalized_csv,write_metadata,read_metadata
contains
  subroutine split_tokens(line,tokens,nt)
    character(len=*),intent(in)::line
    character(len=96),intent(out)::tokens(:)
    integer,intent(out)::nt
    integer::i,j,n
    tokens=''; nt=0; n=len_trim(line); i=1
    do while(i<=n)
      if(line(i:i)==' '.or.line(i:i)==','.or.iachar(line(i:i))==9) then
        i=i+1; cycle
      endif
      j=i
      do while(j<=n)
        if(line(j:j)==' '.or.line(j:j)==','.or.iachar(line(j:j))==9) exit
        j=j+1
      enddo
      nt=nt+1
      if(nt>size(tokens)) return
      tokens(nt)=line(i:j-1)
      i=j+1
    enddo
  end subroutine split_tokens

  subroutine next_line(unit,line,ios)
    integer,intent(in)::unit
    character(len=*),intent(out)::line
    integer,intent(out)::ios
    do
      read(unit,'(A)',iostat=ios) line
      if(ios/=0) return
      if(len_trim(adjustl(line))==0) cycle
      line=adjustl(line)
      if(line(1:1)=='#') cycle
      return
    enddo
  end subroutine next_line

  subroutine token_int(tokens,k,value,ok)
    character(len=*),intent(in)::tokens(:)
    integer,intent(in)::k
    integer,intent(out)::value
    logical,intent(out)::ok
    integer::ios
    value=0; ok=.false.
    if(k<1.or.k>size(tokens)) return
    read(tokens(k),*,iostat=ios)value
    ok=ios==0
  end subroutine token_int

  subroutine token_real(tokens,k,value,ok)
    character(len=*),intent(in)::tokens(:)
    integer,intent(in)::k
    real(wp),intent(out)::value
    logical,intent(out)::ok
    integer::ios
    value=0._wp; ok=.false.
    if(k<1.or.k>size(tokens)) return
    read(tokens(k),*,iostat=ios)value
    ok=ios==0
    if(ok) ok=ieee_is_finite(value)
  end subroutine token_real

  subroutine read_grid(path,i0,i1,values,status)
    character(len=*),intent(in)::path
    integer,intent(out)::i0,i1
    real(wp),allocatable,intent(out)::values(:)
    character(len=*),intent(out)::status
    integer::u,ios,k
    character(len=4096)::line
    i0=0;i1=-1
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='cannot open '//trim(path);return;endif
    call next_line(u,line,ios)
    if(ios/=0) then;close(u);status='missing first grid index';return;endif
    read(line,*,iostat=ios)i0
    if(ios/=0) then;close(u);status='invalid first grid index';return;endif
    call next_line(u,line,ios)
    if(ios/=0) then;close(u);status='missing last grid index';return;endif
    read(line,*,iostat=ios)i1
    if(ios/=0.or.i1<i0) then;close(u);status='invalid last grid index';return;endif
    allocate(values(i1-i0+1))
    do k=1,size(values)
      call next_line(u,line,ios)
      if(ios/=0) then;close(u);status='grid has fewer values than declared';return;endif
      read(line,*,iostat=ios)values(k)
      if(ios/=0.or..not.ieee_is_finite(values(k))) then
        close(u);status='grid contains malformed or nonfinite value';return
      endif
    enddo
    close(u);status='ok'
  end subroutine read_grid

  subroutine read_compose(directory,eos,status)
    character(len=*),intent(in)::directory
    type(eos_data_t),intent(out)::eos
    character(len=*),intent(out)::status
    integer::nb0,nb1,t0,t1,y0,y1,u,ios,i,k,nt,it,ib,iy,nadd,extraid,npairs,nquad,pos,pid,qid,phase
    integer::cold_t_idx,cold_yq_idx,declared_t_idx,declared_yq_idx,selected_t,selected_yq
    integer::iy_local,iy_min,iy_max
    logical::t_idx_resolved,yq_idx_resolved
    character(len=4096)::line
    character(len=96)::tok(256)
    real(wp),allocatable::nb(:),tv(:),yv(:)
    real(wp)::qv(7),v,a,z,yi,tval,residual
    logical::ok,file_exists
    logical,allocatable::seen(:)
    logical,allocatable::seen_2d(:,:)
    character(len=1024)::path
    eos=eos_data_t();status='uninitialized'
    path=trim(directory)//'/eos.nb'
    call read_grid(trim(path),nb0,nb1,nb,status);if(status/='ok')return
    call read_grid(trim(directory)//'/eos.t',t0,t1,tv,status);if(status/='ok')return
    call read_grid(trim(directory)//'/eos.yq',y0,y1,yv,status);if(status/='ok')return
    if(size(nb)<3.or.any(nb<=0._wp).or.any(nb(2:)<=nb(:size(nb)-1))) then
      status='eos.nb must contain at least three positive strictly increasing values';return
    endif
    ! M10: store the full T and Yq axes and their declared index bounds on the
    ! eos object, so downstream code can see the generic 1D/2D/3D layout. The
    ! existing reader remains scoped to a single cold T=0 slice (no implicit
    ! flattening, no axis averaging, no nearest-index substitution), but it
    ! now resolves the cold slice by VALUE rather than by literal index label.
    eos%nT=size(tv); eos%nYq=size(yv)
    eos%t_first=t0; eos%t_last=t1; eos%yq_first=y0; eos%yq_last=y1
    allocate(eos%t_axis(size(tv))); eos%t_axis=tv
    allocate(eos%yq_axis(size(yv))); eos%yq_axis=yv
    ! E390-10B phase 3: relax the nYq==1 bound. A cold slice now allows a
    ! multi-point Yq axis (nYq>1) while keeping nT==1 and T=0. Every Yq thermo
    ! row is stored into q_full(1,iy_local,7,k); the 1D q(7,n) surface still
    ! materializes the selected (first) Yq for backward compatibility with
    ! M3-M6 downstream code. The nT==1 + T=0 cold-slice requirement stays.
    if(size(tv)/=1.or.abs(tv(1))>1.e-12_wp) then
      status='only a single cold T=0 eos.t slice is supported';return
    endif
    ! M10: resolve the cold-slice T index by VALUE. CompOSE allows 0-based or
    ! 1-based T indexing: some sources declare eos.t as index 0..0 with value
    ! T=0 (e.g. RG_SLY4, RG_KDE0v1) and write thermo rows with iT=1, while
    ! others declare eos.t as index 1..1 (e.g. GPPVA_FSU2H). Both refer to the
    ! same single cold-T=0 point. Accept either the declared first index or
    ! the declared first index + 1 if the grid value is T=0; record the actual
    ! label used in the thermo rows as selected_iT. This is NOT guessing the
    ! temperature - it is recognizing the index convention for the cold slice.
    declared_t_idx=t0
    t_idx_resolved=.false.
    if(t0>=0.and.t0<=size(tv)) then
      tval=tv(max(1,t0-t0+1))  ! value at the declared first index
      if(abs(tval)<=1.e-12_wp) then
        cold_t_idx=t0; t_idx_resolved=.true.
      endif
    endif
    if(.not.t_idx_resolved) then
      ! fall back: try t0+1 (1-based convention for a 0-based declaration)
      if(t0+1>=1.and.t0+1<=size(tv)+t0) then
        tval=tv(1)
        if(abs(tval)<=1.e-12_wp) then
          cold_t_idx=t0+1; t_idx_resolved=.true.
        endif
      endif
    endif
    if(.not.t_idx_resolved) then
      status='cold T=0 slice could not be resolved from eos.t grid';return
    endif
    declared_yq_idx=y0
    yq_idx_resolved=.false.
    if(y0>=0.and.y0<=size(yv)) then
      cold_yq_idx=y0; yq_idx_resolved=.true.
    endif
    if(.not.yq_idx_resolved) then
      status='cold Yq slice could not be resolved from eos.yq grid';return
    endif
    call allocate_eos(eos,size(nb));eos%nb=nb;eos%inb_start=nb0;eos%temperature_mev=tv(1)
    eos%charge_fraction_grid=yv(1)
    eos%selected_iT=cold_t_idx; eos%selected_iYq=cold_yq_idx
    ! E390-10B phase 3: allocate the multidim q_full(nT,nYq,7,n) tensor when
    ! the table carries more than one Yq point (cold slice: nT==1). Slot
    ! (1,selected_iYq_local,:,k) mirrors the 1D q(7,n) surface so M3-M6
    ! downstream code that reads the 1D fields is unaffected.
    if(eos%nYq>1) then
      if(.not.allocated(eos%q_full)) then
        allocate(eos%q_full(eos%nT,eos%nYq,7,eos%n))
        eos%q_full=0._wp
      endif
    endif
    path=trim(directory)//'/eos.thermo'
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='missing or unreadable eos.thermo';return;endif
    call next_line(u,line,ios)
    if(ios/=0) then;close(u);status='empty eos.thermo';return;endif
    read(line,*,iostat=ios)eos%mn,eos%mp,eos%il
    if(ios/=0.or.eos%mn<=0._wp.or.eos%mp<=0._wp.or.(eos%il/=0.and.eos%il/=1)) then
      close(u);status='invalid eos.thermo header';return
    endif
    ! E390-10B phase 3: track coverage over the full (Yq,density) grid when
    ! nYq>1, or the 1D density grid when nYq==1 (legacy cold slice). The 1D
    ! q(7,n) surface always materializes the selected (first) Yq; the other
    ! Yq rows go into q_full(1,iy_local,7,k).
    if(eos%nYq>1) then
      allocate(seen_2d(eos%nYq,eos%n));seen_2d=.false.
    else
      allocate(seen(eos%n));seen=.false.
    endif
    eos%nadd=-1
    ! accepted Yq label range: CompOSE allows 0-based or 1-based Yq indexing
    ! (the cold-slice value). Accept the declared first index and the 1-based
    ! variant; record the actual label as selected_iYq. For nYq>1 the rows
    ! cover consecutive Yq labels from y0 to y0+eos%nYq-1 (1-based convention).
    iy_min=eos%selected_iYq
    iy_max=eos%selected_iYq+eos%nYq-1
    do
      call next_line(u,line,ios)
      if(ios<0) exit
      if(ios>0) then;close(u);status='error reading eos.thermo';return;endif
      call split_tokens(line,tok,nt)
      if(nt<11) then;close(u);status='short eos.thermo row';return;endif
      call token_int(tok,1,it,ok);if(.not.ok) then;close(u);status='invalid iT';return;endif
      call token_int(tok,2,ib,ok);if(.not.ok) then;close(u);status='invalid inb';return;endif
      call token_int(tok,3,iy,ok);if(.not.ok) then;close(u);status='invalid iYq';return;endif
      ! M10/phase 3: T must still be the single cold slice (T=0); the T label
      ! accepts either the declared first index or the 1-based variant. Yq
      ! may now span nYq consecutive labels (phase 3 relaxation).
      if((it/=eos%selected_iT .and. it/=eos%selected_iT+1)) then
        close(u);status='thermo table is not a single selected T slice';return
      endif
      if(iy<iy_min.or.iy>iy_max) then
        close(u);status='thermo Yq label outside the declared eos.yq grid';return
      endif
      iy_local=iy-iy_min+1
      eos%selected_iT=it
      k=ib-nb0+1
      if(k<1.or.k>eos%n) then;close(u);status='thermo density index outside eos.nb';return;endif
      if(eos%nYq>1) then
        if(seen_2d(iy_local,k)) then;close(u);status='duplicate eos.thermo row key';return;endif
      else
        if(seen(k)) then;close(u);status='duplicate eos.thermo row key';return;endif
      endif
      do i=1,7
        call token_real(tok,3+i,qv(i),ok)
        if(.not.ok) then;close(u);status='invalid/nonfinite Q column';return;endif
      enddo
      call token_int(tok,11,nadd,ok)
      if(.not.ok.or.nadd<0) then;close(u);status='Nadd must be a non-negative integer';return;endif
      if(eos%nadd<0) then
        eos%nadd=nadd
        ! E390-10B phase 2/3: track the declared max extra count and allocate
        ! the full q_extra_full(nadd_max,nT,nYq,n) tensor once we know the
        ! per-row Nadd. q_extra_full holds every extra column while the 1D
        ! q_extra(n) surface keeps the first one for backward compatibility.
        eos%nadd_max=nadd
        if(nadd>=1.and.eos%nT>=1.and.eos%nYq>=1.and.allocated(eos%nb)) then
          if(.not.allocated(eos%q_extra_full)) then
            allocate(eos%q_extra_full(nadd,eos%nT,eos%nYq,eos%n))
            eos%q_extra_full=0._wp
          endif
        endif
      else if(eos%nadd/=nadd) then
        close(u);status='Nadd changes across eos.thermo rows';return
      endif
      if(nt/=11+nadd) then;close(u);status='Nadd does not match thermo row width';return;endif
      ! materialize the selected (first) Yq into the 1D q(7,n) surface
      if(iy_local==1) eos%q(:,k)=qv
      ! store every Yq row into the multidim tensor (cold slice: nT=1)
      if(allocated(eos%q_full)) eos%q_full(1,iy_local,:,k)=qv
      if(nadd>=1) then
        do i=1,nadd
          call token_real(tok,11+i,v,ok)
          if(.not.ok) then;close(u);status='invalid q_extra column';return;endif
          ! store every column in the tensor (slot 1 mirrors 1D q_extra for
          ! the selected Yq only)
          if(allocated(eos%q_extra_full)) eos%q_extra_full(i,1,iy_local,k)=v
          if(i==1.and.iy_local==1) eos%q_extra(k)=v
        enddo
      endif
      if(eos%nYq>1) then
        seen_2d(iy_local,k)=.true.
      else
        seen(k)=.true.
      endif
    enddo
    close(u)
    if(eos%nYq>1) then
      if(.not.all(seen_2d)) then;status='eos.thermo does not cover every (Yq,density) grid point';return;endif
    else
      if(.not.all(seen)) then;status='eos.thermo does not cover every density grid point';return;endif
    endif
    if(eos%nadd<0)eos%nadd=0
    path=trim(directory)//'/eos.compo'
    if(allocated(seen)) deallocate(seen)
    if(allocated(seen_2d)) deallocate(seen_2d)
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    file_exists=ios==0
    eos%composition_file_present=file_exists
    if(eos%composition_file_present) then
      allocate(seen(eos%n));seen=.false.;eos%composition_available=.false.
      ! M2A: per-row presence masks. row_status starts NOT_AVAILABLE (0).
      ! The parser fills presence flags as each field is read; range/C3/C4
      ! checks then promote to CHECKED (2) or FAIL (3), or leave PARTIAL (1).
      eos%row_status=0
      eos%yn_present=.false.; eos%yp_present=.false.; eos%ye_present=.false.; eos%ymu_present=.false.
      eos%abar_present=.false.; eos%zbar_present=.false.; eos%yh_present=.false.
      eos%c3_applicable=.false.; eos%c4_applicable=.false.
      do
        call next_line(u,line,ios)
        if(ios<0) exit
        if(ios>0) then;close(u);status='error reading eos.compo';return;endif
        call split_tokens(line,tok,nt)
        if(nt<6) then;close(u);status='short eos.compo row';return;endif
        call token_int(tok,1,it,ok);if(.not.ok) then;close(u);status='invalid compo iT';return;endif
        call token_int(tok,2,ib,ok);if(.not.ok) then;close(u);status='invalid compo inb';return;endif
        call token_int(tok,3,iy,ok);if(.not.ok) then;close(u);status='invalid compo iYq';return;endif
        k=ib-nb0+1
        if(k<1.or.k>eos%n.or. &
           ((it/=eos%selected_iT .and. it/=eos%selected_iT+1) .or. iy/=eos%selected_iYq)) then
          close(u);status='composition key outside selected thermo slice';return
        endif
        eos%selected_iT=it
        if(seen(k)) then;close(u);status='duplicate eos.compo row key';return;endif
        call token_int(tok,4,phase,ok)
        if(.not.ok) then;close(u);status='invalid compo phase';return;endif
        eos%iphase(k)=phase
        call token_int(tok,5,npairs,ok)
        if(.not.ok.or.npairs<0.or.npairs>16) then;close(u);status='invalid Npairs';return;endif
        pos=6
        do i=1,npairs
          if(pos+1>nt) then;close(u);status='truncated particle pair';return;endif
          call token_int(tok,pos,pid,ok);if(.not.ok) then;close(u);status='invalid particle id';return;endif
          call token_real(tok,pos+1,v,ok);if(.not.ok) then;close(u);status='invalid particle fraction';return;endif
          select case(pid)
          case(10)
            if(eos%yn_present(k))then;close(u);status='duplicate CompOSE neutron id in eos.compo row';return;endif
            eos%yn(k)=v;eos%yn_present(k)=.true.
          case(11)
            if(eos%yp_present(k))then;close(u);status='duplicate CompOSE proton id in eos.compo row';return;endif
            eos%yp(k)=v;eos%yp_present(k)=.true.
          case(0)
            if(eos%ye_present(k))then;close(u);status='duplicate CompOSE electron id in eos.compo row';return;endif
            eos%ye(k)=v;eos%ye_present(k)=.true.
          case(1)
            if(eos%ymu_present(k))then;close(u);status='duplicate CompOSE muon id in eos.compo row';return;endif
            eos%ymu(k)=v;eos%ymu_present(k)=.true.
          case default
            close(u);status='unsupported CompOSE particle id in this EOSData model';return
          end select
          ! Presence, not numerical magnitude, defines availability: an explicit
          ! zero is a measured/written value and must survive round-trip.
          eos%composition_available=.true.
          pos=pos+2
        enddo
        if(pos>nt) then;close(u);status='missing Nquad';return;endif
        call token_int(tok,pos,nquad,ok)
        if(.not.ok.or.nquad<0.or.nquad>1) then;close(u);status='only Nquad=0 or 1 is supported';return;endif
        pos=pos+1
        do i=1,nquad
          if(pos+3>nt) then;close(u);status='truncated nuclear-group quadruple';return;endif
          call token_int(tok,pos,qid,ok);if(.not.ok) then;close(u);status='invalid nuclear group id';return;endif
          call token_real(tok,pos+1,a,ok);if(.not.ok) then;close(u);status='invalid Abar';return;endif
          call token_real(tok,pos+2,z,ok);if(.not.ok) then;close(u);status='invalid Zbar';return;endif
          call token_real(tok,pos+3,yi,ok);if(.not.ok) then;close(u);status='invalid YI';return;endif
          if(qid/=1.or.a<0._wp.or.z<0._wp.or.z>a.or.yi<0._wp) then
            close(u);status='unsupported or invalid nuclear group';return
          endif
          if(eos%abar_present(k))then;close(u);status='duplicate CompOSE nuclear group in eos.compo row';return;endif
          eos%abar(k)=a;eos%zbar(k)=z;eos%yh(k)=a*yi
          eos%abar_present(k)=.true.;eos%zbar_present(k)=.true.;eos%yh_present(k)=.true.
          eos%composition_available=.true.
          pos=pos+4
        enddo
        if(pos/=nt+1) then;close(u);status='trailing fields in eos.compo row';return;endif
        seen(k)=.true.
      enddo
      close(u)
      if(.not.all(seen)) then;status='eos.compo does not cover the thermo grid';return;endif
      ! M2A: per-row status + C3/C4 applicability. Phase-only rows (Npairs=0
      ! and Nquad=0) are PARTIAL (1), not FAIL. The old global range check is
      ! replaced by a per-row range check that sets row_status=FAIL only on the
      ! offending rows; the read still succeeds so partial-coverage sources
      ! (e.g. GM_GM1 outer-crust phase-only rows) do not abort the whole read.
      do k=1,eos%n
        if(.not.seen(k)) then
          eos%row_status(k)=0  ! NOT_AVAILABLE
          cycle
        endif
        ! per-row range: a written fraction outside [-1e-10, 1+5e-3] -> FAIL.
        if((eos%yn_present(k).and.eos%yn(k)< -1.e-10_wp).or. &
           (eos%yp_present(k).and.eos%yp(k)< -1.e-10_wp).or. &
           (eos%ye_present(k).and.eos%ye(k)< -1.e-10_wp).or. &
           (eos%ymu_present(k).and.eos%ymu(k)< -1.e-10_wp).or. &
           (eos%yh_present(k).and.eos%yh(k)>1._wp+5.e-3_wp)) then
          eos%row_status(k)=3  ! FAIL
          cycle
        endif
        ! C3: complete baryon inventory for Yn+Yp+Yh=1 (Yn, Yp, Yh all present).
        eos%c3_applicable(k) = eos%yn_present(k).and.eos%yp_present(k).and.eos%yh_present(k)
        ! C4: complete charge inventory for Yp+(Z/A)Yh=Ye+Ymu (Yp, Abar>0, Zbar, Yh, Ye, Ymu).
        eos%c4_applicable(k) = eos%yp_present(k).and.eos%abar_present(k).and. &
                              eos%zbar_present(k).and.eos%yh_present(k).and. &
                              eos%ye_present(k).and.eos%ymu_present(k).and. &
                              eos%abar(k)>0._wp
        ! C3/C4 applicability is a data-coverage statement, not a closure
        ! verdict. CompOSE permits relevant particle subsets to be omitted;
        ! absent pairs are not zero. Keep row_status as the M2A input-range/
        ! coverage label; the verifier evaluates closure only where applicable.
        eos%c3_applicable(k)=eos%yn_present(k).and.eos%yp_present(k).and.eos%yh_present(k)
        eos%c4_applicable(k)=eos%yp_present(k).and.eos%abar_present(k).and. &
          eos%zbar_present(k).and.eos%yh_present(k).and.eos%ye_present(k).and. &
          eos%ymu_present(k).and.eos%abar(k)>0._wp
        if(eos%c3_applicable(k).and.eos%c4_applicable(k))then
          eos%row_status(k)=2
        else
          eos%row_status(k)=1
        endif
      enddo
    else
      eos%composition_available=.false.
      ! M2A: no compo file -> every row NOT_AVAILABLE (0); stays at default.
    endif
    eos%beta_equilibrated=maxval(abs(eos%q(5,:)))<1.e-8_wp
    status='ok'
  end subroutine read_compose

  subroutine write_internal_csv(path,eos,status)
    character(len=*),intent(in)::path
    type(eos_data_t),intent(in)::eos
    character(len=*),intent(out)::status
    integer::u,ios,i,j,iy,nyq_out,first_iy,last_iy
    real(wp)::v,qv(7)
    character(len=512)::header
    character(len=64)::colname
    open(newunit=u,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0) then;status='cannot create internal EOS CSV';return;endif
    ! E390-10B phase 2: emit one q_extra column per declared extra quantity
    ! (q_extra for slot 1, q_extra_2..q_extra_{nadd_max} for slots 2+). Legacy
    ! readers that only consume q_extra (slot 1) keep working; the additional
    ! columns are optional and default to 0 on read if absent.
    ! E390-10B phase 3: when nYq>1, prefix every row with an iy_local column so
    ! the multidim Yq blocks can be round-tripped. nYq==1 (the legacy cold
    ! slice) omits iy_local entirely, keeping the CSV shape unchanged.
    nyq_out=eos%nYq
    if(nyq_out<1) nyq_out=1
    if(nyq_out>1) then
      header='iy_local,nb,Q1,Q2,Q3,Q4,Q5,Q6,Q7,q_extra'
    else
      header='nb,Q1,Q2,Q3,Q4,Q5,Q6,Q7,q_extra'
    endif
    do j=2,max(1,eos%nadd_max)
      write(colname,'(a,i0)')',q_extra_',j
      header=trim(header)//trim(colname)
    enddo
    header=trim(header)//',Iphase,Yn,Yp,Ye,Ymu,Abar,Zbar,Yh,'// &
      'row_status,yn_present,yp_present,ye_present,ymu_present,abar_present,zbar_present,yh_present,c3_applicable,c4_applicable'
    write(u,'(a)')trim(header)
    first_iy=1
    last_iy=nyq_out
      do iy=first_iy,last_iy
      do i=1,eos%n
        ! E390-10B phase 3: for nYq>1 pull Q1..Q7 from q_full(1,iy,:,i); the 1D
        ! q(7,n) surface only holds the selected (first) Yq. q_extra slot 1 is
        ! taken from q_extra_full(1,1,iy,i) for the same reason.
        if(allocated(eos%q_full).and.eos%nYq>1) then
          qv=eos%q_full(1,iy,:,i)
        else
          qv=eos%q(:,i)
        endif
        if(allocated(eos%q_extra_full).and.eos%nYq>1.and.eos%nadd_max>=1) then
          v=eos%q_extra_full(1,1,iy,i)
        else
          v=eos%q_extra(i)
        endif
        ! base row: (iy_local,) nb, Q1..Q7, q_extra (slot 1)
        if(nyq_out>1) then
          write(u,'(i0,",",es25.16e3,7(",",es25.16e3),",",es25.16e3)',advance='no') &
            iy,eos%nb(i),qv,v
        else
          write(u,'(es25.16e3,7(",",es25.16e3),",",es25.16e3)',advance='no') &
            eos%nb(i),qv,v
        endif
        ! additional extra slots 2..nadd_max (from the tensor; 0 if unallocated).
        ! For nYq>1 pull from the tensor's iy slot; the 1D q_extra holds slot 1
        ! of the selected (first) Yq only.
        do j=2,max(1,eos%nadd_max)
          if(allocated(eos%q_extra_full).and.size(eos%q_extra_full,1)>=j) then
            write(u,'(",",es25.16e3)',advance='no')eos%q_extra_full(j,1,iy,i)
          else
            write(u,'(",",es25.16e3)',advance='no')0._wp
          endif
        enddo
        ! composition/masks are only meaningful for the selected (first) Yq;
        ! the other Yq blocks emit zeros (they carry thermo only).
        if(iy==1) then
          write(u,'(",",i0,7(",",es25.16e3),",",i0,9(",",i0))') &
            eos%iphase(i),eos%yn(i),eos%yp(i),eos%ye(i), &
            eos%ymu(i),eos%abar(i),eos%zbar(i),eos%yh(i),eos%row_status(i), &
            merge(1,0,eos%yn_present(i)),merge(1,0,eos%yp_present(i)), &
            merge(1,0,eos%ye_present(i)),merge(1,0,eos%ymu_present(i)), &
            merge(1,0,eos%abar_present(i)),merge(1,0,eos%zbar_present(i)), &
            merge(1,0,eos%yh_present(i)),merge(1,0,eos%c3_applicable(i)),merge(1,0,eos%c4_applicable(i))
        else
          write(u,'(",",i0,7(",",es25.16e3),",",i0,9(",",i0))') &
            0,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0,0,0,0,0,0,0,0,0,0
        endif
      enddo
    enddo
    close(u);status='ok'
  end subroutine write_internal_csv

  subroutine read_internal_csv(path,meta_path,eos,status)
    character(len=*),intent(in)::path,meta_path
    type(eos_data_t),intent(out)::eos
    character(len=*),intent(out)::status
    integer::u,ios,n,i,nt,k,ival,j,nqextra,eos_n
    integer::iy_local,kk
    character(len=4096)::line
    character(len=96)::tok(64)
    character(len=96)::hcol
    integer::col_idx(64)  ! maps our known-field code -> token index; 0 if absent
    character(len=96)::hdr(64)
    integer::nhdr
    logical::ok,found
    real(wp)::v
    integer::phase
    logical::composition_file_present,composition_available,beta_equilibrated
    integer::il,iT,iYq,inb_start,nadd,nadd_max,naxis_t,naxis_yq
    real(wp)::mn,mp,temperature_mev,charge_fraction_grid
    namelist /metadata/ mn,mp,il,iT,iYq,inb_start,nadd,nadd_max,naxis_t,naxis_yq,temperature_mev,charge_fraction_grid, &
       composition_file_present,composition_available,beta_equilibrated
    open(newunit=u,file=trim(meta_path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='cannot open internal EOS metadata';return;endif
    mn=939._wp;mp=939._wp;il=1;iT=1;iYq=1;inb_start=1;nadd=0;nadd_max=0;naxis_t=1;naxis_yq=1
    temperature_mev=0._wp;charge_fraction_grid=0._wp
    composition_file_present=.false.;composition_available=.false.;beta_equilibrated=.false.
    read(u,nml=metadata,iostat=ios);close(u)
    if(ios/=0) then;status='invalid internal EOS metadata namelist';return;endif
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='cannot open internal EOS CSV';return;endif
    call next_line(u,line,ios)
    if(ios/=0) then;close(u);status='missing internal EOS CSV header';return;endif
    ! E390-10B phase 2: header-driven column lookup. Parse the header line into
    ! names and build col_idx(code) = token position. This supports a variable
    ! number of q_extra_N columns and is backward compatible with legacy 17/27
    ! field CSVs (absent optional columns default to 0 / .false.).
    call split_tokens(line,hdr,nhdr)
    col_idx=0
    do k=1,nhdr
      hcol=adjustl(hdr(k))
      ! normalize to lowercase so the header matches the case-insensitive
      ! Fortran namelist convention and the lowercase mask column names the
      ! Python writer emits (yn_present, ...). Q1..Q7 and Abar/Zbar/Yh are
      ! also matched in lowercase here for robustness.
      do j=1,len_trim(hcol)
        if(hcol(j:j)>='A'.and.hcol(j:j)<='Z') then
          hcol(j:j)=achar(iachar(hcol(j:j))+32)
        endif
      enddo
      select case(trim(hcol))
      case('nb');col_idx(1)=k
      case('q1');col_idx(2)=k
      case('q2');col_idx(3)=k
      case('q3');col_idx(4)=k
      case('q4');col_idx(5)=k
      case('q5');col_idx(6)=k
      case('q6');col_idx(7)=k
      case('q7');col_idx(8)=k
      case('q_extra');col_idx(9)=k
      case('iphase');col_idx(10)=k
      case('yn');col_idx(11)=k
      case('yp');col_idx(12)=k
      case('ye');col_idx(13)=k
      case('ymu');col_idx(14)=k
      case('abar');col_idx(15)=k
      case('zbar');col_idx(16)=k
      case('yh');col_idx(17)=k
      case('row_status');col_idx(18)=k
      case('yn_present');col_idx(19)=k
      case('yp_present');col_idx(20)=k
      case('ye_present');col_idx(21)=k
      case('ymu_present');col_idx(22)=k
      case('abar_present');col_idx(23)=k
      case('zbar_present');col_idx(24)=k
      case('yh_present');col_idx(25)=k
      case('c3_applicable');col_idx(26)=k
      case('c4_applicable');col_idx(27)=k
      case('iy_local');col_idx(28)=k
      end select
    enddo
    ! required columns: nb + 7 Q + q_extra + Iphase + 7 compo (codes 1..17)
    do k=1,17
      if(col_idx(k)==0)then;close(u);status='internal EOS CSV missing required column';return;endif
    enddo
    ! count extra q_extra_N columns present in the header (N>=2)
    nqextra=1
    do j=2,64
      write(hcol,'(a,i0)')'q_extra_',j
      found=.false.
      do k=1,nhdr
        if(trim(adjustl(hdr(k)))==trim(hcol))then;found=.true.;exit;endif
      enddo
      if(.not.found)exit
      nqextra=j
    enddo
    n=0
    do
      call next_line(u,line,ios)
      if(ios<0)exit
      if(ios>0)then;close(u);status='error counting internal EOS rows';return;endif
      n=n+1
    enddo
    close(u)
    ! E390-10B phase 3: when iy_local is present the CSV carries nYq*n rows (one
    ! per (Yq,density) pair); the density grid size is n/nYq. Otherwise (legacy
    ! cold slice) the row count is already the density count.
    if(col_idx(28)/=0) then
      if(max(naxis_yq,1)<=0) then;status='naxis_yq must be positive for iy_local CSV';return;endif
      if(mod(n,max(naxis_yq,1))/=0) then;status='iy_local CSV row count is not nYq*n';return;endif
      eos_n=n/max(naxis_yq,1)
    else
      eos_n=n
    endif
    call allocate_eos(eos,eos_n)
    eos%mn=mn;eos%mp=mp;eos%il=il;eos%iT=iT;eos%iYq=iYq;eos%inb_start=inb_start;eos%nadd=nadd
    eos%nadd_max=max(nadd_max,nqextra)
    eos%nT=max(naxis_t,1);eos%nYq=max(naxis_yq,1)
    eos%temperature_mev=temperature_mev;eos%charge_fraction_grid=charge_fraction_grid
    eos%composition_file_present=composition_file_present;eos%composition_available=composition_available
    eos%beta_equilibrated=beta_equilibrated
    ! E390-10B phase 2/3: allocate q_extra_full when the table carries extra
    ! columns (any nadd_max>=1). E390-10B phase 3: allocate q_full when nYq>1.
    if(eos%nadd_max>=1) then
      if(.not.allocated(eos%q_extra_full)) then
        allocate(eos%q_extra_full(eos%nadd_max,eos%nT,eos%nYq,eos_n))
        eos%q_extra_full=0._wp
      endif
    endif
    if(eos%nYq>1) then
      if(.not.allocated(eos%q_full)) then
        allocate(eos%q_full(eos%nT,eos%nYq,7,eos_n))
        eos%q_full=0._wp
      endif
    endif
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='cannot reopen internal EOS CSV';return;endif
    call next_line(u,line,ios)
    do i=1,n
      call next_line(u,line,ios)
      if(ios/=0)then;close(u);status='unexpected EOF in internal EOS CSV';return;endif
      call split_tokens(line,tok,nt)
      ! E390-10B phase 3: when iy_local is present, each row is a (Yq,density)
      ! pair; compute the density index kk and the Yq index iy_local. The 1D
      ! surface (q/q_extra/compo/masks) only materializes the selected (first)
      ! Yq (iy_local==1); the other Yq blocks fill the q_full/q_extra_full
      ! tensors. Without iy_local, kk==i (legacy cold slice).
      if(col_idx(28)/=0.and.col_idx(28)<=nt) then
        call token_int(tok,col_idx(28),iy_local,ok);if(.not.ok)then;close(u);status='invalid iy_local';return;endif
        kk=mod(i-1,eos_n)+1
      else
        iy_local=1;kk=i
      endif
      ! required thermo+compo columns 1..17 via header lookup
      do k=1,8
        if(col_idx(k)==0.or.col_idx(k)>nt)cycle
        call token_real(tok,col_idx(k),v,ok);if(.not.ok)then;close(u);status='invalid internal EOS value';return;endif
        if(k==1)then
          if(iy_local==1) eos%nb(kk)=v
        else
          if(iy_local==1) eos%q(k-1,kk)=v
          if(allocated(eos%q_full)) eos%q_full(1,iy_local,k-1,kk)=v
        endif
      enddo
      if(col_idx(9)/=0.and.col_idx(9)<=nt) then
        call token_real(tok,col_idx(9),v,ok);if(.not.ok)then;close(u);status='invalid q_extra';return;endif
        if(iy_local==1) eos%q_extra(kk)=v
        if(allocated(eos%q_extra_full)) eos%q_extra_full(1,1,iy_local,kk)=v
      endif
      ! additional extra slots 2..nqextra (header-named q_extra_2..)
      do j=2,nqextra
        write(hcol,'(a,i0)')'q_extra_',j
        found=.false.
        do k=1,nhdr
          if(trim(adjustl(hdr(k)))==trim(hcol))then;found=.true.;exit;endif
        enddo
        if(found.and.k<=nt) then
          call token_real(tok,k,v,ok);if(.not.ok)then;close(u);status='invalid q_extra_N';return;endif
          if(allocated(eos%q_extra_full).and.size(eos%q_extra_full,1)>=j) &
            eos%q_extra_full(j,1,iy_local,kk)=v
        endif
      enddo
      if(col_idx(10)/=0.and.col_idx(10)<=nt.and.iy_local==1) then
        call token_int(tok,col_idx(10),phase,ok);if(.not.ok)then;close(u);status='invalid Iphase';return;endif
        eos%iphase(kk)=phase
      endif
      do k=11,17
        if(col_idx(k)==0.or.col_idx(k)>nt)cycle
        if(iy_local/=1) cycle
        call token_real(tok,col_idx(k),v,ok);if(.not.ok)then;close(u);status='invalid internal composition field';return;endif
        select case(k)
        case(11);eos%yn(kk)=v
        case(12);eos%yp(kk)=v
        case(13);eos%ye(kk)=v
        case(14);eos%ymu(kk)=v
        case(15);eos%abar(kk)=v
        case(16);eos%zbar(kk)=v
        case(17);eos%yh(kk)=v
        end select
      enddo
      ! optional M2A mask columns 18..27 (absent in legacy 17-field CSV; only
      ! meaningful for the selected Yq block)
      if(col_idx(18)/=0.and.col_idx(18)<=nt.and.iy_local==1) then
        call token_int(tok,col_idx(18),ival,ok);if(.not.ok)then;close(u);status='invalid row_status';return;endif
        eos%row_status(kk)=ival
        call token_int(tok,col_idx(19),ival,ok);if(.not.ok)then;close(u);status='invalid yn_present';return;endif
        eos%yn_present(kk)=ival/=0
        call token_int(tok,col_idx(20),ival,ok);if(.not.ok)then;close(u);status='invalid yp_present';return;endif
        eos%yp_present(kk)=ival/=0
        call token_int(tok,col_idx(21),ival,ok);if(.not.ok)then;close(u);status='invalid ye_present';return;endif
        eos%ye_present(kk)=ival/=0
        call token_int(tok,col_idx(22),ival,ok);if(.not.ok)then;close(u);status='invalid ymu_present';return;endif
        eos%ymu_present(kk)=ival/=0
        call token_int(tok,col_idx(23),ival,ok);if(.not.ok)then;close(u);status='invalid abar_present';return;endif
        eos%abar_present(kk)=ival/=0
        call token_int(tok,col_idx(24),ival,ok);if(.not.ok)then;close(u);status='invalid zbar_present';return;endif
        eos%zbar_present(kk)=ival/=0
        call token_int(tok,col_idx(25),ival,ok);if(.not.ok)then;close(u);status='invalid yh_present';return;endif
        eos%yh_present(kk)=ival/=0
        call token_int(tok,col_idx(26),ival,ok);if(.not.ok)then;close(u);status='invalid c3_applicable';return;endif
        eos%c3_applicable(kk)=ival/=0
        call token_int(tok,col_idx(27),ival,ok);if(.not.ok)then;close(u);status='invalid c4_applicable';return;endif
        eos%c4_applicable(kk)=ival/=0
      endif
    enddo
    close(u);status='ok'
  end subroutine read_internal_csv

  subroutine write_metadata(path,eos,status)
    character(len=*),intent(in)::path
    type(eos_data_t),intent(in)::eos
    character(len=*),intent(out)::status
    integer::u,ios
    open(newunit=u,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write metadata';return;endif
    write(u,'(a)')'&metadata'
    write(u,'(a,es25.16e3)')' mn=',eos%mn
    write(u,'(a,es25.16e3)')' mp=',eos%mp
    write(u,'(a,i0)')' il=',eos%il
    write(u,'(a,i0)')' iT=',eos%iT
    write(u,'(a,i0)')' iYq=',eos%iYq
    write(u,'(a,i0)')' inb_start=',eos%inb_start
    write(u,'(a,i0)')' nadd=',eos%nadd
    write(u,'(a,i0)')' nadd_max=',eos%nadd_max
    write(u,'(a,i0)')' naxis_t=',eos%nT
    write(u,'(a,i0)')' naxis_yq=',eos%nYq
    write(u,'(a,es25.16e3)')' temperature_mev=',eos%temperature_mev
    write(u,'(a,es25.16e3)')' charge_fraction_grid=',eos%charge_fraction_grid
    write(u,'(a,l1)')' composition_file_present=',eos%composition_file_present
    write(u,'(a,l1)')' composition_available=',eos%composition_available
    write(u,'(a,l1)')' beta_equilibrated=',eos%beta_equilibrated
    write(u,'(a)')'/'
    close(u);status='ok'
  end subroutine write_metadata

  subroutine read_metadata(path,eos,status)
    character(len=*),intent(in)::path
    type(eos_data_t),intent(inout)::eos
    character(len=*),intent(out)::status
    integer::u,ios
    integer::il,iT,iYq,inb_start,nadd
    real(wp)::mn,mp,temperature_mev,charge_fraction_grid
    logical::composition_file_present,composition_available,beta_equilibrated
    namelist /metadata/ mn,mp,il,iT,iYq,inb_start,nadd,temperature_mev,charge_fraction_grid, &
       composition_file_present,composition_available,beta_equilibrated
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0)then;status='cannot open metadata';return;endif
    mn=eos%mn;mp=eos%mp;il=1;iT=1;iYq=1;inb_start=1;nadd=0
    temperature_mev=0._wp;charge_fraction_grid=0._wp
    composition_file_present=.false.;composition_available=.false.;beta_equilibrated=.false.
    read(u,nml=metadata,iostat=ios);close(u)
    if(ios/=0)then;status='invalid metadata namelist';return;endif
    eos%mn=mn;eos%mp=mp;eos%il=il;eos%iT=iT;eos%iYq=iYq;eos%inb_start=inb_start;eos%nadd=nadd
    eos%temperature_mev=temperature_mev;eos%charge_fraction_grid=charge_fraction_grid
    eos%composition_file_present=composition_file_present;eos%composition_available=composition_available
    eos%beta_equilibrated=beta_equilibrated;status='ok'
  end subroutine read_metadata

  subroutine write_normalized_csv(path,eos,status)
    character(len=*),intent(in)::path
    type(eos_data_t),intent(in)::eos
    character(len=*),intent(out)::status
    integer::u,ios,i
    real(wp),parameter::fm3_to_m3=1.e45_wp,mev_fm3_to_si=1.602176634e32_wp,mev_to_j=1.602176634e-13_wp
    real(wp)::p,e,mu
    open(newunit=u,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot create normalized EOS CSV';return;endif
    write(u,'(a)')'source_row,i_t,i_nb,i_yq,n_b_fm3,n_b_m3,pressure_mev_fm3,pressure_pa,energy_mev_fm3,energy_j_m3,mu_b_mev,mu_b_j,q1_p_over_nb_mev,q2_entropy_per_baryon,q3_mu_b_over_mn_minus_1,q4_mu_q_over_mn,q5_mu_l_over_mn,q6_f_over_nb_mn_minus_1,q7_e_over_nb_mn_minus_1'
    do i=1,eos%n
      p=eos%q(1,i)*eos%nb(i);e=eos%nb(i)*(eos%q(7,i)+1._wp)*eos%mn
      mu=(eos%q(3,i)+1._wp)*eos%mn
      write(u,'(4(i0,","),14(es25.16e3,","),es25.16e3)') i,eos%iT,eos%inb_start+i-1,eos%iYq, &
        eos%nb(i),eos%nb(i)*fm3_to_m3,p,p*mev_fm3_to_si,e,e*mev_fm3_to_si, &
        mu,mu*mev_to_j,eos%q(:,i)
    enddo
    close(u);status='ok'
  end subroutine write_normalized_csv

  subroutine write_compose(directory,eos,status)
    character(len=*),intent(in)::directory
    type(eos_data_t),intent(in)::eos
    character(len=*),intent(out)::status
    integer::u,ios,i,nquad,npairs,t0_out,t1_out,yq0_out,yq1_out
    logical::any_composition,group_any,group_complete
    character(len=1024)::path
    real(wp)::yi,baryon_residual,charge_residual,max_baryon_residual,max_charge_residual
    ! M10 safety gate: the current writer serializes one selected T/Yq slice.
    ! Refuse a multidimensional EOS rather than silently dropping tensor cells.
    if(eos%nT>1.or.eos%nYq>1)then
      status='multidimensional CompOSE writer is not implemented; refusing lossy serialization'
      return
    endif
    if(eos%n<1.or..not.allocated(eos%iphase)) then
      status='cannot write invalid EOS phase metadata'
      return
    endif
    if(size(eos%iphase)/=eos%n) then
      status='cannot write inconsistent EOS phase metadata'
      return
    endif
    if(.not.allocated(eos%yn_present).or..not.allocated(eos%yp_present).or. &
       .not.allocated(eos%ye_present).or..not.allocated(eos%ymu_present).or. &
       .not.allocated(eos%abar_present).or..not.allocated(eos%zbar_present).or. &
       .not.allocated(eos%yh_present).or..not.allocated(eos%c3_applicable).or. &
       .not.allocated(eos%c4_applicable).or..not.allocated(eos%row_status))then
      status='row-level composition masks are missing; refusing lossy CompOSE export'
      return
    endif
    if(size(eos%yn_present)/=eos%n.or.size(eos%yp_present)/=eos%n.or. &
       size(eos%ye_present)/=eos%n.or.size(eos%ymu_present)/=eos%n.or. &
       size(eos%abar_present)/=eos%n.or.size(eos%zbar_present)/=eos%n.or. &
       size(eos%yh_present)/=eos%n.or.size(eos%c3_applicable)/=eos%n.or. &
       size(eos%c4_applicable)/=eos%n.or.size(eos%row_status)/=eos%n)then
      status='row-level composition mask shape mismatch'
      return
    endif
    any_composition=any(eos%yn_present.or.eos%yp_present.or.eos%ye_present.or. &
                        eos%ymu_present.or.eos%abar_present.or.eos%zbar_present.or.eos%yh_present)
    if(any_composition.and..not.eos%composition_file_present)then
      status='present composition fields require an eos.compo sidecar'
      return
    endif
    if(eos%composition_available.and..not.any_composition)then
      status='global composition_available conflicts with all row-level presence masks being false'
      return
    endif
    if(any(eos%iphase/=0).and..not.eos%composition_file_present)then
      status='cannot export a phase-discontinuous barotrope without phase composition'
      return
    endif
    if(eos%composition_file_present)then
      max_baryon_residual=0._wp;max_charge_residual=0._wp
      do i=1,eos%n
        if(eos%row_status(i)==3)then
          status='row-level composition status is FAIL; refusing CompOSE export'
          return
        endif
        group_any=eos%abar_present(i).or.eos%zbar_present(i).or.eos%yh_present(i)
        group_complete=eos%abar_present(i).and.eos%zbar_present(i).and.eos%yh_present(i)
        if(group_any.neqv.group_complete)then
          status='partial nuclear-group tuple cannot be represented in CompOSE eos.compo'
          return
        endif
        if(eos%yn_present(i))then
          if(.not.ieee_is_finite(eos%yn(i)).or.eos%yn(i)<-1.e-10_wp.or.eos%yn(i)>1._wp+5.e-3_wp)then
            status='neutron fraction outside physical range';return
          endif
        else if(abs(eos%yn(i))>1.e-14_wp)then
          status='unmarked neutron value would be lost during CompOSE export';return
        endif
        if(eos%yp_present(i))then
          if(.not.ieee_is_finite(eos%yp(i)).or.eos%yp(i)<-1.e-10_wp.or.eos%yp(i)>1._wp+5.e-3_wp)then
            status='proton fraction outside physical range';return
          endif
        else if(abs(eos%yp(i))>1.e-14_wp)then
          status='unmarked proton value would be lost during CompOSE export';return
        endif
        if(eos%ye_present(i))then
          if(.not.ieee_is_finite(eos%ye(i)).or.eos%ye(i)<-1.e-10_wp.or.eos%ye(i)>1._wp+5.e-3_wp)then
            status='electron fraction outside physical range';return
          endif
        else if(abs(eos%ye(i))>1.e-14_wp)then
          status='unmarked electron value would be lost during CompOSE export';return
        endif
        if(eos%ymu_present(i))then
          if(.not.ieee_is_finite(eos%ymu(i)).or.eos%ymu(i)<-1.e-10_wp.or.eos%ymu(i)>1._wp+5.e-3_wp)then
            status='muon fraction outside physical range';return
          endif
        else if(abs(eos%ymu(i))>1.e-14_wp)then
          status='unmarked muon value would be lost during CompOSE export';return
        endif
        if(group_complete)then
          if(.not.all(ieee_is_finite([eos%abar(i),eos%zbar(i),eos%yh(i)])))then
            status='nuclear-group composition contains nonfinite values';return
          endif
          if(eos%yh(i)<-1.e-10_wp.or.eos%yh(i)>1.002_wp.or. &
             eos%zbar(i)<0._wp.or.eos%zbar(i)>eos%abar(i).or. &
             (eos%abar(i)>0._wp.and.eos%abar(i)<1._wp))then
            status='nuclear-group composition outside physical range';return
          endif
          if(eos%abar(i)>0._wp)then
            yi=eos%yh(i)/eos%abar(i)
          else
            if(abs(eos%zbar(i))>1.e-10_wp.or.abs(eos%yh(i))>1.e-10_wp)then
              status='nuclear-group values are nonzero although Abar is zero';return
            endif
            yi=0._wp
          endif
        else
          yi=0._wp
          if(abs(eos%abar(i))+abs(eos%zbar(i))+abs(eos%yh(i))>1.e-14_wp)then
            status='unmarked nuclear-group values would be lost during CompOSE export';return
          endif
        endif
        if(eos%c3_applicable(i))then
          if(.not.(eos%yn_present(i).and.eos%yp_present(i).and.eos%yh_present(i)))then
            status='C3 marked applicable without complete baryon inventory';return
          endif
          baryon_residual=eos%yn(i)+eos%yp(i)+eos%yh(i)-1._wp
          max_baryon_residual=max(max_baryon_residual,abs(baryon_residual))
          if(abs(baryon_residual)>2.e-3_wp)then
            status='C3 baryon abundance closure failed on an applicable row';return
          endif
        endif
        if(eos%c4_applicable(i))then
          if(.not.(eos%yp_present(i).and.eos%abar_present(i).and.eos%zbar_present(i).and. &
             eos%yh_present(i).and.eos%ye_present(i).and.eos%ymu_present(i)).or.eos%abar(i)<=0._wp)then
            status='C4 marked applicable without complete charge inventory';return
          endif
          charge_residual=eos%yp(i)+(eos%zbar(i)/eos%abar(i))*eos%yh(i)-eos%ye(i)-eos%ymu(i)
          max_charge_residual=max(max_charge_residual,abs(charge_residual))
          if(abs(charge_residual)>2.e-3_wp)then
            status='C4 charge abundance closure failed on an applicable row';return
          endif
        endif
      enddo
    else
      if(any(abs(eos%yn)>1.e-14_wp).or.any(abs(eos%yp)>1.e-14_wp).or. &
         any(abs(eos%ye)>1.e-14_wp).or.any(abs(eos%ymu)>1.e-14_wp).or. &
         any(abs(eos%abar)>1.e-14_wp).or.any(abs(eos%zbar)>1.e-14_wp).or. &
         any(abs(eos%yh)>1.e-14_wp))then
        status='composition values are present without an eos.compo sidecar'
        return
      endif
    endif
    t0_out=eos%t_first;t1_out=eos%t_last
    yq0_out=eos%yq_first;yq1_out=eos%yq_last
    ! Preserve source axis labels when known; internal-only tables fall back
    ! to their selected cold-slice indices.
    if(t1_out<t0_out)then;t0_out=eos%iT;t1_out=eos%iT;endif
    if(yq1_out<yq0_out)then;yq0_out=eos%iYq;yq1_out=eos%iYq;endif
    path=trim(directory)//'/eos.nb'
    open(newunit=u,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write eos.nb';return;endif
    write(u,'(i0)')eos%inb_start;write(u,'(i0)')eos%inb_start+eos%n-1
    do i=1,eos%n;write(u,'(es25.16e3)')eos%nb(i);enddo
    close(u)
    open(newunit=u,file=trim(directory)//'/eos.t',status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write eos.t';return;endif
    write(u,'(i0)')t0_out;write(u,'(i0)')t1_out;write(u,'(es25.16e3)')eos%temperature_mev;close(u)
    open(newunit=u,file=trim(directory)//'/eos.yq',status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write eos.yq';return;endif
    write(u,'(i0)')yq0_out;write(u,'(i0)')yq1_out;write(u,'(es25.16e3)')eos%charge_fraction_grid;close(u)
    open(newunit=u,file=trim(directory)//'/eos.thermo',status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write eos.thermo';return;endif
    write(u,'(2(es25.16e3,1x),i0)')eos%mn,eos%mp,eos%il
    do i=1,eos%n
      ! E390-10B phase 2: write every extra column. The 1D q_extra(n) surface
      ! holds only the first column for backward compatibility, so for Nadd>=2
      ! pull the additional columns from q_extra_full (cold slice: nT=nYq=1).
      if(eos%nadd==0) then
        write(u,'(3(i0,1x),7(es25.16e3,1x),i0)')eos%iT,eos%inb_start+i-1,eos%iYq,eos%q(:,i),0
      else if(eos%nadd==1) then
        write(u,'(3(i0,1x),7(es25.16e3,1x),i0,1x,es25.16e3)')eos%iT,eos%inb_start+i-1,eos%iYq,eos%q(:,i),1,eos%q_extra(i)
      else
        ! Nadd>=2: emit the full row from the q_extra_full tensor (cold slice).
        block
          integer::j
          character(len=4096)::row
          character(len=64)::buf
          write(row,'(3(i0,1x),7(es25.16e3,1x),i0)')eos%iT,eos%inb_start+i-1,eos%iYq,eos%q(:,i),eos%nadd
          do j=1,eos%nadd
            if(allocated(eos%q_extra_full).and.size(eos%q_extra_full,1)>=j) then
              write(buf,'(1x,es25.16e3)')eos%q_extra_full(j,eos%nT,eos%nYq,i)
            else
              write(buf,'(1x,es25.16e3)')0._wp
            endif
            row=trim(row)//trim(buf)
          enddo
          write(u,'(a)')trim(row)
        end block
      endif
    enddo
    close(u)
    if(eos%composition_file_present)then
      open(newunit=u,file=trim(directory)//'/eos.compo',status='replace',action='write',iostat=ios)
      if(ios/=0)then;status='cannot write eos.compo';return;endif
      do i=1,eos%n
        npairs=0
        if(eos%yn_present(i))npairs=npairs+1
        if(eos%yp_present(i))npairs=npairs+1
        if(eos%ye_present(i))npairs=npairs+1
        if(eos%ymu_present(i))npairs=npairs+1
        group_complete=eos%abar_present(i).and.eos%zbar_present(i).and.eos%yh_present(i)
        nquad=merge(1,0,group_complete)
        write(u,'(5(i0,1x))',advance='no')eos%iT,eos%inb_start+i-1,eos%iYq,eos%iphase(i),npairs
        if(eos%yn_present(i))write(u,'(i0,1x,es25.16e3,1x)',advance='no')10,eos%yn(i)
        if(eos%yp_present(i))write(u,'(i0,1x,es25.16e3,1x)',advance='no')11,eos%yp(i)
        if(eos%ye_present(i))write(u,'(i0,1x,es25.16e3,1x)',advance='no')0,eos%ye(i)
        if(eos%ymu_present(i))write(u,'(i0,1x,es25.16e3,1x)',advance='no')1,eos%ymu(i)
        write(u,'(i0)',advance='no')nquad
        if(nquad==1)then
          if(eos%abar(i)>0._wp)then
            yi=eos%yh(i)/eos%abar(i)
          else
            yi=0._wp
          endif
          write(u,'(1x,i0,3(1x,es25.16e3))')1,eos%abar(i),eos%zbar(i),yi
        else
          write(u,'(a)')''
        endif
      enddo
      close(u)
    else
      ! Avoid a stale composition sidecar surviving a thermo-only rewrite.
      open(newunit=u,file=trim(directory)//'/eos.compo',status='old',action='read',iostat=ios)
      if(ios==0)close(u,status='delete')
    endif
    status='ok'
  end subroutine write_compose
end module compose_eos_io
