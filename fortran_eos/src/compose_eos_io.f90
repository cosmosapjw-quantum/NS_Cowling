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
    integer::nb0,nb1,t0,t1,y0,y1,u,ios,i,k,nt,it,ib,iy,nadd,extraid,npairs,nquad,pos,pid,qid
    character(len=4096)::line
    character(len=96)::tok(256)
    real(wp),allocatable::nb(:),tv(:),yv(:)
    real(wp)::qv(7),v,a,z,yi
    logical::ok,file_exists
    logical,allocatable::seen(:)
    character(len=1024)::path
    eos=eos_data_t();status='uninitialized'
    path=trim(directory)//'/eos.nb'
    call read_grid(trim(path),nb0,nb1,nb,status);if(status/='ok')return
    call read_grid(trim(directory)//'/eos.t',t0,t1,tv,status);if(status/='ok')return
    call read_grid(trim(directory)//'/eos.yq',y0,y1,yv,status);if(status/='ok')return
    if(size(nb)<3.or.any(nb<=0._wp).or.any(nb(2:)<=nb(:size(nb)-1))) then
      status='eos.nb must contain at least three positive strictly increasing values';return
    endif
    if(size(tv)/=1.or.size(yv)/=1.or.abs(tv(1))>1.e-12_wp) then
      status='only a single cold T=0 eos.t/eos.yq slice is supported';return
    endif
    call allocate_eos(eos,size(nb));eos%nb=nb;eos%inb_start=nb0;eos%temperature_mev=tv(1)
    eos%charge_fraction_grid=yv(1)
    path=trim(directory)//'/eos.thermo'
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='missing or unreadable eos.thermo';return;endif
    call next_line(u,line,ios)
    if(ios/=0) then;close(u);status='empty eos.thermo';return;endif
    read(line,*,iostat=ios)eos%mn,eos%mp,eos%il
    if(ios/=0.or.eos%mn<=0._wp.or.eos%mp<=0._wp.or.(eos%il/=0.and.eos%il/=1)) then
      close(u);status='invalid eos.thermo header';return
    endif
    allocate(seen(eos%n));seen=.false.;eos%nadd=-1
    do
      call next_line(u,line,ios)
      if(ios<0) exit
      if(ios>0) then;close(u);status='error reading eos.thermo';return;endif
      call split_tokens(line,tok,nt)
      if(nt<11) then;close(u);status='short eos.thermo row';return;endif
      call token_int(tok,1,it,ok);if(.not.ok) then;close(u);status='invalid iT';return;endif
      call token_int(tok,2,ib,ok);if(.not.ok) then;close(u);status='invalid inb';return;endif
      call token_int(tok,3,iy,ok);if(.not.ok) then;close(u);status='invalid iYq';return;endif
      if(it<t0.or.it>t1.or.iy<y0.or.iy>y1.or.it/=t0.or.iy/=y0) then
        close(u);status='thermo table is not a single selected T/Yq slice';return
      endif
      k=ib-nb0+1
      if(k<1.or.k>eos%n) then;close(u);status='thermo density index outside eos.nb';return;endif
      do i=1,7
        call token_real(tok,3+i,qv(i),ok)
        if(.not.ok) then;close(u);status='invalid/nonfinite Q column';return;endif
      enddo
      call token_int(tok,11,nadd,ok)
      if(.not.ok.or.nadd<0.or.nadd>1) then;close(u);status='only Nadd=0 or Nadd=1 is supported';return;endif
      if(eos%nadd<0) then
        eos%nadd=nadd
      else if(eos%nadd/=nadd) then
        close(u);status='Nadd changes across eos.thermo rows';return
      endif
      if(nt/=11+nadd) then;close(u);status='Nadd does not match thermo row width';return;endif
      eos%q(:,k)=qv
      if(nadd==1) then
        call token_real(tok,12,v,ok)
        if(.not.ok) then;close(u);status='invalid q_extra';return;endif
        eos%q_extra(k)=v
      endif
      seen(k)=.true.
    enddo
    close(u)
    if(.not.all(seen)) then;status='eos.thermo does not cover every density grid point';return;endif
    if(eos%nadd<0)eos%nadd=0
    path=trim(directory)//'/eos.compo'
    if(allocated(seen)) deallocate(seen)
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    file_exists=ios==0
    eos%composition_file_present=file_exists
    if(eos%composition_file_present) then
      allocate(seen(eos%n));seen=.false.;eos%composition_available=.false.
      do
        call next_line(u,line,ios)
        if(ios<0) exit
        if(ios>0) then;close(u);status='error reading eos.compo';return;endif
        call split_tokens(line,tok,nt)
        if(nt<6) then;close(u);status='short eos.compo row';return;endif
        call token_int(tok,1,it,ok);if(.not.ok) then;close(u);status='invalid compo iT';return;endif
        call token_int(tok,2,ib,ok);if(.not.ok) then;close(u);status='invalid compo inb';return;endif
        call token_int(tok,3,iy,ok);if(.not.ok) then;close(u);status='invalid compo iYq';return;endif
        call token_int(tok,4,eos%iphase(max(1,min(eos%n,ib-nb0+1))),ok)
        if(.not.ok) then;close(u);status='invalid compo phase';return;endif
        k=ib-nb0+1
        if(k<1.or.k>eos%n.or.it/=t0.or.iy/=y0) then
          close(u);status='composition key outside selected thermo slice';return
        endif
        call token_int(tok,5,npairs,ok)
        if(.not.ok.or.npairs<0.or.npairs>16) then;close(u);status='invalid Npairs';return;endif
        pos=6
        do i=1,npairs
          if(pos+1>nt) then;close(u);status='truncated particle pair';return;endif
          call token_int(tok,pos,pid,ok);if(.not.ok) then;close(u);status='invalid particle id';return;endif
          call token_real(tok,pos+1,v,ok);if(.not.ok) then;close(u);status='invalid particle fraction';return;endif
          select case(pid)
          case(10);eos%yn(k)=v
          case(11);eos%yp(k)=v
          case(0);eos%ye(k)=v
          case(1);eos%ymu(k)=v
          case default
            close(u);status='unsupported CompOSE particle id in this EOSData model';return
          end select
          if(abs(v)>1.e-14_wp)eos%composition_available=.true.
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
          eos%abar(k)=a;eos%zbar(k)=z;eos%yh(k)=a*yi
          if(a>0._wp.and.yi>0._wp)eos%composition_available=.true.
          pos=pos+4
        enddo
        if(pos/=nt+1) then;close(u);status='trailing fields in eos.compo row';return;endif
        seen(k)=.true.
      enddo
      close(u)
      if(.not.all(seen)) then;status='eos.compo does not cover the thermo grid';return;endif
      if(any(eos%yn< -1.e-10_wp).or.any(eos%yp< -1.e-10_wp).or.any(eos%ye< -1.e-10_wp).or. &
         any(eos%ymu< -1.e-10_wp).or.any(eos%yh>1._wp+5.e-3_wp)) then
        status='composition fractions outside physical range';return
      endif
    else
      eos%composition_available=.false.
    endif
    eos%beta_equilibrated=maxval(abs(eos%q(5,:)))<1.e-8_wp
    status='ok'
  end subroutine read_compose

  subroutine write_internal_csv(path,eos,status)
    character(len=*),intent(in)::path
    type(eos_data_t),intent(in)::eos
    character(len=*),intent(out)::status
    integer::u,ios,i
    open(newunit=u,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0) then;status='cannot create internal EOS CSV';return;endif
    write(u,'(a)')'nb,Q1,Q2,Q3,Q4,Q5,Q6,Q7,q_extra,Iphase,Yn,Yp,Ye,Ymu,Abar,Zbar,Yh'
    do i=1,eos%n
      write(u,'(es25.16e3,7(",",es25.16e3),",",es25.16e3,",",i0,7(",",es25.16e3))') &
        eos%nb(i),eos%q(:,i),eos%q_extra(i),eos%iphase(i),eos%yn(i),eos%yp(i),eos%ye(i), &
        eos%ymu(i),eos%abar(i),eos%zbar(i),eos%yh(i)
    enddo
    close(u);status='ok'
  end subroutine write_internal_csv

  subroutine read_internal_csv(path,meta_path,eos,status)
    character(len=*),intent(in)::path,meta_path
    type(eos_data_t),intent(out)::eos
    character(len=*),intent(out)::status
    integer::u,ios,n,i,nt,k
    character(len=4096)::line
    character(len=96)::tok(64)
    logical::ok
    real(wp)::v
    integer::phase
    logical::composition_file_present,composition_available,beta_equilibrated
    integer::il,iT,iYq,inb_start,nadd
    real(wp)::mn,mp,temperature_mev,charge_fraction_grid
    namelist /metadata/ mn,mp,il,iT,iYq,inb_start,nadd,temperature_mev,charge_fraction_grid, &
       composition_file_present,composition_available,beta_equilibrated
    open(newunit=u,file=trim(meta_path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='cannot open internal EOS metadata';return;endif
    mn=939._wp;mp=939._wp;il=1;iT=1;iYq=1;inb_start=1;nadd=0
    temperature_mev=0._wp;charge_fraction_grid=0._wp
    composition_file_present=.false.;composition_available=.false.;beta_equilibrated=.false.
    read(u,nml=metadata,iostat=ios);close(u)
    if(ios/=0) then;status='invalid internal EOS metadata namelist';return;endif
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='cannot open internal EOS CSV';return;endif
    call next_line(u,line,ios)
    if(ios/=0) then;close(u);status='missing internal EOS CSV header';return;endif
    n=0
    do
      call next_line(u,line,ios)
      if(ios<0)exit
      if(ios>0)then;close(u);status='error counting internal EOS rows';return;endif
      n=n+1
    enddo
    close(u)
    call allocate_eos(eos,n)
    eos%mn=mn;eos%mp=mp;eos%il=il;eos%iT=iT;eos%iYq=iYq;eos%inb_start=inb_start;eos%nadd=nadd
    eos%temperature_mev=temperature_mev;eos%charge_fraction_grid=charge_fraction_grid
    eos%composition_file_present=composition_file_present;eos%composition_available=composition_available
    eos%beta_equilibrated=beta_equilibrated
    open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
    if(ios/=0) then;status='cannot reopen internal EOS CSV';return;endif
    call next_line(u,line,ios)
    do i=1,n
      call next_line(u,line,ios)
      if(ios/=0)then;close(u);status='unexpected EOF in internal EOS CSV';return;endif
      call split_tokens(line,tok,nt)
      if(nt/=17)then;close(u);status='internal EOS CSV must have 17 fields';return;endif
      do k=1,9
        if(k==10)cycle
        call token_real(tok,k,v,ok);if(.not.ok)then;close(u);status='invalid internal EOS value';return;endif
        select case(k)
        case(1);eos%nb(i)=v
        case(2:8);eos%q(k-1,i)=v
        case(9);eos%q_extra(i)=v
        end select
      enddo
      call token_int(tok,10,phase,ok);if(.not.ok)then;close(u);status='invalid Iphase';return;endif
      eos%iphase(i)=phase
      do k=11,17
        call token_real(tok,k,v,ok);if(.not.ok)then;close(u);status='invalid internal composition field';return;endif
        select case(k)
        case(11);eos%yn(i)=v
        case(12);eos%yp(i)=v
        case(13);eos%ye(i)=v
        case(14);eos%ymu(i)=v
        case(15);eos%abar(i)=v
        case(16);eos%zbar(i)=v
        case(17);eos%yh(i)=v
        end select
      enddo
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
    integer::u,ios,i,nquad
    character(len=1024)::path
    real(wp)::yi
    if(eos%n<1.or..not.allocated(eos%iphase)) then
      status='cannot write invalid EOS phase metadata'
      return
    endif
    if(size(eos%iphase)/=eos%n) then
      status='cannot write inconsistent EOS phase metadata'
      return
    endif
    if(eos%n>1) then
      if(any(eos%iphase(2:)/=eos%iphase(:eos%n-1)).and. &
         (.not.eos%composition_available.or..not.eos%composition_file_present)) then
        status='cannot export a phase-discontinuous barotrope as a complete CompOSE table without phase composition'
        return
      endif
    endif
    path=trim(directory)//'/eos.nb'
    open(newunit=u,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write eos.nb';return;endif
    write(u,'(i0)')eos%inb_start;write(u,'(i0)')eos%inb_start+eos%n-1
    do i=1,eos%n;write(u,'(es25.16e3)')eos%nb(i);enddo
    close(u)
    open(newunit=u,file=trim(directory)//'/eos.t',status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write eos.t';return;endif
    write(u,'(i0)')eos%iT;write(u,'(i0)')eos%iT;write(u,'(es25.16e3)')eos%temperature_mev;close(u)
    open(newunit=u,file=trim(directory)//'/eos.yq',status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write eos.yq';return;endif
    write(u,'(i0)')eos%iYq;write(u,'(i0)')eos%iYq;write(u,'(es25.16e3)')eos%charge_fraction_grid;close(u)
    open(newunit=u,file=trim(directory)//'/eos.thermo',status='replace',action='write',iostat=ios)
    if(ios/=0)then;status='cannot write eos.thermo';return;endif
    write(u,'(2(es25.16e3,1x),i0)')eos%mn,eos%mp,eos%il
    do i=1,eos%n
      if(eos%nadd==1)then
        write(u,'(3(i0,1x),7(es25.16e3,1x),i0,1x,es25.16e3)')eos%iT,eos%inb_start+i-1,eos%iYq,eos%q(:,i),1,eos%q_extra(i)
      else
        write(u,'(3(i0,1x),7(es25.16e3,1x),i0)')eos%iT,eos%inb_start+i-1,eos%iYq,eos%q(:,i),0
      endif
    enddo
    close(u)
    if(eos%composition_available)then
      open(newunit=u,file=trim(directory)//'/eos.compo',status='replace',action='write',iostat=ios)
      if(ios/=0)then;status='cannot write eos.compo';return;endif
      do i=1,eos%n
        nquad=0;if(eos%abar(i)>0._wp)nquad=1
        write(u,'(5(i0,1x),i0,1x,es25.16e3,1x,i0,1x,es25.16e3,1x,i0,1x,es25.16e3,1x,i0,1x,i0)',advance='no') &
          eos%iT,eos%inb_start+i-1,eos%iYq,eos%iphase(i),4,10,eos%yn(i),11,eos%yp(i),0,eos%ye(i),1
        ! Re-write each row in a simple internal formatted string to avoid format ambiguities.
        if(nquad==1)then
          yi=eos%yh(i)/eos%abar(i)
          write(u,'(es25.16e3,1x,i0,1x,i0,1x,es25.16e3,1x,es25.16e3,1x,es25.16e3)') &
            eos%ymu(i),1,1,eos%abar(i),eos%zbar(i),yi
        else
          write(u,'(es25.16e3,1x,i0)')eos%ymu(i),0
        endif
      enddo
      close(u)
    endif
    status='ok'
  end subroutine write_compose
end module compose_eos_io
