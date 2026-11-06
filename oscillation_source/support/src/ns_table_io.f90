! T14 normalized cold-table reader.  It preserves source row order and performs no sorting,
! smoothing, deduplication, clipping, interpolation, or extrapolation.
module ns_table_io
    use ns_precision, only: wp
    use ns_status
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    character(len=*), parameter, public :: NORMALIZED_TABLE_HEADER = &
        'source_row,i_t,i_nb,i_yq,n_b_fm3,n_b_m3,pressure_mev_fm3,pressure_pa,'// &
        'energy_mev_fm3,energy_j_m3,mu_b_mev,mu_b_j,q1_p_over_nb_mev,'// &
        'q2_entropy_per_baryon,q3_mu_b_over_mn_minus_1,q4_mu_q_over_mn,'// &
        'q5_mu_l_over_mn,q6_f_over_nb_mn_minus_1,q7_e_over_nb_mn_minus_1'
    type, public :: normalized_table_t
        integer :: n=0
        real(wp), allocatable :: n_b_m3(:), pressure_pa(:), energy_j_m3(:), mu_j(:)
        logical :: ready=.false.
    end type normalized_table_t
    public :: clear_normalized_table, load_normalized_table
contains
    subroutine clear_normalized_table(table)
        type(normalized_table_t), intent(out) :: table
        table=normalized_table_t()
    end subroutine clear_normalized_table

    subroutine load_normalized_table(path,table,status)
        character(len=*), intent(in) :: path
        type(normalized_table_t), intent(out) :: table
        type(status_t), intent(out) :: status
        character(len=4096) :: line
        real(wp) :: v(19)
        integer :: u,ios,count_rows,i
        call clear_normalized_table(table)
        if(len_trim(path)==0) then
            call fail_status(status,NS_INVALID_INPUT,'data','empty normalized table path');return
        end if
        open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
        if(ios/=0) then
            call fail_status(status,NS_INVALID_INPUT,'data','cannot open normalized table');return
        end if
        read(u,'(A)',iostat=ios) line
        if(ios/=0 .or. trim(line)/=NORMALIZED_TABLE_HEADER) then
            close(u);call fail_status(status,NS_INVALID_INPUT,'data','unexpected normalized table header');return
        end if
        count_rows=0
        do
            read(u,'(A)',iostat=ios) line
            if(ios<0) exit
            if(ios>0) then
                close(u);call fail_status(status,NS_INVALID_INPUT,'data','normalized table read failure');return
            end if
            if(len_trim(line)==0) then
                close(u);call fail_status(status,NS_INVALID_INPUT,'data','blank normalized table row');return
            end if
            count_rows=count_rows+1
        end do
        close(u)
        if(count_rows<2) then
            call fail_status(status,NS_INVALID_INPUT,'data','normalized table requires at least two rows');return
        end if
        allocate(table%n_b_m3(count_rows),table%pressure_pa(count_rows), &
                 table%energy_j_m3(count_rows),table%mu_j(count_rows))
        open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
        if(ios/=0) then
            call clear_normalized_table(table)
            call fail_status(status,NS_INVALID_INPUT,'data','cannot reopen normalized table');return
        end if
        read(u,'(A)',iostat=ios) line
        do i=1,count_rows
            read(u,'(A)',iostat=ios) line
            if(ios/=0) then
                close(u);call clear_normalized_table(table)
                call fail_status(status,NS_INVALID_INPUT,'data','normalized table row count changed');return
            end if
            read(line,*,iostat=ios) v
            if(ios/=0) then
                close(u);call clear_normalized_table(table)
                call fail_status(status,NS_INVALID_INPUT,'data','malformed normalized table numeric row');return
            end if
            table%n_b_m3(i)=v(6);table%pressure_pa(i)=v(8)
            table%energy_j_m3(i)=v(10);table%mu_j(i)=v(12)
            if(.not.all(ieee_is_finite([table%n_b_m3(i),table%pressure_pa(i), &
                                        table%energy_j_m3(i),table%mu_j(i)]))) then
                close(u);call clear_normalized_table(table)
                call fail_status(status,NS_NONFINITE_INPUT,'data','nonfinite normalized table state');return
            end if
            if(min(table%n_b_m3(i),table%pressure_pa(i),table%energy_j_m3(i),table%mu_j(i))<=0._wp) then
                close(u);call clear_normalized_table(table)
                call fail_status(status,NS_INVALID_INPUT,'data','nonpositive normalized table state');return
            end if
        end do
        close(u)
        table%n=count_rows;table%ready=.true.;call clear_status(status)
    end subroutine load_normalized_table
end module ns_table_io
