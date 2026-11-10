program r16_bridge_cell_audit
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    use ns_precision, only: wp
    use ns_status
    use ns_matching_types, only: mu_bridge_t
    use ns_mu_bridge, only: construct_bridge_endpoints
    implicit none

    character(len=1024) :: input_path
    character(len=4096) :: header
    integer :: arg_status, unit, ios, total, failures, reported
    real(wp) :: mu0, p0, n0, mu1, p1, n1
    real(wp) :: min_d2p, min_causal_margin
    type(mu_bridge_t) :: bridge
    type(status_t) :: status

    input_path = ''
    call get_command_argument(1, input_path, status=arg_status)
    if (arg_status /= 0 .or. len_trim(input_path) == 0) then
        write(error_unit,'(A)') 'usage: r16_bridge_cell_audit <si-cell-endpoints.csv>'
        error stop 2
    end if

    open(newunit=unit, file=trim(input_path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A)') 'could not open SI endpoint input file'
        error stop 2
    end if
    read(unit,'(A)',iostat=ios) header
    if (ios /= 0) then
        write(error_unit,'(A)') 'could not read SI endpoint CSV header'
        error stop 2
    end if

    total = 0
    failures = 0
    reported = 0
    min_d2p = huge(1._wp)
    min_causal_margin = huge(1._wp)

    do
        read(unit,*,iostat=ios) mu0, p0, n0, mu1, p1, n1
        if (ios < 0) exit
        if (ios > 0) then
            failures = failures + 1
            if (reported < 10) then
                write(error_unit,'(A,I0)') 'CSV parse failure after record ', total + 1
                reported = reported + 1
            end if
            cycle
        end if

        total = total + 1
        call construct_bridge_endpoints(mu0,p0,n0,mu1,p1,n1,bridge,status)
        if (.not.is_ok(status)) then
            failures = failures + 1
            if (reported < 10) then
                write(error_unit,'(A,I0,2A)') 'Fortran bridge rejected record ', total, ': ', trim(status%detail)
                reported = reported + 1
            end if
            cycle
        end if
        if (.not.bridge%ready .or. .not.bridge%interval_extrema_checked .or. &
            .not.bridge%stable_interval .or. .not.bridge%causal_interval) then
            failures = failures + 1
            if (reported < 10) then
                write(error_unit,'(A,I0)') 'Fortran bridge returned incomplete acceptance record ', total
                reported = reported + 1
            end if
            cycle
        end if

        min_d2p = min(min_d2p,bridge%min_d2p)
        min_causal_margin = min(min_causal_margin,bridge%min_causal_margin)
    end do
    close(unit)

    write(output_unit,'(A,I0)') 'r16_fortran_cells_total=',total
    write(output_unit,'(A,I0)') 'r16_fortran_cells_pass=',total-failures
    write(output_unit,'(A,I0)') 'r16_fortran_cells_fail=',failures
    write(output_unit,'(A,ES24.16E3)') 'r16_fortran_min_d2p=',min_d2p
    write(output_unit,'(A,ES24.16E3)') 'r16_fortran_min_causal_margin=',min_causal_margin

    if (total == 0 .or. failures /= 0) error stop 1
end program r16_bridge_cell_audit
