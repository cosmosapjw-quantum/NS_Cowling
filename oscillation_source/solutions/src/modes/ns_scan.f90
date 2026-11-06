module ns_scan
    use ns_precision, only: wp
    use ns_status
    use ns_scan_types
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: build_brackets,linear_scan_grid
contains
    subroutine linear_scan_grid(lo,hi,n,grid,status)
        real(wp),intent(in)::lo,hi
        integer,intent(in)::n
        real(wp),allocatable,intent(out)::grid(:)
        type(status_t),intent(out)::status
        integer::i,istat
        real(wp)::t
        call fail_status(status,NS_INVALID_INPUT,'numerics','finite increasing interval and 2..100000 samples required')
        if (.not.all(ieee_is_finite([lo,hi]))) return
        if (lo>=hi .or. n<2 .or. n>100000) return
        allocate(grid(n),stat=istat)
        if(istat/=0) then
            call fail_status(status,NS_RANGE_ERROR,'environment','scan allocation failed');return
        end if
        grid(1)=lo;grid(n)=hi
        do i=2,n-1
            t=real(i-1,wp)/real(n-1,wp)
            grid(i)=(1._wp-t)*lo+t*hi
        end do
        if (.not.all(ieee_is_finite(grid))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','unrepresentable grid');return
        end if
        do i=2,n
            if (grid(i)<=grid(i-1)) then
                call fail_status(status,NS_STAGNATION,'numerics','scan grid not distinguishable');return
            end if
        end do
        call clear_status(status)
    end subroutine

    subroutine build_brackets(trials,zero_tol,result)
        type(trial_record_t),intent(in)::trials(:)
        real(wp),intent(in)::zero_tol
        type(scan_result_t),intent(out)::result
        type(bracket_candidate_t),allocatable::work(:)
        type(bracket_candidate_t)::item
        logical,allocatable::usable(:),covered(:)
        integer::n,i,j,nc,orientation,istat
        result=scan_result_t();n=size(trials);nc=0;orientation=0
        call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid scan records or zero tolerance')
        if (n==0 .or. n>100000 .or. .not.ieee_is_finite(zero_tol)) return
        if (zero_tol<0._wp) return
        allocate(result%trials(n),usable(n),covered(n),work(2*n),stat=istat)
        if(istat/=0) then
            call fail_status(result%status,NS_RANGE_ERROR,'environment','scan allocation failed');return
        end if
        result%trials=trials;covered=.false.
        do i=1,n
            if (.not.ieee_is_finite(trials(i)%parameter)) return
            if (i>1) then
                if (trials(i)%parameter<=trials(i-1)%parameter) return
            end if
            ! A success flag alone cannot make a NaN residual an admissible shot.
            usable(i)=trials(i)%valid .and. is_ok(trials(i)%status) .and. ieee_is_finite(trials(i)%residual)
            if (.not.usable(i)) cycle
            if (abs(trials(i)%orientation)/=1) return
            if (orientation==0) orientation=trials(i)%orientation
            if (orientation/=trials(i)%orientation) then
                call fail_status(result%status,NS_INVALID_INPUT,'numerics','orientation changed within residual scan');return
            end if
        end do
        result%barriers=count(.not.usable);result%all_trials_succeeded=result%barriers==0
        do i=1,n-1
            if (.not.usable(i) .or. .not.usable(i+1)) cycle
            if (opposite(trials(i)%residual,trials(i+1)%residual)) then
                call add(i,i+1,CAND_SIGN)
                covered(i)=.true.;covered(i+1)=.true.
            end if
        end do
        do i=1,n
            if (.not.usable(i) .or. covered(i)) cycle
            if (abs(trials(i)%residual)>zero_tol) cycle
            ! Never jump over a failed or nonfinite neighbor. A small residual
            ! without a sign bracket remains an unverified sampled candidate.
            if (i>1 .and. i<n) then
                if (usable(i-1) .and. usable(i+1)) then
                    if (opposite(trials(i-1)%residual,trials(i+1)%residual)) then
                        call add(i-1,i+1,CAND_NEAR_BRACKET);cycle
                    end if
                end if
            end if
            call add(i,i,CAND_SAMPLE)
        end do
        ! Stable order. Adjacent sign-change cells are not merged: they may
        ! contain distinct roots. Covered near-zero endpoints were removed above.
        do i=2,nc
            item=work(i);j=i-1
            do while(j>=1)
                if(work(j)%lo<=item%lo) exit
                work(j+1)=work(j);j=j-1
            end do
            work(j+1)=item
        end do
        result%candidates=work(:nc)
        result%valid=.true.;call clear_status(result%status)
    contains
        subroutine add(left,right,kind)
            integer,intent(in)::left,right,kind
            nc=nc+1
            work(nc)=bracket_candidate_t(trials(left)%parameter,trials(right)%parameter,left,right,kind)
            if(kind==CAND_SAMPLE) then
                result%sample_candidates=result%sample_candidates+1
            else
                result%sign_brackets=result%sign_brackets+1
            end if
        end subroutine
    end subroutine

    pure logical function opposite(a,b)
        real(wp),intent(in)::a,b
        opposite=(a<0._wp .and. b>0._wp) .or. (a>0._wp .and. b<0._wp)
    end function
end module ns_scan
