module ns_background_access
    use ns_precision, only: wp
    use ns_status
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    character(len=*), parameter, public :: BACKGROUND_INTERPOLATION_ID='linear_q_v1'
    type, public :: background_point_t
        real(wp) :: q=0._wp,r_m=0._wp,m_kg=0._wp,p_pa=0._wp,energy_j_m3=0._wp,n_b_m3=0._wp,mu_j=0._wp
        real(wp) :: chi=0._wp,nu=0._wp,lambda=0._wp,dnu_dr=0._wp
        logical :: valid=.false.
    end type
    type, public :: background_table_t
        real(wp),allocatable :: q(:),r_m(:),m_kg(:),p_pa(:),energy_j_m3(:),n_b_m3(:),mu_j(:),chi(:),nu(:),lambda(:),dnu_dr(:)
        integer :: n=0
        logical :: ready=.false.
    end type
    public :: load_background_csv,background_at_r,background_at_q,background_bounds
contains
    subroutine clear_table(t)
        type(background_table_t),intent(out)::t
        t=background_table_t()
    end subroutine
    subroutine load_background_csv(path,t,status)
        character(len=*),intent(in)::path
        type(background_table_t),intent(out)::t
        type(status_t),intent(out)::status
        character(len=1024)::line,header
        integer::u,ios,n,i,astat
        real(wp)::vals(11)
        call clear_table(t)
        call fail_status(status,NS_INVALID_INPUT,'data','valid fixed-background CSV required')
        if(len_trim(path)==0) return
        open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
        if(ios/=0) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'data','background file unavailable');return
        end if
        read(u,'(A)',iostat=ios) header
        if(ios/=0 .or. trim(header)/='q,r_m,m_kg,p_pa,energy_j_m3,n_b_m3,mu_j,chi,nu,lambda,dnu_dr') then
            close(u);call fail_status(status,NS_INVALID_INPUT,'data','unexpected background columns');return
        end if
        n=0
        do
            read(u,'(A)',iostat=ios) line
            if(ios<0) exit
            if(ios>0) then
                close(u);call fail_status(status,NS_INVALID_INPUT,'data','background read failure');return
            end if
            if(len_trim(line)>0) n=n+1
        end do
        if(n<2 .or. n>100000) then
            close(u);call fail_status(status,NS_INVALID_INPUT,'data','background row count outside supported range');return
        end if
        rewind(u);read(u,'(A)') header
        allocate(t%q(n),t%r_m(n),t%m_kg(n),t%p_pa(n),t%energy_j_m3(n), &
                 t%n_b_m3(n),t%mu_j(n),t%chi(n),t%nu(n),t%lambda(n), &
                 t%dnu_dr(n),stat=astat)
        if(astat/=0) then
            close(u);call fail_status(status,NS_RANGE_ERROR,'environment','background allocation failed');return
        end if
        do i=1,n
            read(u,'(A)',iostat=ios) line
            if(ios/=0) then
                close(u);call fail_status(status,NS_INVALID_INPUT,'data','background row missing');return
            end if
            read(line,*,iostat=ios) vals
            if(ios/=0 .or. .not.all(ieee_is_finite(vals))) then
                close(u);call fail_status(status,NS_INVALID_INPUT,'data','malformed/nonfinite background row');return
            end if
            t%q(i)=vals(1);t%r_m(i)=vals(2);t%m_kg(i)=vals(3);t%p_pa(i)=vals(4);t%energy_j_m3(i)=vals(5)
            t%n_b_m3(i)=vals(6);t%mu_j(i)=vals(7);t%chi(i)=vals(8);t%nu(i)=vals(9);t%lambda(i)=vals(10);t%dnu_dr(i)=vals(11)
            if(t%r_m(i)<=0._wp.or.t%m_kg(i)<0._wp.or.t%p_pa(i)<=0._wp.or.t%energy_j_m3(i)<=0._wp.or. &
               t%n_b_m3(i)<=0._wp.or.t%mu_j(i)<=0._wp.or.t%chi(i)<=0._wp.or.t%chi(i)>1._wp.or.t%dnu_dr(i)<=0._wp) then
                close(u);call fail_status(status,NS_OUT_OF_DOMAIN,'physics','background violates smooth barotropic contract');return
            end if
            if(i>1) then
                if(t%q(i)<=t%q(i-1).or.t%r_m(i)<=t%r_m(i-1).or.t%m_kg(i)<t%m_kg(i-1).or.t%p_pa(i)>=t%p_pa(i-1)) then
                    close(u);call fail_status(status,NS_INVALID_INPUT,'data','background order is not monotone');return
                end if
            end if
        end do
        close(u);t%n=n;t%ready=.true.;call clear_status(status)
    end subroutine
    subroutine background_bounds(t,q0,q1,r0,r1,status)
        type(background_table_t),intent(in)::t
        real(wp),intent(out)::q0,q1,r0,r1
        type(status_t),intent(out)::status
        q0=0;q1=0;r0=0;r1=0
        call fail_status(status,NS_INVALID_INPUT,'data','background not ready')
        if(.not.t%ready.or.t%n<2) return
        q0=t%q(1);q1=t%q(t%n);r0=t%r_m(1);r1=t%r_m(t%n);call clear_status(status)
    end subroutine
    subroutine background_at_q(t,q,point,status,nu_shift)
        type(background_table_t),intent(in)::t
        real(wp),intent(in)::q
        type(background_point_t),intent(out)::point
        type(status_t),intent(out)::status
        real(wp),intent(in),optional::nu_shift
        integer::lo,hi,mid
        real(wp)::a,c
        point=background_point_t();call fail_status(status,NS_INVALID_INPUT,'data','invalid background query')
        if(.not.t%ready.or..not.ieee_is_finite(q)) return
        if(q<t%q(1).or.q>t%q(t%n)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','background q extrapolation forbidden');return
        end if
        if(q==t%q(1)) then;lo=1;hi=1
        else if(q==t%q(t%n)) then;lo=t%n;hi=t%n
        else
            lo=1;hi=t%n
            do while(hi-lo>1)
                mid=(lo+hi)/2
                if(t%q(mid)<=q) then;lo=mid;else;hi=mid;end if
            end do
        end if
        if(lo==hi) then;a=0._wp
        else;a=(q-t%q(lo))/(t%q(hi)-t%q(lo))
        end if
        c=0._wp;if(present(nu_shift)) c=nu_shift
        point%q=q
        point%r_m=lerp(t%r_m,lo,hi,a);point%m_kg=lerp(t%m_kg,lo,hi,a);point%p_pa=lerp(t%p_pa,lo,hi,a)
        point%energy_j_m3=lerp(t%energy_j_m3,lo,hi,a);point%n_b_m3=lerp(t%n_b_m3,lo,hi,a);point%mu_j=lerp(t%mu_j,lo,hi,a)
        point%chi=lerp(t%chi,lo,hi,a);point%nu=lerp(t%nu,lo,hi,a)+c
        point%lambda=lerp(t%lambda,lo,hi,a);point%dnu_dr=lerp(t%dnu_dr,lo,hi,a)
        if(.not.all(ieee_is_finite([point%r_m,point%m_kg,point%p_pa,point%energy_j_m3, &
             point%n_b_m3,point%mu_j,point%chi,point%nu,point%lambda,point%dnu_dr]))) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','background interpolation nonfinite');return
        end if
        if(point%r_m<=0._wp.or.point%chi<=0._wp.or.point%chi>1._wp.or.point%dnu_dr<=0._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','interpolated background not admissible');return
        end if
        point%valid=.true.;call clear_status(status)
    end subroutine
    subroutine background_at_r(t,r,point,status,nu_shift)
        type(background_table_t),intent(in)::t
        real(wp),intent(in)::r
        type(background_point_t),intent(out)::point
        type(status_t),intent(out)::status
        real(wp),intent(in),optional::nu_shift
        integer::lo,hi,mid
        real(wp)::a,q
        point=background_point_t();call fail_status(status,NS_INVALID_INPUT,'data','invalid background query')
        if(.not.t%ready.or..not.ieee_is_finite(r)) return
        if(r<t%r_m(1).or.r>t%r_m(t%n)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','background r extrapolation forbidden');return
        end if
        if(r==t%r_m(1)) then;q=t%q(1)
        else if(r==t%r_m(t%n)) then;q=t%q(t%n)
        else
            lo=1;hi=t%n
            do while(hi-lo>1)
                mid=(lo+hi)/2
                if(t%r_m(mid)<=r) then;lo=mid;else;hi=mid;end if
            end do
            a=(r-t%r_m(lo))/(t%r_m(hi)-t%r_m(lo));q=(1._wp-a)*t%q(lo)+a*t%q(hi)
        end if
        if(present(nu_shift)) then
            call background_at_q(t,q,point,status,nu_shift)
        else
            call background_at_q(t,q,point,status)
        end if
        if(is_ok(status)) point%r_m=r
    end subroutine
    pure real(wp) function lerp(v,lo,hi,a)
        real(wp),intent(in)::v(:),a
        integer,intent(in)::lo,hi
        if(lo==hi) then;lerp=v(lo);else;lerp=(1._wp-a)*v(lo)+a*v(hi);end if
    end function
end module ns_background_access
