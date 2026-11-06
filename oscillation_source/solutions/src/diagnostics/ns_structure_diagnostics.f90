! T15 read-only diagnostics on an already accepted static background.
! The central 0<r<r0 diagnostic segment holds the first accepted row fixed;
! this does not modify the TOV solution or create a new stellar background.
module ns_structure_diagnostics
    use ns_precision, only: wp
    use ns_physical_constants, only: c_si,g_si,two_pi
    use ns_status
    use ns_background_access, only: background_table_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    character(len=*),parameter,public :: STRUCTURE_DIAGNOSTIC_ID='background_trapezoid_center_hold_v1'
    type,public :: structure_diagnostics_t
        real(wp) :: radius_m=0._wp,mass_kg=0._wp,proper_radius_m=0._wp,acoustic_time_s=0._wp
        real(wp) :: r_split_m=0._wp,outer_areal_thickness_m=0._wp,outer_proper_thickness_m=0._wp
        real(wp) :: outer_acoustic_time_s=0._wp,sound_speed_min_m_s=0._wp,sound_speed_max_m_s=0._wp
        real(wp) :: f_hat=0._wp,p1_hat=0._wp
        logical :: valid=.false.
    end type structure_diagnostics_t
    public :: compute_structure_diagnostics
contains
    subroutine compute_structure_diagnostics(bg,split_pressure_pa,f_hz,p1_hz,result,status)
        type(background_table_t),intent(in)::bg
        real(wp),intent(in)::split_pressure_pa,f_hz,p1_hz
        type(structure_diagnostics_t),intent(out)::result
        type(status_t),intent(out)::status
        integer::i,k
        real(wp)::dr,fp0,fp1,fa0,fa1,a,rs,fps,fas,dyn
        logical::found
        result=structure_diagnostics_t()
        call fail_status(status,NS_INVALID_INPUT,'data','valid accepted background and diagnostic inputs required')
        if(.not.bg%ready.or.bg%n<2) return
        if(.not.all(ieee_is_finite([split_pressure_pa,f_hz,p1_hz]))) return
        if(split_pressure_pa<=0._wp.or.f_hz<=0._wp.or.p1_hz<=0._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','positive pressure and frequencies required');return
        end if
        if(split_pressure_pa>bg%p_pa(1).or.split_pressure_pa<bg%p_pa(bg%n)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','split pressure outside stored background');return
        end if
        if(any(bg%chi<=0._wp).or.any(bg%chi>1._wp)) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','background sound-speed response outside admitted range');return
        end if
        result%radius_m=bg%r_m(bg%n);result%mass_kg=bg%m_kg(bg%n)
        if(result%radius_m<=0._wp.or.result%mass_kg<=0._wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','positive terminal radius and mass required');return
        end if

        ! Diagnostic closure for the omitted center series interval 0<r<r0.
        fp0=exp(bg%lambda(1));fa0=exp(bg%lambda(1)-bg%nu(1))/(c_si*sqrt(bg%chi(1)))
        result%proper_radius_m=bg%r_m(1)*fp0
        result%acoustic_time_s=bg%r_m(1)*fa0
        do i=1,bg%n-1
            dr=bg%r_m(i+1)-bg%r_m(i)
            fp0=exp(bg%lambda(i));fp1=exp(bg%lambda(i+1))
            fa0=exp(bg%lambda(i)-bg%nu(i))/(c_si*sqrt(bg%chi(i)))
            fa1=exp(bg%lambda(i+1)-bg%nu(i+1))/(c_si*sqrt(bg%chi(i+1)))
            result%proper_radius_m=result%proper_radius_m+0.5_wp*dr*(fp0+fp1)
            result%acoustic_time_s=result%acoustic_time_s+0.5_wp*dr*(fa0+fa1)
        end do

        found=.false.;k=0;a=0._wp
        if(split_pressure_pa==bg%p_pa(bg%n)) then
            result%r_split_m=bg%r_m(bg%n);found=.true.;k=bg%n
        else
            do i=1,bg%n-1
                if(bg%p_pa(i)>=split_pressure_pa.and.split_pressure_pa>=bg%p_pa(i+1)) then
                    a=(bg%p_pa(i)-split_pressure_pa)/(bg%p_pa(i)-bg%p_pa(i+1))
                    result%r_split_m=(1._wp-a)*bg%r_m(i)+a*bg%r_m(i+1)
                    found=.true.;k=i;exit
                end if
            end do
        end if
        if(.not.found) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'numerics','split pressure could not be located');return
        end if
        result%outer_areal_thickness_m=result%radius_m-result%r_split_m
        result%outer_proper_thickness_m=0._wp;result%outer_acoustic_time_s=0._wp
        if(k<bg%n) then
            rs=result%r_split_m
            fp0=exp(bg%lambda(k));fp1=exp(bg%lambda(k+1));fps=(1._wp-a)*fp0+a*fp1
            fa0=exp(bg%lambda(k)-bg%nu(k))/(c_si*sqrt(bg%chi(k)))
            fa1=exp(bg%lambda(k+1)-bg%nu(k+1))/(c_si*sqrt(bg%chi(k+1)));fas=(1._wp-a)*fa0+a*fa1
            dr=bg%r_m(k+1)-rs
            result%outer_proper_thickness_m=0.5_wp*dr*(fps+fp1)
            result%outer_acoustic_time_s=0.5_wp*dr*(fas+fa1)
            do i=k+1,bg%n-1
                dr=bg%r_m(i+1)-bg%r_m(i)
                fp0=exp(bg%lambda(i));fp1=exp(bg%lambda(i+1))
                fa0=exp(bg%lambda(i)-bg%nu(i))/(c_si*sqrt(bg%chi(i)))
                fa1=exp(bg%lambda(i+1)-bg%nu(i+1))/(c_si*sqrt(bg%chi(i+1)))
                result%outer_proper_thickness_m=result%outer_proper_thickness_m+0.5_wp*dr*(fp0+fp1)
                result%outer_acoustic_time_s=result%outer_acoustic_time_s+0.5_wp*dr*(fa0+fa1)
            end do
        end if
        result%sound_speed_min_m_s=c_si*sqrt(minval(bg%chi))
        result%sound_speed_max_m_s=c_si*sqrt(maxval(bg%chi))
        dyn=sqrt(g_si*result%mass_kg/result%radius_m**3)
        result%f_hat=two_pi*f_hz/dyn;result%p1_hat=two_pi*p1_hz/dyn
        if(.not.all(ieee_is_finite([result%proper_radius_m,result%acoustic_time_s,result%r_split_m, &
            result%outer_proper_thickness_m,result%outer_acoustic_time_s,result%f_hat,result%p1_hat]))) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','nonfinite structure diagnostic');return
        end if
        result%valid=.true.;call clear_status(status)
    end subroutine compute_structure_diagnostics
end module ns_structure_diagnostics
