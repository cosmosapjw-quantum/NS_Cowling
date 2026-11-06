! T15 controlled low-density perturbation of an admitted native-knot P(mu) table.
! The perturbation has compact support below split_index and vanishes with its
! first derivative at both endpoints.  Rows split_index:n are copied exactly.
module ns_controlled_table_variant
    use ns_precision, only: wp
    use ns_status
    use ns_table_io, only: normalized_table_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    character(len=*),parameter,public :: CONTROLLED_VARIANT_ID='compact_relative_Pmu_bump_v1'
    type,public :: controlled_variant_audit_t
        integer::split_index=0
        real(wp)::alpha=0._wp,split_n_b_m3=0._wp,split_mu_j=0._wp
        real(wp)::max_low_density_pressure_rel_change=0._wp
        real(wp)::max_low_density_density_rel_change=0._wp
        real(wp)::perturbation_endpoint_pressure_change_pa=0._wp
        real(wp)::perturbation_endpoint_density_change_m3=0._wp
        logical::core_function_unchanged=.false.,valid=.false.
    end type controlled_variant_audit_t
    public :: build_controlled_variant
contains
    subroutine build_controlled_variant(raw,alpha,split_index,variant,audit,status)
        type(normalized_table_t),intent(in)::raw
        real(wp),intent(in)::alpha
        integer,intent(in)::split_index
        type(normalized_table_t),intent(out)::variant
        type(controlled_variant_audit_t),intent(out)::audit
        type(status_t),intent(out)::status
        integer::i
        real(wp)::mu0,muspan,x,b,db,pnew,nnew,enew,pr,nr
        variant=normalized_table_t();audit=controlled_variant_audit_t()
        call fail_status(status,NS_INVALID_INPUT,'physics','valid table, amplitude and interior split index required')
        if(.not.raw%ready.or.raw%n<4.or..not.ieee_is_finite(alpha)) return
        if(split_index<3.or.split_index>=raw%n.or.abs(alpha)>0.5_wp) then
            call fail_status(status,NS_OUT_OF_DOMAIN,'physics','controlled variant split/amplitude outside teaching domain')
            return
        end if
        mu0=raw%mu_j(1);muspan=raw%mu_j(split_index)-mu0
        if(mu0<=0._wp.or.muspan<=0._wp) return
        variant=raw
        do i=1,split_index-1
            x=(raw%mu_j(i)-mu0)/muspan
            b=x*x*(1._wp-x)*(1._wp-x)
            db=2._wp*x*(1._wp-x)*(1._wp-2._wp*x)/muspan
            pnew=raw%pressure_pa(i)*(1._wp+alpha*b)
            nnew=raw%n_b_m3(i)*(1._wp+alpha*b)+alpha*raw%pressure_pa(i)*db
            enew=raw%mu_j(i)*nnew-pnew
            if(.not.all(ieee_is_finite([pnew,nnew,enew])).or.min(pnew,nnew,enew)<=0._wp) then
                variant=normalized_table_t()
                call fail_status(status,NS_MATCHING_INCONSISTENT,'physics','controlled low-density perturbation is not positive')
                return
            end if
            variant%pressure_pa(i)=pnew;variant%n_b_m3(i)=nnew;variant%energy_j_m3(i)=enew
            pr=abs(pnew/raw%pressure_pa(i)-1._wp);nr=abs(nnew/raw%n_b_m3(i)-1._wp)
            audit%max_low_density_pressure_rel_change=max(audit%max_low_density_pressure_rel_change,pr)
            audit%max_low_density_density_rel_change=max(audit%max_low_density_density_rel_change,nr)
        end do
        audit%split_index=split_index;audit%alpha=alpha
        audit%split_n_b_m3=raw%n_b_m3(split_index);audit%split_mu_j=raw%mu_j(split_index)
        audit%perturbation_endpoint_pressure_change_pa=variant%pressure_pa(split_index)-raw%pressure_pa(split_index)
        audit%perturbation_endpoint_density_change_m3=variant%n_b_m3(split_index)-raw%n_b_m3(split_index)
        audit%core_function_unchanged= &
            all(variant%pressure_pa(split_index:)==raw%pressure_pa(split_index:)).and. &
            all(variant%n_b_m3(split_index:)==raw%n_b_m3(split_index:)).and. &
            all(variant%energy_j_m3(split_index:)==raw%energy_j_m3(split_index:)).and. &
            all(variant%mu_j(split_index:)==raw%mu_j(split_index:))
        if(.not.audit%core_function_unchanged) then
            variant=normalized_table_t()
            call fail_status(status,NS_MATCHING_INCONSISTENT,'implementation','controlled variant changed frozen core rows')
            return
        end if
        variant%ready=.true.;audit%valid=.true.;call clear_status(status)
    end subroutine build_controlled_variant
end module ns_controlled_table_variant
