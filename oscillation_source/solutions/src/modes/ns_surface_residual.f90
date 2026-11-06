module ns_surface_residual
    use ns_precision, only: wp
    use ns_physical_constants, only: c_si
    use ns_status
    use ns_safe_math, only: safe_norm2
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: cowling_surface_residual
contains
    subroutine cowling_surface_residual(r,omega,nu,lambda,dnu,w,v,residual,q_over_r2,status)
        real(wp),intent(in)::r,omega,nu,lambda,dnu,w,v
        real(wp),intent(out)::residual,q_over_r2
        type(status_t),intent(out)::status
        real(wp)::z(2),b(2),nb,nz,k,expo
        type(status_t)::s1,s2
        residual=0._wp;q_over_r2=0._wp
        call fail_status(status,NS_INVALID_INPUT,'physics','surface residual requires finite r>0, omega>0 and nonzero vectors')
        if(.not.all(ieee_is_finite([r,omega,nu,lambda,dnu,w,v])).or.r<=0._wp.or.omega<=0._wp) return
        k=omega/c_si;expo=lambda-2._wp*nu
        if(.not.ieee_is_finite(expo).or.expo>log(huge(1._wp)).or.expo<log(tiny(1._wp))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','surface exponential unrepresentable');return
        end if
        z=[w/r**3,v/r**2];b=[r*dnu,(k*r)**2*exp(expo)]
        if(.not.all(ieee_is_finite(z)).or..not.all(ieee_is_finite(b))) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','surface vector nonfinite');return
        end if
        call safe_norm2(b,nb,s1);call safe_norm2(z,nz,s2)
        if(.not.is_ok(s1).or..not.is_ok(s2).or.nb<=0._wp.or.nz<=0._wp) then
            call fail_status(status,NS_INVALID_INPUT,'physics','surface residual normalization vector is zero/invalid');return
        end if
        residual=(b(1)/nb)*(z(1)/nz)+(b(2)/nb)*(z(2)/nz)
        q_over_r2=b(1)*z(1)+b(2)*z(2)
        if(.not.all(ieee_is_finite([residual,q_over_r2]))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','surface residual unrepresentable');return
        end if
        call clear_status(status)
    end subroutine
end module ns_surface_residual
