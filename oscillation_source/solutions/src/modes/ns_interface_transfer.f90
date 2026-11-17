module ns_interface_transfer
    use ns_precision, only: wp
    use ns_status
    use ns_physical_constants, only: c_si
    use ns_eos_types, only: SIDE_INNER,SIDE_OUTER
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: material_interface_transfer,rapid_interface_transfer,lagrangian_pressure_jump,identify_interface_family
contains
    subroutine lagrangian_pressure_jump(r,omega,p,eps,nu,lambda,dnu,w,v,delta_p,qterm,status)
        real(wp),intent(in)::r,omega,p,eps,nu,lambda,dnu,w,v
        real(wp),intent(out)::delta_p,qterm
        type(status_t),intent(out)::status
        real(wp)::a,expo,h
        delta_p=0._wp;qterm=0._wp
        call fail_status(status,NS_INVALID_INPUT,'physics', &
            'material-interface diagnostic requires positive radius/frequency/enthalpy')
        if(.not.all(ieee_is_finite([r,omega,p,eps,nu,lambda,dnu,w,v]))) return
        if(r<=0._wp.or.omega<=0._wp.or.p<0._wp.or.eps+p<=0._wp) return
        expo=lambda-2._wp*nu
        if(expo>log(huge(1._wp)).or.expo<log(tiny(1._wp))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','interface exponential unrepresentable');return
        end if
        a=(omega/c_si*r)**2*exp(expo);h=eps+p
        qterm=a*v+dnu*w
        delta_p=-h*exp(-lambda)*qterm/(r*r)
        if(.not.all(ieee_is_finite([a,qterm,delta_p])).or.a<=0._wp) then
            call fail_status(status,NS_NONFINITE_RESULT,'numerics','interface diagnostic nonfinite');return
        end if
        call clear_status(status)
    end subroutine

    subroutine material_interface_transfer(r,omega,p,eps_in,eps_out,nu,lambda,dnu, &
        w_in,v_in,inner_side,outer_side,w_out,v_out,dp_in,dp_out,status)
        real(wp),intent(in)::r,omega,p,eps_in,eps_out,nu,lambda,dnu,w_in,v_in
        integer,intent(in)::inner_side,outer_side
        real(wp),intent(out)::w_out,v_out,dp_in,dp_out
        type(status_t),intent(out)::status
        real(wp)::q_in,q_out,a,h_in,h_out,target,expo
        w_out=0._wp;v_out=0._wp;dp_in=0._wp;dp_out=0._wp
        call fail_status(status,NS_INVALID_INPUT,'physics', &
            'slow material interface requires inner-to-outer sides and omega>0')
        if(inner_side/=SIDE_INNER.or.outer_side/=SIDE_OUTER) return
        call lagrangian_pressure_jump(r,omega,p,eps_in,nu,lambda,dnu,w_in,v_in, &
            dp_in,q_in,status);if(.not.is_ok(status))return
        h_in=eps_in+p;h_out=eps_out+p
        if(h_out<=0._wp) return
        expo=lambda-2._wp*nu;a=(omega/c_si*r)**2*exp(expo)
        if(a<=0._wp.or..not.ieee_is_finite(a)) return
        target=(h_in/h_out)*q_in
        w_out=w_in;v_out=(target-dnu*w_out)/a
        call lagrangian_pressure_jump(r,omega,p,eps_out,nu,lambda,dnu,w_out,v_out, &
            dp_out,q_out,status);if(.not.is_ok(status))return
        if(abs(dp_out-dp_in)>2e-12_wp*max(abs(dp_in),abs(dp_out),1._wp)) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'Delta P not conserved by material-interface transfer');return
        end if
        call clear_status(status)
    end subroutine


    subroutine rapid_interface_transfer(r,omega,p,eps_in,eps_out,nu,lambda,dnu, &
        w_in,v_in,inner_side,outer_side,w_out,v_out,dp_in,dp_out,status)
        real(wp),intent(in)::r,omega,p,eps_in,eps_out,nu,lambda,dnu,w_in,v_in
        integer,intent(in)::inner_side,outer_side
        real(wp),intent(out)::w_out,v_out,dp_in,dp_out
        type(status_t),intent(out)::status
        real(wp)::q_in,q_out,q_check,a,h_in,h_out,expo,xi_in,xi_out,iso_in,iso_out,pprime_in,pprime_out
        w_out=0._wp;v_out=0._wp;dp_in=0._wp;dp_out=0._wp
        call fail_status(status,NS_INVALID_INPUT,'physics', &
            'rapid-conversion interface requires inner-to-outer sides and omega>0')
        if(inner_side/=SIDE_INNER.or.outer_side/=SIDE_OUTER) return
        if(.not.ieee_is_finite(dnu).or.dnu<=0._wp) return
        call lagrangian_pressure_jump(r,omega,p,eps_in,nu,lambda,dnu,w_in,v_in, &
            dp_in,q_in,status);if(.not.is_ok(status))return
        h_in=eps_in+p;h_out=eps_out+p
        if(h_in<=0._wp.or.h_out<=0._wp) return
        expo=lambda-2._wp*nu
        if(expo>log(huge(1._wp)).or.expo<log(tiny(1._wp))) then
            call fail_status(status,NS_RANGE_ERROR,'numerics','rapid-interface exponential unrepresentable');return
        end if
        a=(omega/c_si*r)**2*exp(expo)
        if(a<=0._wp.or..not.ieee_is_finite(a)) return
        q_out=(h_in/h_out)*q_in
        ! [Delta P]=0 and [xi^r-Delta P/P'_0]=0.
        ! With P'_0=-(eps+P)*nu', the latter is [W-Q/nu']=0.
        w_out=w_in+(q_out-q_in)/dnu
        v_out=(q_out-dnu*w_out)/a
        call lagrangian_pressure_jump(r,omega,p,eps_out,nu,lambda,dnu,w_out,v_out, &
            dp_out,q_check,status);if(.not.is_ok(status))return
        if(abs(dp_out-dp_in)>2e-12_wp*max(abs(dp_in),abs(dp_out),1._wp)) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'Delta P not conserved by rapid-interface transfer');return
        end if
        xi_in=exp(-lambda)*w_in/(r*r);xi_out=exp(-lambda)*w_out/(r*r)
        pprime_in=-h_in*dnu;pprime_out=-h_out*dnu
        iso_in=xi_in-dp_in/pprime_in;iso_out=xi_out-dp_out/pprime_out
        if(abs(iso_out-iso_in)>2e-12_wp*max(abs(iso_in),abs(iso_out),1._wp)) then
            call fail_status(status,NS_MATCHING_INCONSISTENT,'physics', &
                'rapid isobar-displacement condition not conserved');return
        end if
        call clear_status(status)
    end subroutine

    subroutine identify_interface_family(reference_label,root_valid,transfer_valid,overlap, &
        candidate_interface_g,label,status)
        character(len=*),intent(in)::reference_label
        logical,intent(in)::root_valid,transfer_valid,candidate_interface_g
        real(wp),intent(in)::overlap
        character(len=*),intent(out)::label
        type(status_t),intent(out)::status
        label='unresolved'
        call fail_status(status,NS_IDENTITY_UNRESOLVED,'physics', &
            'interface family requires smooth reference overlap and valid transfer')
        if(.not.root_valid.or..not.transfer_valid.or.candidate_interface_g) return
        if(.not.ieee_is_finite(overlap).or.overlap<0.8_wp.or.overlap>1.0000001_wp) return
        if(trim(reference_label)/='f'.and.trim(reference_label)/='p1') return
        label=trim(reference_label);call clear_status(status)
    end subroutine
end module ns_interface_transfer
