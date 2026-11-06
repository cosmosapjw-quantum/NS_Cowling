! Dimensionless profile diagnostics; this norm is NOT physical mode energy.
module ns_mode_normalization
    use ns_precision, only: wp
    use ns_status
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    character(len=*),parameter,public :: MODE_NORM_ID='scaled_WV_L2_dx_v1'
    public :: normalize_mode, count_smooth_nodes
contains
    subroutine normalize_mode(r,w,v,x,wn,vn,status)
        real(wp),intent(in)::r(:),w(:),v(:)
        real(wp),allocatable,intent(out)::x(:),wn(:),vn(:)
        type(status_t),intent(out)::status
        real(wp)::radius,scale,integral
        integer::n,i,m,j
        integer,allocatable::keep(:)
        call fail_status(status,NS_INVALID_INPUT,'numerics','finite increasing r and nonzero smooth mode required')
        n=size(r)
        if(n<3.or.size(w)/=n.or.size(v)/=n) return
        if(.not.all(ieee_is_finite(r)).or..not.all(ieee_is_finite(w)).or..not.all(ieee_is_finite(v))) return
        if(any(r<=0._wp)) return
        allocate(keep(n));m=1;keep(1)=1
        do i=2,n
            if(r(i)<r(i-1)) return
            if(r(i)==r(i-1)) then
                ! Only collapse redundant rounded coordinates, never a jump.
                if(abs(w(i)-w(i-1))>64._wp*epsilon(1._wp)*max(abs(w(i)),abs(w(i-1)))) return
                if(abs(v(i)-v(i-1))>64._wp*epsilon(1._wp)*max(abs(v(i)),abs(v(i-1)))) return
                cycle
            end if
            m=m+1;keep(m)=i
        end do
        if(m<3) return
        radius=r(n);x=r(keep(:m))/radius
        ! Divide successively instead of cubing a dimensional scale.
        wn=((w(keep(:m))/radius)/radius)/radius;vn=(v(keep(:m))/radius)/radius
        n=m
        scale=max(maxval(abs(wn)),maxval(abs(vn)))
        if(.not.ieee_is_finite(scale).or.scale<=0._wp) return
        wn=wn/scale;vn=vn/scale
        integral=0._wp
        do i=2,n
            integral=integral+0.5_wp*(x(i)-x(i-1))* &
                (wn(i)**2+wn(i-1)**2+6._wp*(vn(i)**2+vn(i-1)**2))
        end do
        if(.not.ieee_is_finite(integral).or.integral<=0._wp) return
        wn=wn/sqrt(integral);vn=vn/sqrt(integral)
        call clear_status(status)
    end subroutine

    subroutine count_smooth_nodes(x,y,xlo,xhi,relative_zero,nodes,status)
        real(wp),intent(in)::x(:),y(:),xlo,xhi,relative_zero
        integer,intent(out)::nodes
        type(status_t),intent(out)::status
        integer::i,n,last_sign,s,kept
        real(wp)::threshold,scale
        nodes=-1;n=size(x)
        call fail_status(status,NS_INVALID_INPUT,'numerics','node count requires ordered smooth samples and a nonzero interior')
        if(n<3.or.size(y)/=n) return
        if(.not.all(ieee_is_finite(x)).or..not.all(ieee_is_finite(y))) return
        if(.not.all(ieee_is_finite([xlo,xhi,relative_zero]))) return
        if(xlo<0._wp.or.xhi>1._wp.or.xlo>=xhi.or.relative_zero<0._wp.or.relative_zero>=1._wp) return
        do i=2,n
            if(x(i)<=x(i-1)) return
        end do
        scale=0._wp;kept=0
        do i=1,n
            if(x(i)<xlo.or.x(i)>xhi) cycle
            kept=kept+1;scale=max(scale,abs(y(i)))
        end do
        if(kept<3.or.scale<=0._wp) return
        threshold=relative_zero*scale;nodes=0;last_sign=0;kept=0
        do i=1,n
            if(x(i)<xlo.or.x(i)>xhi) cycle
            if(abs(y(i))<=threshold) cycle
            s=1;if(y(i)<0._wp) s=-1
            if(last_sign/=0.and.s/=last_sign) nodes=nodes+1
            last_sign=s;kept=kept+1
        end do
        if(kept<3) then
            nodes=-1;return
        end if
        call clear_status(status)
    end subroutine
end module ns_mode_normalization
