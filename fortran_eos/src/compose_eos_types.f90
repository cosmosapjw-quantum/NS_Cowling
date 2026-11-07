module compose_eos_types
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  integer, parameter, public :: wp=real64
  type, public :: eos_data_t
    integer :: n=0, il=1, iT=1, iYq=1, inb_start=1, nadd=0
    real(wp) :: mn=939._wp, mp=939._wp, temperature_mev=0._wp, charge_fraction_grid=0._wp
    logical :: composition_file_present=.false., composition_available=.false., beta_equilibrated=.false.
    real(wp), allocatable :: nb(:), q(:,:), q_extra(:)
    integer, allocatable :: iphase(:)
    real(wp), allocatable :: yn(:), yp(:), ye(:), ymu(:), abar(:), zbar(:), yh(:)
  end type eos_data_t
  public :: allocate_eos, validate_eos, eos_pressure, eos_energy, eos_mu_b, eos_mu_q, eos_mu_l, eos_free_energy
contains
  subroutine allocate_eos(eos,n)
    type(eos_data_t),intent(inout)::eos
    integer,intent(in)::n
    if(allocated(eos%nb)) deallocate(eos%nb)
    if(allocated(eos%q)) deallocate(eos%q)
    if(allocated(eos%q_extra)) deallocate(eos%q_extra)
    if(allocated(eos%iphase)) deallocate(eos%iphase)
    if(allocated(eos%yn)) deallocate(eos%yn)
    if(allocated(eos%yp)) deallocate(eos%yp)
    if(allocated(eos%ye)) deallocate(eos%ye)
    if(allocated(eos%ymu)) deallocate(eos%ymu)
    if(allocated(eos%abar)) deallocate(eos%abar)
    if(allocated(eos%zbar)) deallocate(eos%zbar)
    if(allocated(eos%yh)) deallocate(eos%yh)
    eos%n=n
    allocate(eos%nb(n),eos%q(7,n),eos%q_extra(n),eos%iphase(n))
    allocate(eos%yn(n),eos%yp(n),eos%ye(n),eos%ymu(n),eos%abar(n),eos%zbar(n),eos%yh(n))
    eos%nb=0._wp; eos%q=0._wp; eos%q_extra=0._wp; eos%iphase=0
    eos%yn=0._wp; eos%yp=0._wp; eos%ye=0._wp; eos%ymu=0._wp
    eos%abar=0._wp; eos%zbar=0._wp; eos%yh=0._wp
  end subroutine allocate_eos

  subroutine validate_eos(eos,ok,message)
    type(eos_data_t),intent(in)::eos
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    integer::i
    ok=.false.; message='EOS arrays are not allocated'
    if(eos%n<3) then
      message='EOS requires at least three density points';return
    endif
    if(.not.allocated(eos%nb).or..not.allocated(eos%q)) return
    if(size(eos%nb)/=eos%n.or.any(shape(eos%q)/=[7,eos%n])) then
      message='EOS array shapes are inconsistent';return
    endif
    if(.not.all(ieee_finite(eos%nb)).or..not.all(ieee_finite(eos%q))) then
      message='EOS primary arrays contain nonfinite values';return
    endif
    if(any(eos%nb<=0._wp).or.any(eos%nb(2:)<=eos%nb(:eos%n-1))) then
      message='EOS density grid must be positive and strictly increasing';return
    endif
    if(eos%mn<=0._wp.or.eos%mp<=0._wp.or..not.ieee_finite(eos%mn+eos%mp)) then
      message='invalid nucleon masses';return
    endif
    do i=1,eos%n
      if(eos%q(1,i)*eos%nb(i)<=0._wp) then
        message='nonpositive pressure encountered';return
      endif
    enddo
    ok=.true.;message='ok'
  end subroutine validate_eos

  pure function eos_pressure(eos) result(p)
    type(eos_data_t),intent(in)::eos
    real(wp)::p(eos%n)
    p=eos%q(1,:)*eos%nb
  end function
  pure function eos_mu_b(eos) result(mu)
    type(eos_data_t),intent(in)::eos
    real(wp)::mu(eos%n)
    mu=(eos%q(3,:)+1._wp)*eos%mn
  end function
  pure function eos_mu_q(eos) result(mu)
    type(eos_data_t),intent(in)::eos
    real(wp)::mu(eos%n)
    mu=eos%q(4,:)*eos%mn
  end function
  pure function eos_mu_l(eos) result(mu)
    type(eos_data_t),intent(in)::eos
    real(wp)::mu(eos%n)
    mu=eos%q(5,:)*eos%mn
  end function
  pure function eos_free_energy(eos) result(f)
    type(eos_data_t),intent(in)::eos
    real(wp)::f(eos%n)
    f=(eos%q(6,:)+1._wp)*eos%mn
  end function
  pure function eos_energy(eos) result(e)
    type(eos_data_t),intent(in)::eos
    real(wp)::e(eos%n)
    e=eos%nb*(eos%q(7,:)+1._wp)*eos%mn
  end function
  pure elemental logical function ieee_finite(x)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(wp),intent(in)::x
    ieee_finite=ieee_is_finite(x)
  end function
end module compose_eos_types
