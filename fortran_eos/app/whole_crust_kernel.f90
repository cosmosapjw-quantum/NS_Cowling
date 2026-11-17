! Explicit thermodynamic pieces for the reference-assisted whole-crust baseline.
! Units: n [fm^-3], P and epsilon [MeV fm^-3], mu [MeV], cs2=dP/depsilon.
! The Python frontend records origins and independently checks these formulae.
program whole_crust_kernel
  use, intrinsic :: iso_fortran_env, only: real64, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: wp=real64
  character(len=64) :: mode
  character(len=2048) :: input_path, output_path
  integer :: iu, ou, ios, count, j
  real(wp), allocatable :: nb(:), pr(:), en(:), ff(:), gg(:)
  real(wp) :: n, p, e, mu, cs, gamma, a, cutoff, na, pa, ea
  real(wp) :: nl, pl, el, nr, pp, er, ml, mr, h, c2, c3, t, dndmu
  real(wp) :: lr, z, deltaf, mnmin, dnmin, cqmin, n1, n2, q0, q1, q2, cqbound
  if(command_argument_count()/=3) call fail('usage: whole_crust_kernel MODE INPUT OUTPUT')
  call get_command_argument(1,mode)
  call get_command_argument(2,input_path)
  call get_command_argument(3,output_path)
  open(newunit=iu,file=trim(input_path),status='old',action='read',iostat=ios)
  if(ios/=0) call fail('cannot read input')
  ! Never overwrite an existing output: a failed attempt remains distinguishable.
  open(newunit=ou,file=trim(output_path),status='new',action='write',iostat=ios)
  if(ios/=0) call fail('output exists or cannot be created')
  write(ou,'(a)') 'n_b_fm3,pressure_mev_fm3,epsilon_mev_fm3,mu_mev,cs2_over_c2'
  select case(trim(mode))
  case('pressure-primary')
    read(iu,*,iostat=ios) count
    if(ios/=0.or.count<2) call fail('at least two source knots required')
    allocate(nb(count),pr(count),en(count),ff(count),gg(count-1))
    do j=1,count
      read(iu,*,iostat=ios) nb(j),pr(j),en(j)
      if(ios/=0) call fail('invalid source row')
      if(.not.all(ieee_is_finite([nb(j),pr(j),en(j)]))) call fail('nonfinite source')
      if(min(nb(j),pr(j),en(j))<=0._wp) call fail('positive source quantities required')
      if(j>1) then
        if(nb(j)<=nb(j-1).or.pr(j)<=pr(j-1)) call fail('strictly increasing n and P required')
      endif
    enddo
    ff(1)=en(1)/nb(1)
    do j=1,count-1
      lr=log(nb(j+1)/nb(j))
      gg(j)=log(pr(j+1)/pr(j))/lr
      z=(gg(j)-1._wp)*lr
      deltaf=(pr(j)/nb(j))*lr*exprel(z)
      ff(j+1)=ff(j)+deltaf
      ! cs2 is monotone within a power-law pressure piece; both endpoints suffice.
      cs=gg(j)*pr(j)/(nb(j)*ff(j)+pr(j))
      if(cs<=0._wp.or.cs>1._wp) call fail('noncausal pressure-law left endpoint')
      cs=gg(j)*pr(j+1)/(nb(j+1)*ff(j+1)+pr(j+1))
      if(cs<=0._wp.or.cs>1._wp) call fail('noncausal pressure-law right endpoint')
    enddo
    do j=1,count
      e=nb(j)*ff(j); mu=ff(j)+pr(j)/nb(j)
      cs=gg(min(j,count-1))*pr(j)/(e+pr(j))
      call row(nb(j),pr(j),e,mu,cs)
    enddo
  case('surface-polytrope')
    read(iu,*,iostat=ios) na,pa,ea,gamma,cutoff,count
    if(ios/=0) call fail('invalid surface input')
    if(.not.all(ieee_is_finite([na,pa,ea,gamma,cutoff]))) call fail('nonfinite surface input')
    if(min(na,pa,ea,cutoff)<=0._wp.or.cutoff>=na.or.count<2) call fail('invalid surface range')
    if(gamma<=1._wp.or.gamma>2._wp) call fail('surface gamma must be in (1,2]')
    a=(ea-pa/(gamma-1._wp))/na
    if(a<=0._wp) call fail('nonpositive surface chemical-potential limit')
    do j=1,count
      n=exp(log(cutoff)+(log(na)-log(cutoff))*real(j-1,wp)/real(count-1,wp))
      if(j==1) n=cutoff
      if(j==count) n=na
      p=pa*(n/na)**gamma
      e=a*n+p/(gamma-1._wp)
      mu=a+gamma*p/((gamma-1._wp)*n)
      cs=gamma*p/(e+p)
      call row(n,p,e,mu,cs)
    enddo
  case('hermite-bridge')
    read(iu,*,iostat=ios) nl,pl,el,nr,pp,er,count
    if(ios/=0) call fail('invalid bridge input')
    if(.not.all(ieee_is_finite([nl,pl,el,nr,pp,er]))) call fail('nonfinite bridge endpoints')
    if(min(nl,pl,el,nr,pp,er)<=0._wp.or.count<2) call fail('positive bridge endpoints required')
    if(nr<=nl.or.pp<=pl.or.er<=el) call fail('bridge endpoints must increase')
    ml=(el+pl)/nl; mr=(er+pp)/nr; h=mr-ml
    if(h<=0._wp) call fail('chemical potential must increase')
    c2=3._wp*(pp-pl)-h*(2._wp*nl+nr)
    c3=2._wp*(pl-pp)+h*(nl+nr)
    n1=2._wp*c2/h; n2=3._wp*c3/h
    mnmin=quadratic_min(nl,n1,n2)
    dnmin=min(2._wp*c2,2._wp*c2+6._wp*c3)/h**2
    q0=ml*2._wp*c2/h**2-nl
    q1=ml*6._wp*c3/h**2
    q2=3._wp*c3/h
    cqmin=quadratic_min(q0,q1,q2)
    ! At a luminal analytic boundary mu*n_mu and n cancel. Scale the
    ! representation error BEFORE subtraction, never by the cancelled residual.
    cqbound=64._wp*epsilon(1._wp)*(abs(ml*2._wp*c2/h**2)+abs(nl)+ &
        abs(ml*6._wp*c3/h**2)+abs(3._wp*c3/h))
    if(mnmin<=0._wp.or.dnmin<=0._wp.or.cqmin < -cqbound) call fail('bridge exact interval bounds fail')
    do j=1,count
      t=real(j-1,wp)/real(count-1,wp); mu=ml+h*t
      p=pl+h*nl*t+c2*t*t+c3*t*t*t
      n=nl+n1*t+n2*t*t
      dndmu=(2._wp*c2+6._wp*c3*t)/h**2
      e=mu*n-p; cs=n/(mu*dndmu)
      ! Preserve the supplied physical endpoints, not roundoff in polynomial evaluation.
      if(j==1) then
        n=nl; p=pl; e=el; mu=ml
      elseif(j==count) then
        n=nr; p=pp; e=er; mu=mr
      endif
      call row(n,p,e,mu,cs)
    enddo
  case default
    call fail('unknown construction mode')
  end select
  close(iu); close(ou)
  if(trim(mode)=='hermite-bridge') then
    write(*,'(a,es25.16e3,a,es25.16e3,a)') &
      '{"status":"ok","backend":"fortran","operation":"whole-crust-piece","minimum_causal_margin":', &
      cqmin,',"causal_roundoff_bound":',cqbound,'}'
  else
    write(*,'(a)') '{"status":"ok","backend":"fortran","operation":"whole-crust-piece"}'
  endif
contains
  subroutine fail(message)
    character(len=*),intent(in) :: message
    write(error_unit,'(a)') trim(message)
    stop 1
  end subroutine fail
  subroutine row(n,p,e,mu,cs)
    real(wp),intent(in) :: n,p,e,mu,cs
    if(.not.all(ieee_is_finite([n,p,e,mu,cs]))) call fail('nonfinite constructed row')
    if(min(n,p,e,mu,cs)<=0._wp.or.cs>1._wp+64._wp*epsilon(1._wp)) call fail('nonphysical constructed row')
    write(ou,'(es26.17e3,4(",",es26.17e3))') n,p,e,mu,cs
  end subroutine row
  real(wp) function exprel(z) result(v)
    real(wp),intent(in) :: z
    if(abs(z)<1.e-5_wp) then
      v=1._wp+z*(.5_wp+z*(1._wp/6._wp+z*(1._wp/24._wp+z*(1._wp/120._wp+z/720._wp))))
    else
      v=(exp(z)-1._wp)/z
    endif
  end function exprel
  real(wp) function quadratic_min(a,b,c) result(v)
    real(wp),intent(in) :: a,b,c
    real(wp) :: x
    v=min(a,a+b+c)
    if(c>0._wp) then
      x=-b/(2._wp*c)
      if(x>0._wp.and.x<1._wp) v=min(v,a+b*x+c*x*x)
    endif
  end function quadratic_min
end program whole_crust_kernel
