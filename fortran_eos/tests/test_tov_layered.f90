program test_tov_layered
  use compose_eos_types, only: wp
  use compose_eos_tov, only: tov_love_result_t, solve_tov_love_smooth, solve_tov_love_maxwell, &
      solve_tov_love_layered, solve_tov_layered_target
  implicit none
  integer, parameter :: n=32
  integer :: i, phases(2*n)
  real(wp) :: p(2*n), e(2*n), triple_p(2*n+1), triple_e(2*n+1), scale
  integer :: triple_phase(2*n+1)
  real(wp), parameter :: pt=10._wp, pc=50._wp, step=.01_wp
  type(tov_love_result_t) :: reference, star, previous, exact, target

  do i=1,n
    p(i)=pt*exp(log(.05_wp/pt)*real(n-i,wp)/(n-1))
    p(n+i)=pt*exp(log(100._wp/pt)*real(i-1,wp)/(n-1))
  enddo
  p(n)=pt; p(n+1)=pt
  phases(:n)=1; phases(n+1:)=2
  e=100._wp+3*p
  call solve_tov_love_smooth(n,p(:n),e(:n),8._wp,reference)
  call solve_tov_love_layered(p(:n),e(:n),phases(:n),8._wp,step,star)
  call compare('smooth',star,reference,1.e-5_wp)
  if(star%transition_count/=0)error stop 'unexpected smooth interface'
  e(n+1:)=e(n+1:)+60._wp
  call solve_tov_love_maxwell(n,p(n+1:),e(n+1:),n,p(:n),e(:n),pt,60._wp,pc,reference)
  call solve_tov_love_layered(p,e,phases,pc,step,star)
  call compare('single jump',star,reference,1.e-5_wp)
  if(star%transition_count/=1.or.reference%delta_y>=0)error stop 'jump orientation'
  if(star%surface_pressure/=p(1).or.star%surface_epsilon/=e(1))error stop 'source cutoff'
  scale=star%surface_jump_scale
  if(scale<=0)error stop 'finite-density endpoint diagnostic'
  target=star
  call solve_tov_layered_target(p,e,phases,target%m_msun,5._wp,100._wp,step,star)
  if(.not.star%converged.or.abs(star%m_msun-target%m_msun)>1.e-8_wp)error stop 'target mass'
  call solve_tov_layered_target(p,e,phases,100._wp,5._wp,100._wp,step,star)
  if(star%converged.or.index(star%status,'no sampled')/=1)error stop 'missing target'

  ! A triple coexistence point splits the same total jump into two at equal P.
  triple_p=[p(:n),pt,p(n+1:)]
  triple_e=[e(:n),e(n)+15._wp,e(n+1:)]
  triple_phase=[phases(:n),3,phases(n+1:)]
  call solve_tov_love_layered(triple_p,triple_e,triple_phase,pc,step,star)
  call compare('coincident double jump',star,target,1.e-10_wp)
  if(star%transition_count/=2)error stop 'lost coincident jump'
  call solve_tov_love_layered(p(:n),e(:n),phases(:n),pt,step,exact)
  call solve_tov_love_layered(triple_p,triple_e,triple_phase,pt,step,star)
  call compare('central triple point',star,exact,1.e-10_wp)
  if(star%transition_count/=0.or.star%central_eps/=e(n))error stop 'central triple point'
  call solve_tov_love_layered(triple_p,triple_e,triple_phase,nearest(pt,1._wp),step,star)
  if(.not.star%converged.or.star%transition_count/=2)error stop 'just above crossing'
  call solve_tov_love_layered(triple_p,triple_e,triple_phase,nearest(pt,-1._wp),step,star)
  if(.not.star%converged.or.star%transition_count/=0)error stop 'just below crossing'

  previous=target
  p(n+1)=nearest(pt,-1._wp)
  call solve_tov_love_layered(p,e,phases,pc,step,star)
  call compare('roundoff on labelled plateau',star,previous,1.e-10_wp)
  p(n+1)=pt
  phases(n+1)=1
  call solve_tov_love_layered(p,e,phases,pc,step,star)
  if(star%converged)error stop 'unlabelled plateau accepted'
  ! A close increasing *same-phase* pair is a continuous cell.
  p(n+1)=pt+4*spacing(pt)
  e(n+1)=e(n)+4*spacing(e(n))
  call solve_tov_love_layered(p,e,phases,pc,step,star)
  if(.not.star%converged.or.star%transition_count/=0)error stop 'mistook continuous interval for jump'
  p(n+1)=pt
  phases(n+1)=2
  e(n+1)=e(n)-1._wp
  call solve_tov_love_layered(p,e,phases,pc,step,star)
  if(star%converged)error stop 'negative density jump accepted'
  call solve_tov_love_layered(p(:n),.5_wp*p(:n),phases(:n),8._wp,step,star)
  if(star%converged)error stop 'acausal segment accepted'
  call solve_tov_love_layered(p,e,phases,2*p(2*n),step,star)
  if(star%converged)error stop 'out of range accepted'
  write(*,'(a)') 'layered TOV controls passed'
contains
  subroutine compare(label,a,b,tol)
    character(*),intent(in) :: label
    type(tov_love_result_t),intent(in) :: a,b
    real(wp),intent(in) :: tol
    real(wp) :: errors(3)
    if(.not.a%converged.or..not.b%converged)then
      print *,label,trim(a%status),' / ',trim(b%status)
      error stop 'reference or layered solve did not converge'
    endif
    errors=abs([a%m_msun/b%m_msun-1,a%r_km/b%r_km-1,a%lambda_tidal/b%lambda_tidal-1])
    write(*,'(a,3(1x,es14.6))') label//' relative differences (M,R,Lambda)', errors
    if(any(errors>tol))error stop 'regression outside tolerance'
  end subroutine
end program
