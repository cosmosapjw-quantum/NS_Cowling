! Explicit source-prescription diagnostics and pressure-integrated BSk22 alternative.
! All energies include the SAME neutron rest-energy convention as the source fit.
! This module neither shifts a branch to match another nor claims ETFSI minimization.
module compose_eos_crust_attachment
 use compose_eos_types, only: wp
 use compose_eos_lower_crust_native, only: lower_crust_state_t,bsk_anchor_t, &
  bsk22_state,bsk22_anchored_state,BSK_DRIP,BSK_LOWER_MAX,BSK_EXT_MAX
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 type,public :: crust_attachment_t
  real(wp) :: boundary=0,left_nb=0,right_nb=0,delta_nb=0
  real(wp) :: delta_pressure=0,delta_mu=0,delta_epsilon=0
  real(wp) :: left_euler_residual=0,right_euler_residual=0
  logical :: source_prescription_valid=.false.,density_overlap=.false.,density_gap=.false.
  logical :: cross_interval_supported=.false.,equilibrium_established=.false.
 end type
 type,public :: crust_integral_budget_t
  real(wp) :: anchor_nb=0,anchor_e=0,anchor_error=0,relative_pressure_envelope=0
  real(wp) :: pressure_integral=0,quadrature_estimate=0
  real(wp) :: energy_per_baryon_model_bound=0,mu_model_bound=0,epsilon_model_bound=0
  real(wp) :: epsilon_quadrature_estimate=0,mu_quadrature_estimate=0
  real(wp) :: source_delta_energy_per_baryon=0,composition_mu_residual=0
  logical :: bound_is_conditional=.true.
 end type
 public :: source_attachment_diagnostics,source_attachment_supports
 public :: bsk22_reference_integral
contains
 ! Each endpoint is [n (fm^-3), epsilon (MeV fm^-3), P (MeV fm^-3), mu (MeV)].
 ! Equality of P and mu is reported, never required to accept a source prescription.
 ! The unsampled open interval remains unsupported, even if endpoints agree.
 subroutine source_attachment_diagnostics(left,right,boundary,diagnostic,status)
  real(wp),intent(in)::left(4),right(4),boundary
  type(crust_attachment_t),intent(out)::diagnostic
  character(len=*),intent(out)::status
  diagnostic=crust_attachment_t();status='invalid source attachment endpoints'
  if(.not.all(ieee_is_finite(left)).or..not.all(ieee_is_finite(right)))return
  if(.not.ieee_is_finite(boundary))return
  if(any(left<=0).or.any(right<=0).or.boundary<=0)return
  if(boundary<min(left(1),right(1)).or.boundary>max(left(1),right(1)))then
   status='source endpoints do not bracket nominal density boundary';return
  endif
  diagnostic%boundary=boundary;diagnostic%left_nb=left(1);diagnostic%right_nb=right(1)
  diagnostic%delta_nb=right(1)-left(1)
  diagnostic%density_overlap=right(1)<left(1);diagnostic%density_gap=right(1)>left(1)
  diagnostic%delta_epsilon=right(2)-left(2)
  diagnostic%delta_pressure=right(3)-left(3);diagnostic%delta_mu=right(4)-left(4)
  diagnostic%left_euler_residual=left(2)+left(3)-left(1)*left(4)
  diagnostic%right_euler_residual=right(2)+right(3)-right(1)*right(4)
  diagnostic%source_prescription_valid=.true.;status='ok'
 end subroutine

 ! This query concerns the cross-seam interval only; it does not certify remote support.
 logical function source_attachment_supports(diagnostic,n) result(supported)
  type(crust_attachment_t),intent(in)::diagnostic
  real(wp),intent(in)::n
  supported=.false.
  if(.not.diagnostic%source_prescription_valid.or..not.ieee_is_finite(n))return
  supported=n==diagnostic%left_nb.or.n==diagnostic%right_nb
 end function

 ! Fixed anchor n0: e_I(n)=e_source(n0)+integral[n0,n] P_source(s)/s^2 ds.
 ! epsilon_I=n e_I, mu_I=e_I+P_source/n. No fitted matching constant is introduced.
 ! The supplied envelope r means |P_source-P_true| <= r P_true, 0<=r<1,
 ! at every SAME baryon density. This must be justified by the caller, not inferred
 ! from the paper's empirical errors at equal mass-energy density.
 ! Conditional model bounds exclude quadrature error, tabulated separately.
 subroutine bsk22_reference_integral(n0,n,anchor_error,pressure_envelope,state,budget,status,epsabs,extended)
  real(wp),intent(in)::n0,n,anchor_error,pressure_envelope
  type(lower_crust_state_t),intent(out)::state
  type(crust_integral_budget_t),intent(out)::budget
  character(len=*),intent(out)::status
  real(wp),intent(in),optional::epsabs
  logical,intent(in),optional::extended
  type(bsk_anchor_t)::anchor
  type(lower_crust_state_t)::source_anchor
  real(wp)::upper,tolerance,factor,increment
  logical::ext
  state=lower_crust_state_t();budget=crust_integral_budget_t()
  status='invalid pressure-integral reference or uncertainty'
  if(.not.all(ieee_is_finite([n0,n,anchor_error,pressure_envelope])))return
  if(anchor_error<0.or.pressure_envelope<0.or.pressure_envelope>=1)return
  tolerance=1e-10_wp
  if(present(epsabs))tolerance=epsabs
  if(.not.ieee_is_finite(tolerance).or.tolerance<=0)return
  ext=.false.;if(present(extended))ext=extended
  upper=BSK_LOWER_MAX;if(ext)upper=BSK_EXT_MAX
  if(n0<BSK_DRIP.or.n0>upper.or.n<n0.or.n>upper)return
  call bsk22_state(n0,'source_fit',source_anchor,status,ext)
  if(status/='ok')return
  anchor%n0=n0;anchor%e0=source_anchor%epsilon/n0
  anchor%epsabs=tolerance;anchor%valid=.true.
  call bsk22_anchored_state(anchor,n,state,status,ext)
  if(status/='ok')return
  state%representation='source-energy-anchored pressure integral; analytic composition'
  increment=state%epsilon/n-anchor%e0
  ! Positive pressure and n>=n0 make the exact integral nonnegative.
  factor=pressure_envelope/(1-pressure_envelope)
  budget%anchor_nb=n0;budget%anchor_e=anchor%e0;budget%anchor_error=anchor_error
  budget%relative_pressure_envelope=pressure_envelope
  budget%pressure_integral=increment;budget%quadrature_estimate=state%quadrature_error
  budget%energy_per_baryon_model_bound=anchor_error+factor*max(0._wp,increment)
  budget%mu_model_bound=budget%energy_per_baryon_model_bound+factor*state%pressure/n
  budget%epsilon_model_bound=n*budget%energy_per_baryon_model_bound
  ! Include integral uncertainty in both the reconstructed state and model-bound estimate.
  ! Adaptive Simpson estimates are numerical estimates, not certified enclosures.
  budget%epsilon_quadrature_estimate=n*(1+factor)*state%quadrature_error
  budget%mu_quadrature_estimate=(1+factor)*state%quadrature_error
  budget%source_delta_energy_per_baryon=(state%epsilon-state%source_epsilon)/n
  budget%composition_mu_residual=state%composition_mu_residual
 end subroutine
end module compose_eos_crust_attachment
