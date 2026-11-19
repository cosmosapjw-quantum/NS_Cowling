"""Stage-05 thermodynamic reconciliation controls; no source mutation."""
import numpy as np
import pytest

from compose_eos.rmf_source_reconciliation import (
    cold_pressure_from_energy,
    compare_phase_homogeneous_states,
    rest_mass_shift_pressure_invariant,
)


def test_exact_cold_polytrope_derivative_pressure():
    n=np.geomspace(.025,.45,44)
    mass=939.
    K=400.
    eps=mass*n+K*n*n
    p=K*n*n
    got=cold_pressure_from_energy(n,eps,p)
    np.testing.assert_allclose(got["P_from_energy"],p,rtol=5e-9,atol=2e-10)
    assert got["max_inner_abs_P_error"]<2e-8
    assert got["max_provider_Q3_Euler_discrepancy_MeV"] is None


def test_rest_mass_zero_point_does_not_change_pressure():
    n=np.geomspace(.03,.36,51)
    eps=939*n+100*n*n
    p=100*n*n
    res=rest_mass_shift_pressure_invariant(n,eps,p,delta_mass_MeV=1.)
    np.testing.assert_allclose(res["P_after_shift"],p,atol=0,rtol=0)
    np.testing.assert_allclose(res["epsilon_after_shift"],eps+n,atol=0,rtol=0)
    assert res["max_abs_pressure_shift"]==0


def test_thermodynamic_reconstruction_catches_nonmonotone_phase():
    n=[.02,.04,.05,.07]
    e=np.array(n)*939
    with pytest.raises(ValueError,match="strictly increasing"):
        cold_pressure_from_energy(n[::-1],e[::-1],np.array(n[::-1])**2)


def test_composition_compatible_is_distinct_from_pressure_identity():
    n=np.array([.11,.18,.3])
    src=dict(n=n,P=np.array([1.,5.,40.]),
             epsilon=939*n+2*n*n,mu_Q3=np.array([940.,945.,1000.]),
             Yp=np.array([.04,.09,.14]),Ye=np.array([.04,.07,.09]),
             Ymu=np.array([0,.02,.05]))
    model=[]
    for i,ni in enumerate(n):
        model.append(dict(n=float(ni),P=float(src["P"][i]*1.04),
            epsilon=float(src["epsilon"][i]),mu_B=float(src["mu_Q3"][i]+1.),
            Yp=float(src["Yp"][i]+2e-6),Ye=float(src["Ye"][i]+3e-6),
            Ymu=float(src["Ymu"][i]-1e-6)))
    out=compare_phase_homogeneous_states(src,model,rounding_half_width=5.01e-6)
    assert out["composition_compatible_with_printed_rounding"] is True
    assert abs(out["pressure_max_abs_rel"]-.04)<2e-14
    assert out["bulk_pressure_source_identical"] is False
    assert out["pressure_max_abs_rel"]>.01


def test_reject_physical_output_at_other_density():
    n=np.array([.1,.2,.3])
    src=dict(n=n,P=n*n,epsilon=939*n,mu_Q3=np.full(3,939.),
             Yp=np.full(3,.1),Ye=np.full(3,.1),Ymu=np.zeros(3))
    model=[dict(n=float(ni+.01),P=.1,epsilon=100.,mu_B=960.,
                Yp=.1,Ye=.1,Ymu=0.) for ni in n]
    with pytest.raises(ValueError,match="density"):
        compare_phase_homogeneous_states(src,model)


def test_compose_q4_free_electron_consistency_and_rounding_gate():
    from compose_eos.cold_leptons import ELECTRON_MASS, free_fermion
    from compose_eos.rmf_source_reconciliation import charge_mu_consistency
    n=np.array([.09,.11,.16,.20])
    ye=np.array([.035,.045,.055,.065])
    mn=938.
    table_mu_e=np.array([free_fermion(float(ni*yi),ELECTRON_MASS)[2]
                          for ni,yi in zip(n,ye)])
    q4=-table_mu_e/mn
    q5=np.zeros_like(q4)
    exact=charge_mu_consistency(n,ye,q4,q5,mn=mn,
                                source_fraction_rounding_half_width=5e-6)
    assert exact["all_rows_rounding_consistent"] is True
    assert exact["max_abs_mu_e_discrepancy_MeV"]<2e-12
    altered=charge_mu_consistency(n,ye,q4+2./mn,q5,mn=mn,
                                 source_fraction_rounding_half_width=5e-6)
    assert altered["all_rows_rounding_consistent"] is False
    assert altered["max_abs_mu_e_discrepancy_MeV"]>1.99


def test_compose_q5_conservation_term_in_mu_e_accepted():
    from compose_eos.cold_leptons import ELECTRON_MASS, free_fermion
    from compose_eos.rmf_source_reconciliation import charge_mu_consistency
    n=np.array([.09,.11,.16,.20]);ye=np.array([.035,.045,.055,.065])
    mn=939.
    exact=np.array([free_fermion(float(ni*yi),ELECTRON_MASS)[2]
                    for ni,yi in zip(n,ye)])
    q5=np.full(4,10./mn)
    q4=q5-exact/mn
    found=charge_mu_consistency(n,ye,q4,q5,mn=mn)
    assert found["all_rows_rounding_consistent"] is True


def test_pressure_mismatch_decomposition_separates_energy_derivative_from_source_failure():
    from compose_eos.rmf_source_reconciliation import decompose_pressure_difference
    n=np.linspace(.025,.45,64)
    m=939.;ks=200.;km=240.;unphysical_extra=.03
    energy_provider=m*n+ks*n**2
    pressure_provider=ks*n**2+unphysical_extra*n
    energy_model=m*n+km*n**2
    pressure_model=km*n**2
    actual=decompose_pressure_difference(
        n,pressure_provider,energy_provider,pressure_model,energy_model)
    physical=actual["pressure_from_energy_difference"]
    np.testing.assert_allclose(physical,(km-ks)*n**2,rtol=8e-12,atol=1e-10)
    np.testing.assert_allclose(actual["provider_identity_residual"],
                               unphysical_extra*n,rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(actual["model_identity_residual"],
                               np.zeros_like(n),rtol=0,atol=1e-10)
    np.testing.assert_allclose(actual["pressure_observed_difference"],
                               physical-actual["provider_identity_residual"]
                               +actual["model_identity_residual"],rtol=2e-12,atol=1e-10)
    assert actual["max_abs_balance_error_MeV_fm3"]<2e-10
    assert actual["rms_model_source_pressure_difference_MeV_fm3"]>0


def test_pressure_mismatch_decomposition_rejects_a_phase_junction():
    from compose_eos.rmf_source_reconciliation import decompose_pressure_difference
    n=np.array([.02,.03,.04,.04,.05,.06])
    P=np.array([.01,.02,.03,.03,.04,.06])
    eps=939*n+P
    with pytest.raises(ValueError,match="strictly increasing"):
        decompose_pressure_difference(n,P,eps,P,eps)
