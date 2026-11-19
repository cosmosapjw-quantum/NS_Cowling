"""Source-calibration contracts: original bytes are immutable, no phase smearing."""
import numpy as np
import pytest

from compose_eos.source_calibration import (
    first_law_segment_audit, source_model_row_audit, phase_contiguous_slices,
)


def test_first_law_exact_gamma2_polytrope():
    n=np.linspace(.02,.4,52)
    k,m=200.,939.
    P=k*n*n
    eps=m*n+P
    mu=m+2*k*n
    report=first_law_segment_audit(n,P,eps,mu)
    assert report["valid_pressure_primary"]
    assert report["max_abs_delta_epsilon_per_baryon_MeV"]<2.e-10
    assert report["max_abs_mu_Q3_minus_Euler_MeV"]<1e-10
    np.testing.assert_allclose(report["corrected_pressure"],P,rtol=0,atol=0)
    np.testing.assert_allclose(report["corrected_epsilon"],eps,rtol=2e-14)
    assert report["max_cs2_over_c2"]<1


def test_first_law_requires_explicit_segment_split_not_across_jump():
    n=np.array([.01,.02,.03,.05,.06,.08])
    P=np.array([.001,.002,.003,.003,.004,.005])
    e=939*n+P
    with pytest.raises(ValueError,match="strictly increasing pressure"):
        first_law_segment_audit(n,P,e,(e+P)/n)
    phase=np.array([7,7,7,6,6,6])
    assert phase_contiguous_slices(phase)==[(7,0,3),(6,3,6)]


def test_phase_contiguous_slices_preserves_disconnected_same_id():
    assert phase_contiguous_slices([0,0,7,7,0,0])==[
        (0,0,2),(7,2,4),(0,4,6)]


def test_same_density_comparison_separates_pressure_from_energy():
    n=np.array([.05,.1,.2])
    P=np.array([.01,.1,.5])
    eps=np.array([47.,95.,195.])
    mu=(eps+P)/n
    provider=dict(n=n,P=P,epsilon=eps,mu_Q3=mu)
    # Measured differences at exact source density, no grid interpolation.
    model=[dict(P=float(p*1.01),epsilon=float(e*1.001),
                mu_B=float(v+2.)) for p,e,v in zip(P,eps,mu)]
    x=source_model_row_audit(provider,model)
    assert abs(x["pressure_abs_rel_max"]-.01)<1e-14
    assert abs(x["epsilon_abs_rel_max"]-.001)<1e-14
    assert x["mu_model_minus_provider_q3_abs_max_MeV"]==2.
    assert x["n_evaluated"]==3
    assert all(y["density_fm3"]==z for y,z in zip(x["samples"],n))


def test_rejects_nonfinite_and_shape_inconsistent_source_inputs():
    n=np.array([.01,.02,.03])
    P=np.array([.01,.03,.04])
    eps=np.array([9.4,19.,28.5])
    with pytest.raises(ValueError):
        first_law_segment_audit(n,np.array([.01,float("nan"),.04]),eps,(P+eps)/n)
    with pytest.raises(ValueError):
        source_model_row_audit(dict(n=n,P=P,epsilon=eps,mu_Q3=(P+eps)/n),[dict(P=.1,epsilon=1.,mu_B=939.)])
