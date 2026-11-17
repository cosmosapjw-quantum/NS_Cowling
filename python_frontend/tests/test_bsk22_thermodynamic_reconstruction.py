"""Optional first-law reconstruction; no microscopic ETFSI claim."""
import math
import numpy as np
import pytest
from compose_eos import bsk22_analytic as source
from compose_eos.bsk22_thermodynamic_reconstruction import BoundaryAnchoredBSk22

@pytest.fixture
def anchored():
    return BoundaryAnchoredBSk22(outer_density=2.694320064010902e-4,
        outer_pressure=4.992059167689819e-4,outer_mu_b=939.5654205200001,
        outer_epsilon=.2526497904790187)

def test_anchor_has_correct_local_density_order_and_gibbs_continuity(anchored):
    a=anchored.anchor_diagnostics()
    assert a["inner_density"]>a["outer_density"]
    assert a["delta_dmu_dP"]<0
    assert abs(a["pressure_residual"])<1e-15
    assert abs(a["chemical_potential_residual"])<1e-12
    assert a["inner_epsilon"]>a["outer_epsilon"]
    assert a["inner_free_neutron_fraction"]>0  # not a newly solved drip state

def test_refuses_inconsistent_outer_thermodynamic_anchor():
    with pytest.raises(ValueError):
        BoundaryAnchoredBSk22(outer_density=.00027,outer_pressure=.0005,
            outer_mu_b=939.5654,outer_epsilon=1.)

def test_refuses_reversed_density_jump():
    with pytest.raises(ValueError):
        BoundaryAnchoredBSk22(outer_density=.001,outer_pressure=.0009,outer_mu_b=940.)

@pytest.mark.parametrize("bad",[float("nan"),float("inf"),-1.,.00027,.002001])
def test_refuses_outside_supported_anchored_interval(anchored,bad):
    with pytest.raises(ValueError):anchored.state(bad)

def test_firstlaw_and_reconstructed_derivative_metadata(anchored):
    for n in np.geomspace(anchored.n0*1.001,.001999,21):
        s=anchored.state(n);h=1e-4
        lo=anchored.state(n*math.exp(-h));hi=anchored.state(n*math.exp(h))
        p=n*(hi["energy_per_baryon"]-lo["energy_per_baryon"])/(2*h)
        assert p==pytest.approx(s["pressure"],rel=5e-8)
        cs=(hi["pressure"]-lo["pressure"])/(hi["epsilon"]-lo["epsilon"])
        assert cs==pytest.approx(s["sound_speed_squared_over_c_squared"],rel=5e-8)
        chi=math.log(hi["pressure"]/lo["pressure"])/math.log(hi["epsilon"]/lo["epsilon"])
        assert chi==pytest.approx(s["dlogP_dlogrho"],rel=5e-8)
        assert s["sound_speed_squared_energy_derivative"]==s["sound_speed_squared_over_c_squared"]
        assert s["energy_minus_neutron_rest_MeV"]==s["energy_per_baryon"]-source.NEUTRON_REST_MEV
        assert s["pressure"]==source.state(n)["pressure"]
        assert s["epsilon"]!=s["source_epsilon"]

def test_composition_is_preserved_and_energy_export_is_reconstructed(anchored):
    grid=np.geomspace(anchored.n0,.002,9)
    data,rows=anchored.eos_data(grid)
    old,_=source.eos_data(grid)
    for key in ("Yn","Yp","Ye","Ymu","Aav","Zav","Yh"):
        assert np.array_equal(getattr(data,key),getattr(old,key))
    assert not np.array_equal(data.Q6,old.Q6)
    assert np.allclose(data.Q1,old.Q1,rtol=0,atol=0)
    assert np.max(abs(data.Yn+data.Yp+data.Yh-1))<1e-14
