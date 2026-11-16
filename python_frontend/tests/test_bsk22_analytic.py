"""Published-fit implementation checks, independent of CompOSE target matching."""
import math
import numpy as np
import pytest
from compose_eos import bsk22_analytic as bsk

@pytest.mark.parametrize("n",[0.,-1.,float("nan"),float("inf"),bsk.NEUTRON_DRIP_DENSITY*(1-1e-10),.002*(1+1e-10)])
def test_refuses_outside_declared_lower_inner_support(n):
    with pytest.raises(ValueError):bsk.state(n)

def test_unknown_thermodynamic_representation_is_not_silently_selected():
    with pytest.raises(ValueError):bsk.state(.001,thermodynamic_mode="automatic")

def test_composition_preserves_cell_cluster_and_free_neutron_distinction():
    for n in np.geomspace(bsk.NEUTRON_DRIP_DENSITY,.002,97):
        s=bsk.state(n)
        assert s["Y_p_free"]==0 and s["Y_mu"]==0 and s["Z_cluster"]==40
        assert 0<=s["Y_n_free"]<1
        assert s["A_cell"]>=s["A_cluster"]>s["Z_cluster"]
        assert abs(s["Y_n_free"]+s["Y_p_free"]+s["Y_h"]-1)<1e-14
        assert abs(s["Y_e"]-s["Y_p_free"]-s["Z_cluster"]/s["A_cluster"]*s["Y_h"])<1e-14
    drip=bsk.state(bsk.NEUTRON_DRIP_DENSITY)
    assert drip["Y_n_free"]==0
    assert drip["A_cell"]==drip["A_cluster"]

@pytest.mark.parametrize("mode",["source_fit","energy_derivative"])
def test_analytic_derivatives_against_finite_differences(mode):
    for n in np.geomspace(bsk.NEUTRON_DRIP_DENSITY*1.001,.002*.999,31):
        h=1e-5;s=bsk.state(n,thermodynamic_mode=mode)
        lo=bsk.state(n*math.exp(-h),thermodynamic_mode=mode)
        hi=bsk.state(n*math.exp(h),thermodynamic_mode=mode)
        p=n*(hi["energy_minus_neutron_rest_MeV"]-lo["energy_minus_neutron_rest_MeV"])/(2*h)
        assert p==pytest.approx(s["pressure_from_energy_derivative"],rel=3e-8)
        cs=(hi["pressure"]-lo["pressure"])/(hi["epsilon"]-lo["epsilon"])
        assert cs==pytest.approx(s["sound_speed_squared_over_c_squared"],rel=3e-8)
        assert 0<s["sound_speed_squared_over_c_squared"]<1

def test_two_representations_report_the_source_fit_thermodynamic_defect():
    s=bsk.state(.0003)
    t=bsk.state(.0003,thermodynamic_mode="energy_derivative")
    assert s["epsilon"]==t["epsilon"] and s["Y_e"]==t["Y_e"]
    assert s["pressure"]!=t["pressure"]
    assert t["pressure"]==t["pressure_from_energy_derivative"]
    assert t["mu_b_euler"]==pytest.approx(t["mu_b_energy_derivative"],abs=2e-12)
    assert abs(s["pressure_fit_thermodynamic_relative_residual"])>.01

def test_independent_author_fortran_anchor():
    # bskfit18.f(2023-02-13), SHA256 recorded in module; compiled without
    # changing its coefficients. Paper7/6 vs source1.16667 and default-real
    # literals explain the deliberately bounded non-bitwise difference.
    s=bsk.state(.002)
    assert s["pressure"]==pytest.approx(.0021738420487500206,rel=1e-6)
    assert s["energy_minus_neutron_rest_MeV"]==pytest.approx(.7579807622156363,abs=1.1e-4)
    assert s["A_cell"]==pytest.approx(518.5899964582201,rel=1e-7)
    assert s["A_cluster"]==pytest.approx(150.29675385687733,rel=2e-7)
    assert s["Y_n_free"]==pytest.approx(.7101819262165696,abs=3e-8)

def test_compose_constructor_uses_free_fractions_and_bound_cluster_mass():
    data,rows=bsk.eos_data(np.geomspace(bsk.NEUTRON_DRIP_DENSITY,.002,9))
    assert np.all(data.Yp==0)
    assert data.Yn[-1]>.7
    assert np.allclose(data.Yn+data.Yp+data.Yh,1,rtol=0,atol=2e-15)
    assert np.allclose(data.Ye,data.Yp+data.Zav/data.Aav*data.Yh,rtol=0,atol=2e-15)
    assert np.all(data.Yn_present) and np.all(data.Yp_present)
    assert np.all(data.row_status==2)
