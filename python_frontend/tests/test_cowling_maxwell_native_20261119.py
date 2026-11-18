"""Native Fortran sharp Maxwell TOV + Cowling, slow and rapid conversion."""
import json
from pathlib import Path

import numpy as np
import pytest

from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import (
    SharpContractEnthalpyEOS, solve_sharp_background, transfer_interface,
)
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp_native import (
    solve_sharp_background_native, NativeSharpCowlingDVODE,
    native_interface_transfer,
)

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/cubic_quintic_matching_20261118/native_stellar/SPG_M1__PK1_CANONICAL/native_maxwell/REPRESENTATION.json"


@pytest.fixture(scope="module")
def eos():
    contract=json.loads(SOURCE.read_text())
    return SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(contract))


@pytest.fixture(scope="module")
def stars(eos):
    kw=dict(surface_h=3e-8,rtol=2e-9,atol=2e-11,max_step_h=.005)
    py=solve_sharp_background(eos,.3,**kw)
    native=solve_sharp_background_native(eos,.3,native_samples=4097,**kw)
    return py,native


def test_native_two_segment_tov_parity(stars):
    py,native=stars
    assert len(native.interfaces)==len(py.interfaces)==1
    assert abs(py.mass_msun-native.mass_msun)<1e-4
    assert abs(py.radius_km-native.radius_km)<1e-3
    assert abs(py.interfaces[0]["r_km"]-native.interfaces[0]["r_km"])<1e-3
    t=native.interfaces[0]
    assert t["n_high_fm3"]>t["n_low_fm3"]


def test_native_transfer_matches_independent_python_both_conversions(stars):
    _,star=stars
    h=star.interfaces[0]["h"]
    inner=star.state_at_enthalpy(h,"high")
    outer=star.state_at_enthalpy(h,"low")
    omega=2*np.pi*2300/299792.458
    y=np.array([0.8,-0.43],dtype=float)
    for conv in ("slow","rapid"):
        actual,diag=native_interface_transfer(y,inner,outer,omega,conv)
        expected=transfer_interface(y,inner,outer,omega,conv)
        np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-12)
        assert diag["traction_relative_residual"]<2e-12
        assert diag["kinematic_relative_residual"]<2e-12
    with pytest.raises(ValueError,match="conversion"):
        native_interface_transfer(y,inner,outer,omega,"invalid")


def test_native_sharp_cowling_determinant_sign_matches_python(stars):
    # Determinant *normalization* is match-point dependent; zeros and
    # local signs are physically comparable, not their magnitudes.
    from compose_eos.stellar_sharp import matching_residual
    py,native=stars
    settings=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)
    solver=NativeSharpCowlingDVODE(native,settings,samples=4097)
    for conv in ("slow","rapid"):
        for freq in (1800.,2400.,5600.):
            z=solver.matching_residual(freq,conversion=conv)
            expected=matching_residual(py,freq,settings,conversion=conv)
            assert np.isfinite(z) and abs(z)<=1
            assert np.sign(z)==np.sign(expected),(conv,freq,z,expected)
    assert solver.max_interface_residual<1e-8


def test_native_sharp_crossing_mode_nodes_and_surface(stars):
    _,native=stars
    settings=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)
    solver=NativeSharpCowlingDVODE(native,settings,samples=4097)
    for conv in ("slow","rapid"):
        brackets,_=solver.scan_brackets(np.linspace(1000,9000,50),conversion=conv)
        assert brackets,conv
        modes=[solver.solve_mode(br,conversion=conv,profile_points=160)
               for br in brackets[:2]]
        assert any(abs(x["matching_residual"])<1e-5 for x in modes)
        for x in modes:
            assert x["surface_BC_relative"]<1e-12
            assert x["interface_junctions"]
            assert all(z["traction_relative_residual"]<1e-8 for z in x["interface_junctions"])
            assert all(z["kinematic_relative_residual"]<1e-8 for z in x["interface_junctions"])
