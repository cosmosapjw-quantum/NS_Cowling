"""C10: Fortran bracketed eigenfrequency polishing, preserving Cowling physics."""
from pathlib import Path
import json
import numpy as np
import pytest

from compose_eos.stellar_background import ContractEnthalpyEOS
from compose_eos.stellar_cowling_tov_native import solve_background_native
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_cowling_native import NativeCowlingDVODE

ROOT=Path(__file__).resolve().parents[2]

@pytest.fixture(scope="module")
def star():
    eos=ContractEnthalpyEOS.from_path(
        ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json")
    central=json.loads((ROOT/"data/stellar_consumers_20261117/tov/CANONICAL_1P4.json").read_text())
    return solve_background_native(eos,central["central_density_fm3"],
        surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002,native_samples=4097)


@pytest.mark.parametrize("method",("adams","bdf"))
@pytest.mark.parametrize("label,br,nodes",(
    ("f",(2250.,2500.),0),("p1",(6500.,6750.),1)))
def test_native_root_real_source_matches_brent(star,method,label,br,nodes):
    solver=NativeCowlingDVODE(star,
        CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12,
                        max_step_h=.002,root_xtol_hz=1e-6),
        method=method,samples=8193)
    expected=solver.solve_mode(br,profile_points=200,root_backend="python_brentq")
    actual=solver.solve_mode(br,profile_points=200,root_backend="fortran")
    assert abs(expected["frequency_hz"]-actual["frequency_hz"])<0.002,(label,method)
    assert actual["W_nodes"]==nodes
    assert actual["V_nodes"]==nodes
    assert abs(actual["matching_residual"])<5e-5
    assert actual["relative_join_mismatch"]<2e-5
    assert actual["root_converged"] is True
    assert actual["root_iterations"]>0
    assert actual["root_backend"]=="fortran_bracketed_brent"


def test_native_root_does_not_call_python_residual(star,monkeypatch):
    solver=NativeCowlingDVODE(star,CowlingSettings(surface_h=1e-4),
                              samples=4097)
    def forbid(*args,**kwargs):
        raise AssertionError("Fortran root polishing must own function evaluations")
    monkeypatch.setattr(solver,"matching_residual",forbid)
    actual=solver.solve_mode((2250.,2500.),profile_points=160,
                             root_backend="fortran")
    assert 2350.<actual["frequency_hz"]<2450.


@pytest.mark.parametrize("bad",((2250.,2250.),(-1.,2250.),(2200.,2250.)))
def test_native_root_invalid_or_unbracketed_fails_closed(star,bad):
    solver=NativeCowlingDVODE(star,CowlingSettings(surface_h=1e-4),
                              samples=2049)
    with pytest.raises((ValueError,RuntimeError),match="bracket"):
        solver.solve_mode(bad,profile_points=100,root_backend="fortran")


def test_native_root_rejects_unknown_backend(star):
    solver=NativeCowlingDVODE(star,CowlingSettings(surface_h=1e-4),
                              samples=2049)
    with pytest.raises(ValueError,match="root_backend"):
        solver.solve_mode((2250.,2500.),profile_points=100,
                          root_backend="some_unverified_method")


def test_native_root_rejects_active_sharp_callback_without_invoking_it(star):
    """C ABI must reject unsegmented sharp-source root searches."""
    import ctypes as ct
    from compose_eos import stellar_cowling_native as native
    solver=NativeCowlingDVODE(star,CowlingSettings(surface_h=1e-4),samples=1025)
    invoked=[]
    cbtype=ct.CFUNCTYPE(None,ct.c_double,ct.POINTER(ct.c_double),ct.POINTER(ct.c_int))

    @cbtype
    def callback(s,values,status):
        invoked.append(float(s))
        status[0]=1

    with native._LOCK:
        solver._setup()
        dll=solver._lib
        cs=ct.c_int(-1)
        dll.cowling_native_set_source_callback(
            ct.cast(callback,ct.c_void_p),ct.byref(cs))
        assert cs.value==0
        try:
            root,residual=ct.c_double(42.),ct.c_double(42.)
            niter,nfev,status=ct.c_int(-1),ct.c_int(-1),ct.c_int(-1)
            dll.cowling_native_root(
                2250.,2500.,2,1e-9,1e-12,.0001,10,
                solver.match_s,solver.surface_s,1e-5,
                ct.byref(root),ct.byref(residual),
                ct.byref(niter),ct.byref(nfev),ct.byref(status))
            assert status.value==2
            assert root.value==0. and residual.value==0.
            assert niter.value==0 and nfev.value==0
            assert invoked==[]
        finally:
            dll.cowling_native_set_source_callback(None,ct.byref(cs))
            assert cs.value==0
