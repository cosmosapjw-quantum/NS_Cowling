"""C13: native geometry PCHIP coefficients and profile radius from Fortran.

The SciPy PchipInterpolator is an independent reference only. Production
native_source_cs2 must not invoke Python PCHIP for its 5 background columns.
Original EOS and Maxwell one-sided identities remain immutable.
"""
from __future__ import annotations
import ctypes as ct
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from compose_eos import stellar_cowling_native as native
from compose_eos import stellar_sharp_native as sharp
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
CASES=(
 ("SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.),5143.719608336112),
 ("SPG_M3","css",.192510767556007,"slow",(3710.,3760.),3736.210698739710)
)

@pytest.fixture(scope="module",params=CASES)
def physical(request):
    model,join,nc,conversion,br,ref=request.param
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
    with zipfile.ZipFile(SOURCE) as z:
        with zipfile.ZipFile(io.BytesIO(z.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            raw=inner.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    bg=solve_sharp_background_native(eos,nc,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    opt=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,
                        max_step_h=.002,root_xtol_hz=1e-5)
    return model,bg,opt,conversion,br,ref


@pytest.mark.parametrize("case",["monotone","flat","nonmonotone"])
def test_native_pchip_matches_scipy_nonuniform_endpoint_and_flat(case):
    s=np.ascontiguousarray([0.,1e-5,2e-4,2e-3,.125,.3,.9,1.1])
    y=np.zeros((5,len(s)),dtype=np.float64,order="F")
    if case=="monotone":
        for k in range(5):y[k]=np.exp((k+1)*s)
    if case=="flat":
        for k in range(5):y[k]=[1.,1.,1.5,1.5,1.5,2.,2.,2.]
    if case=="nonmonotone":
        for k in range(5):y[k]=np.sin((k+1)*s)+.03*s
    dll=native._native_library()
    slopes=np.empty_like(y,order="F")
    st=ct.c_int(-1)
    dll.cowling_native_pchip_slopes(len(s),s,y,slopes,ct.byref(st))
    assert st.value==0
    oracle=PchipInterpolator(s,y.T,axis=0,extrapolate=False)
    np.testing.assert_allclose(slopes.T,oracle(s,nu=1),
       atol=1e-10,rtol=2e-12)
    xs=np.ascontiguousarray(np.linspace(s[0],s[-1],92))
    values=np.empty((5,len(xs)),order="F")
    dll.cowling_native_pchip_values(len(s),s,y,slopes,len(xs),xs,values,ct.byref(st))
    assert st.value==0
    np.testing.assert_allclose(values.T,oracle(xs),atol=3e-12,rtol=1e-11)


@pytest.mark.parametrize("mutate",["duplicate","nonfinite","outside_query"])
def test_native_pchip_fails_closed_on_invalid_knots(mutate):
    s=np.ascontiguousarray([0.,.1,.4,.7])
    y=np.asfortranarray(np.vstack([s+1 for _ in range(5)]))
    if mutate=="duplicate":s[2]=s[1]
    if mutate=="nonfinite":y[0,1]=np.nan
    slopes=np.full_like(y,19.,order="F")
    dll=native._native_library();st=ct.c_int(-1)
    dll.cowling_native_pchip_slopes(len(s),s,y,slopes,ct.byref(st))
    if mutate!="outside_query":
        assert st.value!=0
        assert np.all(slopes==0.)
    else:
        assert st.value==0
        xs=np.ascontiguousarray([-.001,.2])
        v=np.full((5,len(xs)),88.,order="F")
        dll.cowling_native_pchip_values(len(s),s,y,slopes,len(xs),xs,v,ct.byref(st))
        assert st.value!=0
        assert np.all(v==0.)


def test_native_geometry_slopes_match_scipy_on_source_locked_star(physical):
    model,bg,opt,conversion,br,ref=physical
    solver=NativeSharpCowlingDVODE(bg,opt,samples=2049,
           coefficient_backend="native_source_cs2",
           source_restart_relative_jump=.01)
    assert all(p["pchip"] is None for p in solver._phases)
    assert solver.geometry_interpolation_backend=="fortran_pchip"
    assert solver.source_restart_count in (39,40)
    for phase in solver._phases:
        s=phase["s"];v=phase["vals"]
        baseline=PchipInterpolator(s,v.T,axis=0,extrapolate=False)
        expected=np.asfortranarray(baseline(s,nu=1).T)
        np.testing.assert_allclose(phase["slopes"],expected,
             rtol=2e-10,atol=5e-9)
        coords=np.ascontiguousarray(np.linspace(s[0],s[-1],83))
        sample=solver._native_phase_profile_values(phase,coords)
        np.testing.assert_allclose(sample.T,baseline(coords),
             rtol=2e-10,atol=1e-10)


@pytest.mark.slow
def test_source_exact_native_root_remains_identical_after_fortran_pchip(physical):
    model,bg,opt,conversion,br,ref=physical
    solver=NativeSharpCowlingDVODE(bg,opt,samples=4097,
           coefficient_backend="native_source_cs2",
           source_restart_relative_jump=.01)
    assert all(p["pchip"] is None for p in solver._phases)
    m=solver.solve_mode(br,conversion=conversion,
           profile_points=250,root_backend="fortran_source_exact")
    assert abs(m["frequency_hz"]-ref)<.005
    assert (m["W_nodes_per_phase"],m["V_nodes_per_phase"])==(1,1)
    assert m["relative_join_mismatch"]<1e-5
    assert m["source_callback_rhs_calls"]==0
    assert m["root_fortran_restart_count"]==solver.source_restart_count


def test_native_solver_constructs_without_python_pchip(physical,monkeypatch):
    model,bg,opt,conversion,br,ref=physical
    def disallowed(*args,**kw):
        raise AssertionError("Python SciPy PCHIP invoked in source-exact native mode")
    monkeypatch.setattr(sharp,"PchipInterpolator",disallowed)
    solver=NativeSharpCowlingDVODE(bg,opt,samples=1025,
           coefficient_backend="native_source_cs2",
           source_restart_relative_jump=.01)
    assert all(p["pchip"] is None for p in solver._phases)
    assert solver.source_restart_count in (39,40)


def test_pchip_endpoint_flat_to_decreasing_is_exactly_zero():
    """The endpoint limiter must treat sign(0) as 0, not Fortran SIGN(+1,0)."""
    s=np.ascontiguousarray([0.,.1,.2,.5])
    y=np.asfortranarray(np.vstack([
      [1.,1.,0.,-1.],
      [0.,0.,-1.,-2.],
      [3.,3.,3.,2.],
      [1.,1.,1.,1.],
      [2.,2.,1.,1.]]).astype(float))
    expected=PchipInterpolator(s,y.T,axis=0,extrapolate=False)(s,nu=1).T
    v=np.full_like(y,-42.,order="F");status=ct.c_int(-1)
    native._native_library().cowling_native_pchip_slopes(
        len(s),s,y,v,ct.byref(status))
    assert status.value==0
    np.testing.assert_allclose(v,expected,rtol=1e-12,atol=1e-12)
    assert np.all(v[:,0]==0.)
