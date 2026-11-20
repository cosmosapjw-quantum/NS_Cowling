"""C12: Fortran-native source-exact sharp Maxwell roots and knot scheduling.

Original EOS bytes and one-sided scientific invariants are immutable.
A successful native root must use no Python bracket function evaluations
and no Python EOS RHS callback, preserve M2/M3 p1 node identity and
compare with the previously verified native-source SciPy Brent baseline.
"""
from __future__ import annotations
import ctypes as ct
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pytest

from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos import stellar_cowling_native as native

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
CASES=(
 ("SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.),5143.719721563306),
 ("SPG_M3","css",.192510767556007,"slow",(3710.,3760.),3736.2088951581254),
)

@pytest.fixture(scope="module",params=CASES)
def case(request):
    model,join,nc,conversion,bracket,reference=request.param
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
    with zipfile.ZipFile(SOURCE) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            raw=inner.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    bg=solve_sharp_background_native(eos,nc,surface_h=3e-8,rtol=3e-10,
       atol=3e-11,max_step_h=.004,native_samples=2049)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
          atol=1e-12,max_step_h=.002,root_xtol_hz=1e-5)
    solver=NativeSharpCowlingDVODE(bg,settings,samples=4097,
          coefficient_backend="native_source_cs2",source_restart_relative_jump=.01)
    return model,join,solver,conversion,bracket,reference

@pytest.mark.slow
def test_native_source_brent_reproduces_sharp_p1_and_original_reference(case):
    model,join,solver,conversion,br,external_reference=case
    baseline=solver.solve_mode(br,conversion=conversion,
                root_backend="python_brentq",profile_points=300)
    result=solver.solve_mode(br,conversion=conversion,
                root_backend="fortran_source_exact",profile_points=300)
    assert result["root_backend"]=="fortran_source_exact"
    assert result["root_native_rhs_calls"]>0
    assert result["root_iterations"]>0
    assert result["root_converged"] is True
    assert abs(result["frequency_hz"]-baseline["frequency_hz"])<.005
    assert abs(result["frequency_hz"]-external_reference)<.05
    assert (result["W_nodes_per_phase"],result["V_nodes_per_phase"])==(1,1)
    assert result["relative_join_mismatch"]<1e-5
    assert abs(result["matching_residual"])<1e-5
    assert result["source_callback_rhs_calls"]==0
    assert result["root_fortran_restart_count"]==solver.source_restart_count


@pytest.mark.slow
def test_native_fortran_brent_does_not_invoke_python_residual(case,monkeypatch):
    model,join,solver,conversion,br,external_reference=case
    def raise_on_python_brent(*args,**kw):
        raise AssertionError("Python residual callback was invoked during Fortran root")
    monkeypatch.setattr(solver,"matching_residual",raise_on_python_brent)
    root,iterations,nfev,cutcount=solver._native_source_exact_root(
             br,conversion=conversion)
    assert root==pytest.approx(external_reference,abs=.05)
    assert iterations>0 and nfev>0
    assert cutcount==solver.source_restart_count


def test_fortran_knot_selector_matches_original_strong_jumps(case):
    model,join,solver,conversion,br,reference=case
    with native._LOCK:
        dll=native._native_library()
        solver._native_exact_eos_setup(dll)
        for phase in solver._phases:
            expected=phase["source_restart_s"]
            capacity=solver._native_source_eos[0].shape[1]+solver._native_source_eos[1].shape[1]+5
            cuts=np.zeros(capacity,dtype=np.float64)
            count,status=ct.c_int(-1),ct.c_int(-1)
            dll.cowling_native_source_select_restarts(
              capacity,solver.background.h_start,
              phase["h0"],phase["h1"],.01,
              2 if phase is solver._phases[0] else 1,
              cuts,ct.byref(count),ct.byref(status))
            assert status.value==0
            got=cuts[:count.value]
            assert len(got)==len(expected),(model,len(got),len(expected))
            np.testing.assert_allclose(got,expected,rtol=0,atol=2e-12)


def test_root_rejects_invalid_onejump_input(case):
    model,join,solver,conversion,br,reference=case
    with pytest.raises(ValueError,match="bracket"):
        solver._native_source_exact_root((br[0],br[0]),conversion=conversion)
    with pytest.raises(ValueError,match="conversion"):
        solver._native_source_exact_root(br,conversion="invalid")
    # A fixed interval may contain a legitimate root for a different EOS.
    # Verify the premise independently with the original scalar determinant.
    intervals=[(br[0]+150.,br[0]+175.),(br[1]+175.,br[1]+200.),
               (5600.,5700.),(2000.,2100.)]
    nonbracket=next(((a,b) for a,b in intervals
                     if solver.matching_residual(a,conversion=conversion)*
                        solver.matching_residual(b,conversion=conversion)>0),None)
    assert nonbracket is not None,"test input must be source-verified nonbracket"
    with pytest.raises((RuntimeError,ValueError),match="bracket"):
        solver._native_source_exact_root(nonbracket,conversion=conversion)


def test_new_native_root_refuses_non_source_backend(case):
    model,join,solver,conversion,br,reference=case
    other=NativeSharpCowlingDVODE(
        solver.background,solver.settings,samples=129,
        coefficient_backend="pchip")
    with pytest.raises(ValueError,match="native_source_cs2"):
        other.solve_mode(br,conversion=conversion,profile_points=100,
                         root_backend="fortran_source_exact")


def test_new_native_root_c_abi_rejects_active_python_source_callback(case):
    model,join,solver,conversion,br,reference=case
    prototype=ct.CFUNCTYPE(None,ct.c_double,ct.POINTER(ct.c_double),
                              ct.POINTER(ct.c_int))
    called=[]

    @prototype
    def callback(s,values,status):
        called.append(float(s))
        status[0]=1

    high,low=solver._phases
    ht=solver.background.interfaces[0]["h"]
    state_low=solver.background.state_at_enthalpy(ht,"low")
    state_high=solver.background.state_at_enthalpy(ht,"high")
    geom=np.ascontiguousarray([
       state_low["r_km"],state_low["Phi"],state_low["lambda_metric"],
       state_low["dPhi_dr_km"],state_low["p"],state_low["epsilon"],
       state_high["epsilon"]],dtype=np.float64)
    max_ds=min(solver._smax/200,
               solver.settings.max_step_h/(2*solver._smax))
    with native._LOCK:
        dll=native._native_library()
        solver._native_exact_eos_setup(dll)
        st=ct.c_int(-1)
        dll.cowling_native_set_source_callback(
            ct.cast(callback,ct.c_void_p),ct.byref(st))
        assert st.value==0
        try:
            root,det=ct.c_double(-1),ct.c_double(-1)
            iterations,nfev,cuts,status=(ct.c_int(-1),ct.c_int(-1),
                                         ct.c_int(-1),ct.c_int(-1))
            dll.cowling_native_sharp_source_brent(
                br[0],br[1],solver.settings.ell,
                solver.settings.rtol,solver.settings.atol,max_ds,
                {"adams":10,"bdf":22}[solver.method],
                len(high["s"]),high["s"],high["vals"],high["slopes"],
                len(low["s"]),low["s"],low["vals"],low["slopes"],
                solver.background.h_start,0.5*high["s"][-1],geom,
                {"slow":1,"rapid":2}[conversion],
                solver.settings.root_xtol_hz,solver.source_restart_relative_jump,
                ct.byref(root),ct.byref(det),ct.byref(iterations),
                ct.byref(nfev),ct.byref(cuts),ct.byref(status))
            assert status.value==2
            assert (root.value,det.value)==(0.,0.)
            assert (iterations.value,nfev.value,cuts.value)==(0,0,0)
            assert called==[]
        finally:
            dll.cowling_native_set_source_callback(None,ct.byref(st))
            native._ACTIVE=None
            assert st.value==0
