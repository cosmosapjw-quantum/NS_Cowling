"""C14: batched Fortran source-exact original-star convergence controls.

TDD: native batch ABI must compute different DVODE tolerances and methods
without Python root/RHS callbacks, preserve independent failing controls,
and leave the original C13 production root path unchanged.
"""
from __future__ import annotations
import ctypes as ct
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pytest
from compose_eos import stellar_cowling_native as native
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT=Path(__file__).resolve().parents[2]
ZIP=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
MODELS=[
 ("SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.),5143.719634196357),
 ("SPG_M3","css",.192510767556007,"slow",(3710.,3760.),3736.210436364205)]

@pytest.fixture(scope="module",params=MODELS)
def source_case(request):
    model,join,nc,conv,br,ref=request.param
    assert hashlib.sha256(ZIP.read_bytes()).hexdigest()==SHA
    with zipfile.ZipFile(ZIP) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inside:
            blob=inside.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(blob)))
    bg=solve_sharp_background_native(eos,nc,surface_h=3e-8,rtol=3e-10,atol=3e-11,
             max_step_h=.004,native_samples=2049)
    opts=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,max_step_h=.002)
    sol=NativeSharpCowlingDVODE(bg,opts,samples=8193,coefficient_backend="native_source_cs2",
              source_restart_relative_jump=.01)
    return model,sol,conv,br,ref

def native_batch(sol,conv,br,rtols,steps,methods):
    hi,lo=sol._phases
    h=sol.background.interfaces[0]["h"]
    low=sol.background.state_at_enthalpy(h,"low")
    high=sol.background.state_at_enthalpy(h,"high")
    geom=np.ascontiguousarray([
       low["r_km"],low["Phi"],low["lambda_metric"],low["dPhi_dr_km"],
       low["p"],low["epsilon"],high["epsilon"]],dtype=np.float64)
    rtol=np.ascontiguousarray(rtols,dtype=np.float64)
    max_ds=np.ascontiguousarray(steps,dtype=np.float64)
    method=np.ascontiguousarray(methods,dtype=np.int32)
    n=len(rtol)
    assert n==len(max_ds)==len(method)
    root=np.full(n,-1.,dtype=np.float64)
    det=np.full(n,-1.,dtype=np.float64)
    iters=np.full(n,-1,dtype=np.int32)
    calls=np.full(n,-1,dtype=np.int32)
    knots=np.full(n,-1,dtype=np.int32)
    statuses=np.full(n,-1,dtype=np.int32)
    with native._LOCK:
        dll=native._native_library()
        sol._native_exact_eos_setup(dll)
        try:
            dll.cowling_native_sharp_source_sweep(
                n,br[0],br[1],sol.settings.ell,sol.settings.atol,
                sol.background.h_start,0.5*float(hi["s"][-1]),geom,
                {"slow":1,"rapid":2}[conv],
                sol.settings.root_xtol_hz,.01,
                len(hi["s"]),hi["s"],hi["vals"],hi["slopes"],
                len(lo["s"]),lo["s"],lo["vals"],lo["slopes"],
                rtol,max_ds,method,
                root,det,iters,calls,knots,statuses)
        finally:
            native._ACTIVE=None
    return root,det,iters,calls,knots,statuses


@pytest.mark.slow
def test_source_exact_native_batch_matches_single_shots_and_preserves_failures(source_case,monkeypatch):
    model,sol,conv,br,ref=source_case
    # The 0 relative tolerance is invalid, and the subsequent case
    # must nevertheless complete independently.
    base_ds=min(sol._smax/200,sol.settings.max_step_h/(2*sol._smax))
    def forbid(*args,**kwargs):
        raise AssertionError("Python root residual was called by native batch")
    monkeypatch.setattr(sol,"matching_residual",forbid)
    r,d,it,nf,kn,st=native_batch(sol,conv,br,
            [1e-9,0.,1e-9],
            [base_ds,base_ds,base_ds],
            [10,10,22])
    assert list(st)==[0,1,0],(model,st)
    assert r[1]==d[1]==0. and it[1]==nf[1]==kn[1]==0
    assert all(abs(r[j]-ref)<.01 for j in (0,2))
    assert all(abs(d[j])<1e-5 for j in (0,2))
    assert all(it[j]>0 and nf[j]>0 for j in (0,2))
    assert kn[0]==kn[2]==sol.source_restart_count


def test_source_exact_native_batch_global_invalid_bracket_zeroes_outputs(source_case):
    model,sol,conv,br,ref=source_case
    ds=min(sol._smax/200,sol.settings.max_step_h/(2*sol._smax))
    r,d,it,nf,kn,st=native_batch(sol,conv,(br[0],br[0]),
                      [1e-9],[ds],[10])
    assert st[0]!=0
    assert r[0]==d[0]==it[0]==nf[0]==kn[0]==0


@pytest.mark.slow
def test_public_frontend_sweep_delegates_entire_batch_to_fortran(source_case,monkeypatch):
    model,sol,conversion,br,ref=source_case
    ds=min(sol._smax/200,sol.settings.max_step_h/(2*sol._smax))
    controls=[
      dict(rtol=1e-9,max_step_s=ds,method="adams"),
      dict(rtol=0.,max_step_s=ds,method="adams"),
      dict(rtol=1e-9,max_step_s=ds,method="bdf")]
    def forbid(*args,**kw):
        raise AssertionError("Python matching_residual invoked by numerical batch")
    monkeypatch.setattr(sol,"matching_residual",forbid)
    results=sol.sweep_native_numerical_controls(
      br,conversion=conversion,controls=controls)
    assert len(results)==3
    assert [r["status"] for r in results]==[0,1,0]
    assert results[1]["frequency_hz"] is None
    assert all(abs(results[j]["frequency_hz"]-ref)<.01 for j in (0,2))
    assert results[0]["source_restart_count"]==sol.source_restart_count
    assert results[2]["source_restart_count"]==sol.source_restart_count
    assert [r["method"] for r in results]==["adams","adams","bdf"]
