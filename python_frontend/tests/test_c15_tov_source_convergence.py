"""C15 original source-locked two-phase native Fortran TOV control batch.

TDD contract: one source/EOS phase sampling, Fortran PCHIP of P and epsilon,
all TOV profiles and independent per-run DVODE histories in Fortran.
Python only packages the source states and reports resulting stars.
"""
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pytest
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import solve_sharp_background_native
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
CASES=(
 ("SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.)),
 ("SPG_M3","css",.192510767556007,"slow",(3710.,3760.))
)

@pytest.fixture(scope="module",params=CASES)
def case(request):
    model,kind,nc,conversion,bracket=request.param
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SHA
    with zipfile.ZipFile(SOURCE) as full:
      with zipfile.ZipFile(io.BytesIO(full.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as nested:
        raw=nested.read(f"FRESH_5x2_WHOLE_STARS/{model}/{kind}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    return model,eos,nc,conversion,bracket

def _sweep(eos,nc,controls,samples=2049):
    from compose_eos.stellar_sharp_native import solve_sharp_background_native_sweep
    return solve_sharp_background_native_sweep(eos,nc,controls,
       native_samples=samples,surface_h=3e-8)

@pytest.mark.parametrize("method",["adams","bdf"])
def test_tov_sweep_source_exact_matches_independently_rebuilt_legacy_method(case,method):
    model,eos,nc,conversion,bracket=case
    cfg=dict(rtol=3e-10,atol=3e-11,max_step_h=.004,method=method)
    baseline=solve_sharp_background_native(eos,nc,
        surface_h=3e-8,native_samples=2049,**cfg)
    batch=_sweep(eos,nc,[cfg])
    assert len(batch)==1 and batch[0]["status"]==0
    star=batch[0]["star"]
    assert star.config["backend"]=="Fortran DVODE native two-phase TOV sweep"
    assert star.config["tov_eos_slopes_backend"]=="fortran_pchip"
    assert star.mass_msun==pytest.approx(baseline.mass_msun,rel=2e-7)
    assert star.radius_km==pytest.approx(baseline.radius_km,rel=2e-7)
    assert star.mass_geom_km>0 and star.radius_km>2*star.mass_geom_km
    assert len(star.interfaces)==1
    assert star.interfaces[0]["r_jump_km"]==star.interfaces[0]["m_jump_km"]==0.
    assert batch[0]["tov_native_nfev"]>0
    for h in (eos.transitions[0]["h"]*.9,eos.transitions[0]["h"]*1.1):
        if h>=star.h_start:continue
        side="low" if h<eos.transitions[0]["h"] else "high"
        a=star.state_at_enthalpy(h,side)
        b=baseline.state_at_enthalpy(h,side)
        for key in ("r_km","m_geom_km","Phi","lambda_metric","dPhi_dr_km"):
            assert a[key]==pytest.approx(b[key],rel=1e-6,abs=1e-8),(model,side,key)

@pytest.mark.slow
def test_tov_sweep_original_p1_branch_identity(case):
    model,eos,nc,conversion,bracket=case
    cfg=dict(rtol=3e-10,atol=3e-11,max_step_h=.004,method="adams")
    b=_sweep(eos,nc,[cfg])
    assert b[0]["status"]==0
    opts=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002,root_xtol_hz=1e-5)
    solver=NativeSharpCowlingDVODE(b[0]["star"],opts,samples=4097,
         coefficient_backend="native_source_cs2",source_restart_relative_jump=.01)
    mode=solver.solve_mode(bracket,conversion=conversion,
         profile_points=300,root_backend="fortran_source_exact")
    ref=5143.719711814548 if model=="SPG_M2" else 3736.210706737009
    assert abs(mode["frequency_hz"]-ref)<.02
    assert (mode["W_nodes_per_phase"],mode["V_nodes_per_phase"])==(1,1)
    assert mode["relative_join_mismatch"]<1e-5
    assert mode["source_callback_rhs_calls"]==0
    assert mode["root_fortran_restart_count"]==solver.source_restart_count


def test_tov_sweep_bad_control_is_isolated_and_zeroed(case):
    model,eos,nc,conversion,bracket=case
    a=_sweep(eos,nc,[
        dict(rtol=3e-10,atol=3e-11,max_step_h=.004,method="adams"),
        dict(rtol=0.,atol=3e-11,max_step_h=.004,method="adams"),
        dict(rtol=1e-10,atol=3e-11,max_step_h=.002,method="bdf")])
    assert [r["status"] for r in a]==[0,1,0],(model,a)
    assert a[1]["star"] is None
    assert a[1]["tov_native_nfev"]==0
    assert all(a[j]["star"].mass_msun>0 for j in (0,2))


def test_tov_sweep_rejects_invalid_phase_input(case):
    model,eos,nc,conversion,bracket=case
    with pytest.raises(ValueError,match="controls"):
        _sweep(eos,nc,[])
    with pytest.raises(ValueError,match="method"):
        _sweep(eos,nc,[dict(rtol=1e-9,atol=1e-11,max_step_h=.005,method="invalid")])
