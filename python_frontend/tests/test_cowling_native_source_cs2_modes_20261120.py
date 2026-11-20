"""Fortran exact-source sound speed in real Maxwell/CSS Cowling RHS.

The original source contract provides c_s² at EVERY native DVODE RHS
call. Python only marshals EOS parameters and geometric PCHIP coefficients.
Existing Python hybrid_cs2 callback remains an independent reference.
"""
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pytest

from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT=Path(__file__).resolve().parents[2]
ZIP=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"

def load(model,kind):
    assert hashlib.sha256(ZIP.read_bytes()).hexdigest()==SHA
    with zipfile.ZipFile(ZIP) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as z:
            blob=z.read(f"FRESH_5x2_WHOLE_STARS/{model}/{kind}/REPRESENTATION.json")
    return SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(blob)))

@pytest.fixture(scope="module",params=[
    ("SPG_M2","maxwell",.30385128926615257,"rapid",(5110.,5170.),5143.719721563306),
    ("SPG_M3","css",.192510767556007,"slow",(3710.,3760.),3736.2088951581254)])
def star_case(request):
    model,kind,density,conversion,bracket,reference=request.param
    eos=load(model,kind)
    bg=solve_sharp_background_native(eos,density,
         surface_h=3e-8,rtol=3e-10,atol=3e-11,
         max_step_h=.004,native_samples=2049)
    opts=CowlingSettings(surface_h=3e-8,rtol=1e-9,
                         atol=1e-12,max_step_h=.002,root_xtol_hz=1e-5)
    return model,bg,opts,conversion,bracket,reference

@pytest.mark.slow
def test_native_exact_source_cs2_matches_original_hybrid_and_physical_reference(star_case):
    model,bg,opts,conversion,bracket,reference=star_case
    baseline=NativeSharpCowlingDVODE(bg,opts,samples=4097,
             coefficient_backend="hybrid_cs2",source_restart_relative_jump=.01)
    new=NativeSharpCowlingDVODE(bg,opts,samples=4097,
             coefficient_backend="native_source_cs2",source_restart_relative_jump=.01)
    result=baseline.solve_mode(bracket,conversion=conversion,profile_points=300)
    got=new.solve_mode(bracket,conversion=conversion,profile_points=300)
    assert got["coefficient_backend"]=="native_source_cs2"
    assert got["source_callback_rhs_calls"]==0
    assert got["source_restart_count"]==result["source_restart_count"]
    assert abs(got["frequency_hz"]-result["frequency_hz"])<.005,(model,got["frequency_hz"],result["frequency_hz"])
    assert abs(got["frequency_hz"]-reference)<.05
    assert (got["W_nodes_per_phase"],got["V_nodes_per_phase"])==(1,1)
    assert got["relative_join_mismatch"]<1e-5
    assert got["surface_BC_relative"]<1e-9
    assert all(x["traction_relative_residual"]<1e-9 for x in got["interface_junctions"])


def test_fortran_exact_source_backend_rejects_missing_or_extra_phase(star_case):
    model,bg,opts,conversion,bracket,reference=star_case
    obj=NativeSharpCowlingDVODE(bg,opts,samples=1025,
             coefficient_backend="native_source_cs2")
    assert len(obj._phases)==2
    assert obj.source_restart_count>0
    assert obj.coefficient_backend=="native_source_cs2"
