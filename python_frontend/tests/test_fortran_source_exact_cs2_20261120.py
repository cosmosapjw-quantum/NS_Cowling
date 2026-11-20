"""Source-locked direct Fortran c_s² physical-law tests (no PCHIP c_s²).

A red-green implementation gate for a Fortran-only source EOS evaluator:
surface, pressure-primary outer polytropes, causal exponential bridge,
upper cubic-P(mu) Hermite intervals, and one-sided Maxwell endpoint.
"""
from __future__ import annotations
import ctypes as ct
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pytest
from compose_eos import stellar_cowling_native as native
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"

def load_eos(model,join):
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
    with zipfile.ZipFile(SOURCE) as outer:
        with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            raw=inner.read(f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    return eos,hashlib.sha256(raw).hexdigest()

def arrays(eos):
    b=eos.barotrope
    outer=np.asfortranarray(np.vstack([
        b.outer.nb,b.outer.p,b.outer.F,
        np.r_[b.outer.gamma,b.outer.gamma[-1]]]))
    upper=np.asfortranarray(np.vstack([
        b.upper_n,b.upper_p,b.upper_epsilon]))
    surface=np.ascontiguousarray([b.surface.rest_energy_MeV,
              b.surface.n0,b.surface.p0,b.surface.gamma])
    cb=b.bridge
    bridge=np.ascontiguousarray([cb.n_left,cb.p_left,cb.epsilon_left,
               cb.n_right,cb.p_right,cb.epsilon_right,cb.kappa])
    return outer,upper,surface,bridge,int(sorted(b.sharp_interfaces)[0])

def setup(eos):
    arr=arrays(eos)
    outer,upper,surface,bridge,jump=arr
    dll=native._native_library()
    status=ct.c_int(-1)
    dll.cowling_native_source_eos_setup(
        outer.shape[1],outer,upper.shape[1],upper,surface,bridge,
        jump,ct.byref(status))
    assert status.value==0
    return dll

def evaluate(dll,h,side="high"):
    value=ct.c_double(0.)
    region,status=ct.c_int(-1),ct.c_int(-1)
    dll.cowling_native_source_cs2_eval(
        float(h),{"low":1,"high":2}[side],
        ct.byref(value),ct.byref(region),ct.byref(status))
    return value.value,region.value,status.value

@pytest.mark.parametrize("model,join",[
    ("SPG_M2","maxwell"),("SPG_M3","css")])
def test_fortran_exact_cs2_matches_original_source_on_four_analytic_regions(model,join):
    eos,source_hash=load_eos(model,join)
    dll=setup(eos)
    edges=[0.,eos.h_surface_anchor,eos.h_outer_end,
           eos.h_bridge_end,eos.transitions[0]["h"],eos.h_max]
    seen=set()
    for region,(a,b) in enumerate(zip(edges[:-1],edges[1:]),start=1):
        # Midpoint and interior points, not derivative-discontinuous endpoints.
        for h in np.linspace(a,b,26)[1:-1]:
            speed,kind,status=evaluate(dll,h)
            assert status==0,(source_hash,region,h,status)
            assert 0<speed<=1.
            if region<=3:
                assert kind==region,(model,region,h,kind)
            else:
                assert kind==4
            ref=eos.state_at_enthalpy(float(h),"high")["cs2_over_c2"]
            np.testing.assert_allclose(speed,ref,rtol=8e-10,atol=2e-13,
                                       err_msg=f"{model}:{region}:{h}")
            seen.add(kind)
    assert seen=={1,2,3,4}


@pytest.mark.parametrize("model,join",[
    ("SPG_M2","maxwell"),("SPG_M3","css")])
def test_fortran_source_cs2_preserves_maxwell_one_sided_limits(model,join):
    eos,_=load_eos(model,join)
    dll=setup(eos)
    ht=eos.transitions[0]["h"]
    low,region_lo,stat_lo=evaluate(dll,ht,side="low")
    high,region_hi,stat_hi=evaluate(dll,ht,side="high")
    assert (stat_lo,stat_hi)==(0,0)
    assert region_lo==region_hi==4
    assert abs(high-low)>1e-4
    assert low==pytest.approx(eos.state_at_enthalpy(ht,"low")["cs2_over_c2"],rel=1e-10)
    assert high==pytest.approx(eos.state_at_enthalpy(ht,"high")["cs2_over_c2"],rel=1e-10)


@pytest.mark.parametrize("model,join",[
    ("SPG_M2","maxwell"),("SPG_M3","css")])
def test_fortran_source_cs2_one_sided_near_outer_and_upper_knots(model,join):
    eos,_=load_eos(model,join)
    dll=setup(eos)
    for h in eos.break_h[::max(1,len(eos.break_h)//18)]:
        if h<=1e-7 or h>=eos.h_max-1e-7:
            continue
        for shift in (-3e-10,3e-10):
            q=float(h+shift)
            if not 0<q<eos.h_max or eos._transition(q) is not None:
                continue
            got,_,status=evaluate(dll,q)
            ref=eos.state_at_enthalpy(q)["cs2_over_c2"]
            assert status==0,(model,q,status)
            np.testing.assert_allclose(got,ref,rtol=1e-8,atol=2e-12)


def test_fortran_exact_source_evaluator_rejects_invalid_argument_and_zeroes_result():
    eos,_=load_eos("SPG_M3","css")
    dll=setup(eos)
    for h,side in [(float("nan"),2),(-1.,2),(eos.h_max+.1,2),(.1,3)]:
        value=ct.c_double(-4.)
        kind,status=ct.c_int(-1),ct.c_int(-1)
        dll.cowling_native_source_cs2_eval(h,side,
               ct.byref(value),ct.byref(kind),ct.byref(status))
        assert status.value!=0
        assert value.value==0.
        assert kind.value==0


def test_fortran_source_setup_rejects_unphysical_upper_segment_and_invalidates_old_state():
    eos,_=load_eos("SPG_M2","maxwell")
    dll=setup(eos)
    outer,upper,surface,bridge,jump=arrays(eos)
    damaged=upper.copy(order="F")
    damaged[1,5]=damaged[1,4]*.90
    status=ct.c_int(-1)
    dll.cowling_native_source_eos_setup(
       outer.shape[1],outer,damaged.shape[1],damaged,
       surface,bridge,jump,ct.byref(status))
    assert status.value!=0
    speed,kind,csstatus=evaluate(dll,eos.transitions[0]["h"]+.001)
    assert csstatus!=0
    assert speed==0 and kind==0
