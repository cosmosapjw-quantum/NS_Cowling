"""C5/C6: sharp single Maxwell interface, pure-Fortran pchip root evaluations.

Authoritative: original source-locked nested EOS contracts; slow/rapid
one-sided phase transfer unchanged. The legacy source-exact hybrid remains
a different coefficient backend, not silently replaced by PCHIP.
"""
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest

from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
    NativeSharpCowlingDVODE, solve_sharp_background_native)
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"


@pytest.fixture(scope="module")
def maxwell():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
    with zipfile.ZipFile(SOURCE) as z:
        with zipfile.ZipFile(io.BytesIO(
                z.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            raw=inner.read("FRESH_5x2_WHOLE_STARS/SPG_M2/maxwell/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    bg=solve_sharp_background_native(eos,.30385128926615257,
         surface_h=3e-8,rtol=3e-10,atol=3e-11,
         max_step_h=.004,native_samples=2049)
    st=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,
                       max_step_h=.002,root_xtol_hz=1e-6)
    return NativeSharpCowlingDVODE(bg,st,samples=8193,
                                    coefficient_backend="pchip")


@pytest.mark.parametrize("conversion,bracket,nodes",[
    ("slow",(1800.,1900.),(0,0)),
    ("rapid",(1800.,1900.),(0,0)),
    ("slow",(5100.,5200.),(1,1)),
    ("rapid",(5110.,5170.),(1,1)),
])
def test_phasewise_fortran_matches_legacy_brent_source(maxwell,conversion,bracket,nodes):
    reference=maxwell.solve_mode(bracket,conversion=conversion,
        profile_points=200,root_backend="python_brentq")
    actual=maxwell.solve_mode(bracket,conversion=conversion,
        profile_points=200,root_backend="fortran_phase_pchip")
    assert abs(actual["frequency_hz"]-reference["frequency_hz"])<.005
    assert actual["root_backend"]=="fortran_phase_pchip"
    assert actual["root_converged"] is True
    assert actual["root_iterations"]>0
    assert actual["coefficient_backend"]=="pchip"
    assert actual["source_callback_rhs_calls"]==0
    assert (actual["W_nodes_per_phase"],actual["V_nodes_per_phase"])==nodes
    assert abs(actual["matching_residual"])<1e-5
    assert actual["relative_join_mismatch"]<1e-5
    assert all(j["traction_relative_residual"]<1e-8
               for j in actual["interface_junctions"])


def test_fortran_root_does_not_call_python_shooting_during_polish(maxwell,monkeypatch):
    def forbid(*args,**kw):
        raise AssertionError("No Python phasewise shooting is allowed while root polishing")
    monkeypatch.setattr(maxwell,"matching_residual",forbid)
    root,iterations,calls=maxwell._native_root_phase_pchip(
        (5110.,5170.),conversion="rapid")
    assert abs(root-5143.705574557526)<.005
    assert iterations>0 and calls>0


def test_unbracketed_native_sharp_root_rejected_without_fabricating_frequency(maxwell):
    with pytest.raises((RuntimeError,ValueError),match="bracket"):
        maxwell._native_root_phase_pchip((5600.,5700.),conversion="rapid")
    with pytest.raises(ValueError,match="bracket"):
        maxwell._native_root_phase_pchip((5000.,5000.),conversion="rapid")


def test_native_phase_root_rejects_unknown_conversion(maxwell):
    with pytest.raises(ValueError,match="conversion"):
        maxwell._native_root_phase_pchip((5110.,5170.),conversion="unknown")


def test_hybrid_cs2_must_not_silently_relabel_original_source_as_pchip(maxwell):
    hybrid=NativeSharpCowlingDVODE(maxwell.background,maxwell.settings,
                       samples=129,coefficient_backend="hybrid_cs2")
    with pytest.raises(ValueError,match="pchip"):
        hybrid.solve_mode((5110.,5170.),conversion="rapid",
           profile_points=160,root_backend="fortran_phase_pchip")


def test_phasewise_root_restores_smooth_native_setup_identity(maxwell):
    # The fresh Fortran phase setup must not poison the cached smooth table.
    from compose_eos.stellar_background import ContractEnthalpyEOS
    from compose_eos.stellar_cowling_native import NativeCowlingDVODE
    from compose_eos.stellar_cowling_tov_native import solve_background_native
    smooth=ContractEnthalpyEOS.from_path(
        ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json")
    bg=solve_background_native(smooth,.54,surface_h=1e-4,
        rtol=1e-10,atol=1e-12,max_step_h=.002,native_samples=2049)
    native=NativeCowlingDVODE(bg,CowlingSettings(surface_h=1e-4,
        rtol=1e-9,atol=1e-12,max_step_h=.002),samples=1025)
    before=native.matching_residual(2410.)
    maxwell._native_root_phase_pchip((5110.,5170.),conversion="rapid")
    after=native.matching_residual(2410.)
    assert after==before


@pytest.mark.slow
def test_css_slow_pchip_remains_numerically_unresolved_not_admitted():
    with zipfile.ZipFile(SOURCE) as z:
        with zipfile.ZipFile(io.BytesIO(
                z.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
            raw=inner.read("FRESH_5x2_WHOLE_STARS/SPG_M3/css/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    bg=solve_sharp_background_native(eos,.192510767556007,
         surface_h=3e-8,rtol=3e-10,atol=3e-11,
         max_step_h=.004,native_samples=2049)
    st=CowlingSettings(surface_h=3e-8,rtol=1e-9,
        atol=1e-12,max_step_h=.002)
    solver=NativeSharpCowlingDVODE(bg,st,samples=4097,
        coefficient_backend="pchip")
    with pytest.raises(RuntimeError,match="not physically admitted"):
        solver.solve_mode((3710.,3760.),conversion="slow",
            profile_points=350,root_backend="fortran_phase_pchip")


def test_native_abi_rejects_stale_source_callback_before_replacing_tables(maxwell):
    """A live exact-source callback must not silently become PCHIP physics."""
    import ctypes as ct
    from compose_eos import stellar_cowling_native as native
    call_s=[]
    callback_type=ct.CFUNCTYPE(
        None,ct.c_double,ct.POINTER(ct.c_double),ct.POINTER(ct.c_int))

    @callback_type
    def source_callback(s,coeff,status):
        call_s.append(float(s))
        status[0]=1

    high,low=maxwell._phases
    ht=maxwell.background.interfaces[0]["h"]
    state_low=maxwell.background.state_at_enthalpy(ht,"low")
    state_high=maxwell.background.state_at_enthalpy(ht,"high")
    iface=np.ascontiguousarray(
        [state_low["r_km"],state_low["Phi"],state_low["lambda_metric"],
         state_low["dPhi_dr_km"],state_low["p"],
         state_low["epsilon"],state_high["epsilon"]],dtype=np.float64)
    max_ds=min(maxwell._smax/200,
               maxwell.settings.max_step_h/(2*maxwell._smax))
    with native._LOCK:
        dll=native._native_library()
        register=ct.c_int(-1)
        dll.cowling_native_set_source_callback(
            ct.cast(source_callback,ct.c_void_p),ct.byref(register))
        assert register.value==0
        try:
            root,resid=ct.c_double(1.),ct.c_double(1.)
            iters,nfev,status=ct.c_int(-1),ct.c_int(-1),ct.c_int(-1)
            dll.cowling_native_sharp_onejump_brent(
                5110.,5170.,2,maxwell.settings.rtol,maxwell.settings.atol,
                max_ds,10,
                len(high["s"]),high["s"],high["vals"],high["slopes"],
                len(low["s"]),low["s"],low["vals"],low["slopes"],
                0.5*float(high["s"][-1]),iface,2,
                maxwell.settings.root_xtol_hz,ct.byref(root),ct.byref(resid),
                ct.byref(iters),ct.byref(nfev),ct.byref(status))
            assert status.value==2
            assert root.value==0. and resid.value==0.
            assert iters.value==0 and nfev.value==0
            assert not call_s
        finally:
            dll.cowling_native_set_source_callback(None,ct.byref(register))
            native._ACTIVE=None
            assert register.value==0
