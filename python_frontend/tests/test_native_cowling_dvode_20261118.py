"""Native DVODE Cowling integration: physical parity, modes, and failure contracts."""
import math

import numpy as np
import pytest

from compose_eos.stellar_cowling import CowlingSettings, C_KM_S, solve_mode
from compose_eos.stellar_cowling_native import NativeCowlingDVODE


class UniformDensityControl:
    """Schwarzschild interior, incompressible Kelvin/Cowling limit (acausal test only)."""
    def __init__(self, compactness=1e-4, radius=10.0, r0=0.01):
        self.C, self.R = compactness, radius
        self.s = np.sqrt(1 - 2 * self.C)
        self.h_c = float(np.log(2 * self.s / (3 * self.s - 1)))
        self.h_start = self.enthalpy(r0)

    def enthalpy(self, radius):
        return float(np.log(2 * self.s /
                            (3 * self.s - np.sqrt(1 - 2*self.C*(radius/self.R)**2))))

    def state_at_enthalpy(self, h):
        z = self.s * (3 - 2 * np.exp(-h))
        r = self.R * np.sqrt((1 - z*z) / (2*self.C))
        ep = (3*self.s-z)/2
        return dict(r_km=float(r), Phi=float(np.log(ep)),
                    lambda_metric=float(-np.log(z)),
                    dPhi_dr_km=float(self.C*r/(self.R**2*z*ep)),
                    cs2_over_c2=np.inf)


def _settings():
    return CowlingSettings(surface_h=1e-10, rtol=1e-9, atol=1e-12,
                           max_step_h=1e-6, root_xtol_hz=1e-6)


@pytest.fixture(scope="module")
def reference_star():
    bg = UniformDensityControl()
    return bg, C_KM_S*np.sqrt(2*bg.C/bg.R**2)/(2*np.pi)


def test_native_backend_matches_weak_gravity_kelvin_mode(reference_star):
    bg, f0 = reference_star
    settings = _settings()
    native = NativeCowlingDVODE(bg, settings, samples=4097)
    result = native.solve_mode([0.99*f0, 1.01*f0], profile_points=400)
    exact = solve_mode(bg, [0.99*f0, 1.01*f0], settings, profile_points=400)
    assert abs(result["frequency_hz"]/f0 - 1.0) < 5e-4
    assert abs(result["frequency_hz"]-exact["frequency_hz"]) < 1e-2
    assert result["W_nodes"] == 0
    assert abs(result["matching_residual"]) < 2e-7
    assert result["relative_join_mismatch"] < 2e-7
    assert result["surface_BC_relative"] < 1e-10
    assert result["backend"] == "Fortran DVODE Adams"
    assert len(result["profile"]) == 399


def test_native_shot_rejects_nonphysical_inputs(reference_star):
    bg, f0 = reference_star
    native = NativeCowlingDVODE(bg, _settings(), samples=1025)
    for frequency in (0, -1, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="frequency"):
            native.matching_residual(frequency)
    with pytest.raises(ValueError, match="bracket"):
        native.solve_mode([f0, f0])


def test_native_adams_bdf_parity(reference_star):
    bg, f0 = reference_star
    adams = NativeCowlingDVODE(bg, _settings(), samples=4097, method="adams")
    bdf = NativeCowlingDVODE(bg, _settings(), samples=4097, method="bdf")
    for factor in (0.95, 1.05):
        a = adams.matching_residual(factor*f0)
        b = bdf.matching_residual(factor*f0)
        assert abs(a-b) < 2e-4


def test_native_sampling_convergence(reference_star):
    bg, f0 = reference_star
    coarse = NativeCowlingDVODE(bg, _settings(), samples=2049)
    fine = NativeCowlingDVODE(bg, _settings(), samples=8193)
    f_coarse = coarse.solve_mode([0.99*f0,1.01*f0], profile_points=160)["frequency_hz"]
    f_fine = fine.solve_mode([0.99*f0,1.01*f0], profile_points=160)["frequency_hz"]
    assert abs(f_coarse-f_fine) < 2e-2


@pytest.fixture(scope="module")
def canonical_smooth_star():
    """Original direct-contract EOS; exact Python DOP853 results computed at 1e-9."""
    import json
    from pathlib import Path
    from compose_eos.stellar_background import ContractEnthalpyEOS, solve_background
    root=Path(__file__).resolve().parents[2]
    source=root/"data/core_support_20261117/primary_v1/REPRESENTATION.json"
    mass=json.loads((root/"data/stellar_consumers_20261117/tov/CANONICAL_1P4.json").read_text())
    eos=ContractEnthalpyEOS.from_path(source)
    star=solve_background(eos,mass["central_density_fm3"],surface_h=1e-4,
                          rtol=1e-10,atol=1e-12,max_step_h=0.002)
    return star


def test_native_real_eos_f_and_p1_match_exact_contract(canonical_smooth_star):
    star=canonical_smooth_star
    settings=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)
    native=NativeCowlingDVODE(star,settings,samples=8193)
    # Source run on this branch at the same star and controls: Python DOP853,
    # full exact EOS laws at every RHS call, 2026-11-18.
    references=[("f",(2250,2500),2401.381874932175,0),
                ("p1",(6500,6750),6721.957056883109,1)]
    for label,br,frequency,nodes in references:
        found=native.solve_mode(br,profile_points=400)
        assert abs(found["frequency_hz"]-frequency)<0.01, (label,found["frequency_hz"],frequency)
        assert found["W_nodes"]==nodes, label
        assert found["V_nodes"]==nodes, label
        assert found["relative_join_mismatch"]<2e-6


def test_public_python_frontend_can_dispatch_fortran(reference_star):
    from compose_eos.stellar_cowling import solve_mode as frontend_solve
    from compose_eos.stellar_cowling import scan_brackets as frontend_scan
    bg,f0=reference_star
    result=frontend_solve(bg,[0.99*f0,1.01*f0],_settings(),180,
                          backend="fortran_dvode",native_samples=2049)
    assert result["backend"]=="Fortran DVODE Adams"
    assert abs(result["frequency_hz"]-f0)/f0<5e-4
    brackets,scan=frontend_scan(bg,[0.99*f0,1.01*f0],_settings(),
                               backend="fortran_dvode",native_samples=2049)
    assert len(brackets)==1 and len(scan)==2
    with pytest.raises(ValueError,match="backend"):
        frontend_solve(bg,[0.99*f0,1.01*f0],_settings(),backend="unknown")


def test_real_eos_table_resolution_stable(canonical_smooth_star):
    star=canonical_smooth_star
    settings=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)
    coarse=NativeCowlingDVODE(star,settings,samples=8193)
    fine=NativeCowlingDVODE(star,settings,samples=16385)
    for br in ([2250,2500],[6500,6750]):
        a=coarse.solve_mode(br,profile_points=240)
        b=fine.solve_mode(br,profile_points=240)
        assert abs(a["frequency_hz"]-b["frequency_hz"])<0.01
        assert a["W_nodes"]==b["W_nodes"]
        assert a["V_nodes"]==b["V_nodes"]


def test_sharp_interface_rejected_before_native_shooting():
    bg=UniformDensityControl()
    bg.sharp_interfaces=True
    with pytest.raises(ValueError,match="sharp Maxwell"):
        NativeCowlingDVODE(bg,_settings())


def test_native_real_star_adams_bdf_agree(canonical_smooth_star):
    bg=canonical_smooth_star
    st=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)
    adams=NativeCowlingDVODE(bg,st,method="adams",samples=8193)
    bdf=NativeCowlingDVODE(bg,st,method="bdf",samples=8193)
    for lo,hi in [(2250.,2500.),(6500.,6750.)]:
        a=adams.solve_mode([lo,hi],profile_points=300)
        b=bdf.solve_mode([lo,hi],profile_points=300)
        assert abs(a["frequency_hz"]-b["frequency_hz"])<0.01
        assert a["W_nodes"]==b["W_nodes"]
        assert a["V_nodes"]==b["V_nodes"]
