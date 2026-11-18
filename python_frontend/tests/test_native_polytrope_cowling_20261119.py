"""Independent analytic Gamma=2 polytropic EOS for native TOV/Cowling.

This is a numerical control, not a claim about nuclear microphysics.
With P=K n^2 and epsilon=m_b n + K n^2, the first law gives
h=ln(1+2Kn/m_b), cs2/c2 = 2Kn/(m_b+2Kn).
"""
import numpy as np
import pytest

from compose_eos.stellar_background import solve_background
from compose_eos.stellar_cowling import CowlingSettings, solve_mode
from compose_eos.stellar_cowling_native import NativeCowlingDVODE


class AnalyticGamma2Polytrope:
    mb = 939.
    K = 400.
    def enthalpy_at_density(self, n):
        if not 0<=n<=1.2:
            raise ValueError("polytrope density outside test support")
        return float(np.log1p(2*self.K*n/self.mb))
    def state_at_enthalpy(self, h):
        if not 0<=h<=self.enthalpy_at_density(1.2):
            raise ValueError("polytrope enthalpy outside test support")
        n=self.mb*np.expm1(h)/(2*self.K)
        p=self.K*n*n
        e=self.mb*n+p
        cs=2*self.K*n/(self.mb+2*self.K*n)
        return dict(n=float(n),p=float(p),epsilon=float(e),cs2_over_c2=float(cs))


@pytest.mark.parametrize("nc", [0.4,0.6])
def test_exact_polytrope_native_tov_mass_radius(nc):
    eos=AnalyticGamma2Polytrope()
    controls=dict(surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002)
    a=solve_background(eos,nc,**controls)
    b=solve_background(eos,nc,backend="fortran_dvode",native_samples=8193,**controls)
    assert abs(a.mass_msun-b.mass_msun)<2e-6
    assert abs(a.radius_km-b.radius_km)<2e-6


def test_exact_polytrope_native_cowling_f_p1_parity():
    eos=AnalyticGamma2Polytrope()
    control=dict(surface_h=1e-4,rtol=1e-10,atol=1e-12,max_step_h=.002)
    py=solve_background(eos,.4,**control)
    native=solve_background(eos,.4,backend="fortran_dvode",
                            native_samples=8193,**control)
    settings=CowlingSettings(surface_h=1e-4,rtol=1e-9,atol=1e-12)
    for name,br,nodes in (("f",[1594.59,1693.7],0),
                          ("p1",[3576.57,3675.68],1)):
        reference=solve_mode(py,br,settings,profile_points=240)
        found=solve_mode(native,br,settings,profile_points=240,
                         backend="fortran_dvode",native_samples=8193)
        print(name,"python",reference["frequency_hz"],
              "fortran",found["frequency_hz"],flush=True)
        assert reference["W_nodes"]==found["W_nodes"]==nodes
        assert reference["V_nodes"]==found["V_nodes"]==nodes
        assert abs(reference["frequency_hz"]-found["frequency_hz"])<0.05
