"""Boundary-anchored thermodynamic reconstruction of the published BSk22 fit.

This OPTIONAL representation preserves the published P(n) and composition
fits, then integrates d(epsilon/n)/dn=P/n^2 with an explicit independent
outer-drip anchor. It is not the original energy fit or an ETF minimization.
The Ioffe author page recommends thermodynamic reconstruction from selected
fitted functions: https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html .
"""
from __future__ import annotations
import math
from scipy.integrate import quad
from scipy.optimize import brentq
from . import bsk22_analytic as source

class BoundaryAnchoredBSk22:
    """Pressure-matched, explicitly anchored lower-inner analytic representation."""
    def __init__(self,*,outer_density:float,outer_pressure:float,outer_mu_b:float,
                 outer_epsilon:float|None=None,quadrature_epsabs:float=1e-10):
        self.outer_density=float(outer_density)
        self.outer_pressure=float(outer_pressure)
        self.outer_mu_b=float(outer_mu_b)
        self.quadrature_epsabs=float(quadrature_epsabs)
        if (not all(math.isfinite(v) and v>0 for v in
                    (self.outer_density,self.outer_pressure,self.outer_mu_b,self.quadrature_epsabs))
            or self.outer_density>=source.LOWER_INNER_MAX_DENSITY):
            raise ValueError("finite positive outer anchor below lower-inner maximum required")
        self.outer_epsilon=(self.outer_density*self.outer_mu_b-self.outer_pressure
                            if outer_epsilon is None else float(outer_epsilon))
        residual=self.outer_epsilon+self.outer_pressure-self.outer_density*self.outer_mu_b
        if not math.isfinite(self.outer_epsilon) or self.outer_epsilon<=0 or abs(residual)>1e-9:
            raise ValueError("outer anchor violates epsilon+P=n*mu_B")
        lo=max(source.NEUTRON_DRIP_DENSITY,self.outer_density)
        hi=source.LOWER_INNER_MAX_DENSITY
        fn=lambda n:source.state(n)["pressure_source_fit"]-self.outer_pressure
        if fn(lo)>=0 or fn(hi)<=0:
            raise ValueError("no positive-density-jump pressure match inside supported BSk22 interval")
        self.n0=brentq(fn,lo,hi,xtol=1e-16,rtol=1e-14)
        if self.n0<=self.outer_density*(1+1e-12):
            raise ValueError("inner density must exceed outer density for the local stable ordering")
        s=source.state(self.n0)
        if s["gamma_from_function_derivatives"]<=0:
            raise ValueError("nonpositive inner compressibility at matching pressure")
        self.e0=self.outer_mu_b-self.outer_pressure/self.n0

    def state(self,n_b:float)->dict:
        n=float(n_b)
        if not math.isfinite(n) or not self.n0<=n<=source.LOWER_INNER_MAX_DENSITY:
            raise ValueError("density outside boundary-anchored lower-inner support")
        s=source.state(n);p=s["pressure_source_fit"]
        integral,error=quad(lambda x:source.state(self.n0*math.exp(x))["pressure_source_fit"]/
                            (self.n0*math.exp(x)),0.,math.log(n/self.n0),
                            epsabs=self.quadrature_epsabs,epsrel=2e-12,limit=100)
        e=self.e0+integral;eps=n*e;mu=e+p/n
        gamma=s["gamma_from_function_derivatives"]
        cs2=gamma*p/(eps+p)
        # These source values remain explicit diagnostics, not a hidden correction.
        return dict(s,representation="boundary_anchored_published_pressure_analytic",
            thermodynamic_mode="outer_boundary_first_law_integral",
            source_energy_per_baryon=s["energy_per_baryon"],
            source_epsilon=s["epsilon"],source_mu_b_euler=s["mu_b_euler"],
            energy_per_baryon=e,energy_minus_neutron_rest_MeV=e-source.NEUTRON_REST_MEV,
            epsilon=eps,rho_g_cm3=eps*source.MEV_FM3_TO_G_CM3,
            pressure=p,pressure_from_energy_derivative=p,
            original_fit_thermodynamic_relative_residual=s["pressure_fit_thermodynamic_relative_residual"],
            pressure_fit_thermodynamic_relative_residual=0.,
            mu_b_euler=mu,mu_b_energy_derivative=mu,
            sound_speed_squared_over_c_squared=cs2,pressure_density_derivative=gamma*p/n,
            source_dlogP_dlogrho=s["dlogP_dlogrho"],dlogP_dlogrho=gamma*eps/(eps+p),
            source_sound_speed_squared_energy_derivative=s["sound_speed_squared_energy_derivative"],
            sound_speed_squared_energy_derivative=cs2,
            energy_shift_from_source_fit_MeV=e-s["energy_per_baryon"],
            source_composition_mu_n_residual_MeV=s["mu_n_composition_fit"]-mu,
            first_law_quadrature_error_MeV=error,
            anchor_inner_density=self.n0,
            limitation="Boundary-anchored published-pressure/composition approximation; no ETF/Strutinsky stationary solution claim.")

    def anchor_diagnostics(self)->dict:
        s=self.state(self.n0)
        return dict(outer_density=self.outer_density,inner_density=self.n0,
            relative_density_jump=self.n0/self.outer_density-1,
            pressure=self.outer_pressure,mu_b=self.outer_mu_b,
            outer_epsilon=self.outer_epsilon,inner_epsilon=s["epsilon"],
            pressure_residual=s["pressure"]-self.outer_pressure,
            chemical_potential_residual=s["mu_b_euler"]-self.outer_mu_b,
            density_jump_positive=self.n0>self.outer_density,
            inner_free_neutron_fraction=s["Y_n_free"],
            inner_electron_fraction=s["Y_e"],
            delta_dmu_dP=1/self.n0-1/self.outer_density,
            local_ordering="inner branch lower chemical potential immediately above the matched pressure",
            limitation="Local Gibbs ordering of two declared approximations; not a shared microscopic variational proof.")

    def eos_data(self,densities):
        import numpy as np
        grid=np.asarray(densities,dtype=float)
        # Reuse only the already verified composition carrier and source constants.
        data,_=source.eos_data(grid)
        rows=[self.state(n) for n in grid]
        col=lambda k:np.asarray([r[k] for r in rows])
        data.Q1=col("pressure")/grid
        data.Q3=col("mu_b_euler")/data.mn-1
        data.Q6=col("epsilon")/(grid*data.mn)-1
        data.Q7=data.Q6.copy()
        return data,rows
