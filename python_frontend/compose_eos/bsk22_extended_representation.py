"""Explicit bounded extension of the pressure-integrated BSk22 representation.

Pearson et al. 2018 / 2019 erratum, arXiv:1903.04981 Appendix C Eqs.80/85:
the analytic composition fits cover the inner crust. Table8 gives Z=40 through
0.034 fm^-3 and corrected Table12 proton drip at 0.059 fm^-3 for BSk22.
We expose only <=0.01 fm^-3 here. Free-proton zero follows the author's 2023
analytic-code convention. This is neither a new ETFSI minimization nor the
provider's original prescription to switch to a different functional at .002.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.integrate import quad
from . import bsk22_analytic as source
from .bsk22_thermodynamic_reconstruction import BoundaryAnchoredBSk22
from .eos_data import EOSData

MAX_SUPPORTED_DENSITY=.01

def extended_source_state(n_b):
    n=float(n_b)
    if not math.isfinite(n) or not source.NEUTRON_DRIP_DENSITY<=n<=MAX_SUPPORTED_DENSITY:
        raise ValueError("explicit BSk22 extension supports drip <= n <= .01")
    if n<=source.LOWER_INNER_MAX_DENSITY:return source.state(n)
    e,de,dde=source._energy(n)
    eps=n*(source.NEUTRON_REST_MEV+e)
    p,chi=source._pressure(eps*source.MEV_FM3_TO_G_CM3)
    p1,p2,p3,p4=source.PROTON_FRACTION_COEFFICIENTS
    ye=(p1*n**(-.75)-p2)*(1+p3*n)*(1+(p4*n)**4)
    x=(n-source.NEUTRON_DRIP_DENSITY)/(source.SOURCE_CORE_TRANSITION_DENSITY-source.NEUTRON_DRIP_DENSITY)
    a,b,c,d,k,f,g=source.FREE_NEUTRON_COEFFICIENTS
    yn=(a*x+d*x**k)/(1+b*x+c*x*x)+f*x*math.exp(g*(x-1))
    yh=1-yn;acell=40/ye;acluster=acell*yh
    if not 0<ye<.5 or not 0<=yn<1 or not acluster>40:
        raise ArithmeticError("extended analytic composition left physical support")
    mue=.51099893*math.sqrt(1+(386.159268*(3*math.pi**2*n*ye)**(1/3))**2)
    nf=n*yn
    mun=source.NEUTRON_REST_MEV+.211*source.NEUTRON_REST_MEV*nf**(2/3)/(1+9.41*math.sqrt(nf)+7.62*nf-259*nf*nf)
    return dict(n_b=n,epsilon=eps,pressure=p,pressure_source_fit=p,
        energy_per_baryon=source.NEUTRON_REST_MEV+e,mu_b_euler=(eps+p)/n,
        pressure_from_energy_derivative=n*de,mu_b_energy_derivative=source.NEUTRON_REST_MEV+e+de,
        gamma_from_function_derivatives=chi*(1+n*de/eps),dlogP_dlogrho=chi,
        pressure_fit_thermodynamic_relative_residual=(p-n*de)/p,
        mu_e=mue,mu_n_composition_fit=mun,Y_e=ye,Y_proton_total=ye,Y_neutron_total=1-ye,
        Y_mu=0.,Y_n_free=yn,Y_p_free=0.,Y_h=yh,A_cell=acell,Z_cell=40.,
        A_cluster=acluster,Z_cluster=40.,composition_baryon_residual=yn+yh-1,
        composition_charge_residual=ye-40/acluster*yh,
        source=source.SOURCE_URL,author_code_sha256=source.AUTHOR_CODE_SHA256,
        representation="explicit_published_fullcrust_analytic_extension")

class ExtendedAnchoredBSk22(BoundaryAnchoredBSk22):
    """First-law extension retaining the existing physical outer boundary."""
    def state(self,n_b):
        n=float(n_b)
        if n<=source.LOWER_INNER_MAX_DENSITY:
            return super().state(n)
        if not math.isfinite(n) or n>MAX_SUPPORTED_DENSITY:
            raise ValueError("bounded BSk22 extension requires n <= .01")
        s=extended_source_state(n)
        edge=super().state(source.LOWER_INNER_MAX_DENSITY)
        integral,error=quad(lambda x:extended_source_state(source.LOWER_INNER_MAX_DENSITY*math.exp(x))['pressure']/
             (source.LOWER_INNER_MAX_DENSITY*math.exp(x)),0.,math.log(n/source.LOWER_INNER_MAX_DENSITY),
             epsabs=self.quadrature_epsabs,epsrel=2e-12)
        e=edge['energy_per_baryon']+integral
        p=s['pressure'];eps=n*e;mu=e+p/n
        return dict(s,epsilon=eps,energy_per_baryon=e,mu_b_euler=mu,
            mu_b_energy_derivative=mu,pressure_from_energy_derivative=p,
            source_epsilon=s['epsilon'],source_energy_per_baryon=s['energy_per_baryon'],
            source_mu_b_euler=s['mu_b_euler'],
            original_fit_thermodynamic_relative_residual=s['pressure_fit_thermodynamic_relative_residual'],
            pressure_fit_thermodynamic_relative_residual=0.,
            energy_shift_from_source_fit_MeV=e-s['energy_per_baryon'],
            source_composition_mu_n_residual_MeV=s['mu_n_composition_fit']-mu,
            sound_speed_squared_over_c_squared=s['gamma_from_function_derivatives']*p/(eps+p),
            pressure_density_derivative=s['gamma_from_function_derivatives']*p/n,
            first_law_quadrature_error_MeV=error+edge['first_law_quadrature_error_MeV'],
            representation="explicit_extended_outer_anchored_pressure_integral",
            limitation="Analytic BSk22 approximation; no ETFSI solve; source composition chemical residual retained.")

    def eos_data(self,densities):
        grid=np.asarray(densities,dtype=float)
        if grid.ndim!=1 or len(grid)<3 or np.any(np.diff(grid)<=0):
            raise ValueError("at least three increasing densities required")
        rows=[self.state(n) for n in grid];col=lambda k:np.array([s[k] for s in rows])
        z=np.zeros(len(grid));mn=source.NEUTRON_REST_MEV
        data=EOSData(nb=grid,mn=mn,mp=mn-1.2933,Q1=col("pressure")/grid,Q2=z.copy(),
          Q3=col("mu_b_euler")/mn-1,Q4=-col("mu_e")/mn,Q5=z.copy(),
          Q6=col("epsilon")/grid/mn-1,Q7=col("epsilon")/grid/mn-1,q_extra=z.copy(),
          Iphase=np.full(len(grid),2),Yn=col("Y_n_free"),Yp=z.copy(),Ye=col("Y_e"),
          Ymu=z.copy(),Aav=col("A_cluster"),Zav=col("Z_cluster"),Yh=col("Y_h"),
          composition_file_present=True,composition_available=True,beta_equilibrated=True)
        return data,rows
