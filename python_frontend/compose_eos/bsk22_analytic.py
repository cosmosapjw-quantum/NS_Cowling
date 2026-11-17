"""Published BSk22 analytic lower-inner-crust representation.

Pearson et al., MNRAS481,2994(2018), corrected MNRAS486,768(2019),
AppendixC Eqs.C1-C9, TablesC1-C6 (HTML Eqs75-86, Tables20-25).
Author reference: http://www.ioffe.ru/astro/NSG/BSk/bskfit18.f, 2023-02-13.
This evaluates published coefficients; it does not solve the ETF/Strutinsky
variational problem or fit/copy a CompOSE target table. Pressure and energy
are separately fitted source functions; their thermodynamic mismatch is
reported explicitly. Support is restricted to the missing BSk22 lower-inner
segment, from source neutron drip through 0.002 fm^-3.
"""
from __future__ import annotations
import math

NEUTRON_DRIP_DENSITY = 2.69403e-4
LOWER_INNER_MAX_DENSITY = .002
SOURCE_CORE_TRANSITION_DENSITY = .0716068
NEUTRON_REST_MEV = 939.5654
ENERGY_GROUND_MEV = -9.1536
MEV_FM3_TO_G_CM3 = 1.782662e12  # exact convention of author fit
SOURCE_URL = "https://arxiv.org/html/1903.04981v1"
AUTHOR_CODE_URL = "http://www.ioffe.ru/astro/NSG/BSk/bskfit18.f"
AUTHOR_CODE_SHA256 = "4bd3b716f04e40c69165fa83b5cd8ecb03c30aa1f2b45342fb070cef09a3a99c"
ENERGY_COEFFICIENTS = (7.02e8,1.133e11,6.19e7,4.54e6,5.46e5,
                       15.24,.0683,8.86,4611.,48.07,2.697,81.7,7.05,1.50)
PRESSURE_COEFFICIENTS = (6.682,5.651,.00459,.14359,2.681,11.972,
    13.993,1.2904,2.665,-27.787,2.0140,4.09,14.135,28.03,-1.921,
    1.08,14.89,.098,11.67,4.75,-.037,14.10,11.9)
PROTON_FRACTION_COEFFICIENTS = (6.08e-4,3.42e-3,133.6,15.89)
FREE_NEUTRON_COEFFICIENTS = (173.2,204.8,57.3,61.4,1.86,.080,18.)

def _density(n_b:float)->float:
    n=float(n_b)
    if not math.isfinite(n) or not NEUTRON_DRIP_DENSITY<=n<=LOWER_INNER_MAX_DENSITY:
        raise ValueError("BSk22 lower-inner representation requires 2.69403e-4 <= n_b <= .002 fm^-3")
    return n

def _energy(n):
    """Energy minus neutron rest mass and analytic d/dlog(n), in MeV."""
    p=(None,)+ENERGY_COEFFICIENTS
    u={j:math.sqrt(p[j]*n) for j in (2,3,4,5)}
    low=(p[1]*n)**(7/6)*(1+u[4])/((1+u[2])*(1+u[3])*(1+u[5]))
    dl=7/6+.5*u[4]/(1+u[4])-.5*sum(u[j]/(1+u[j]) for j in (2,3,5))
    mid=p[6]*n**p[7]*(1+p[8]*n);dm=p[7]+p[8]*n/(1+p[8]*n)
    high=(p[10]*n)**p[11]/(1+p[12]*n);dh=p[11]-p[12]*n/(1+p[12]*n)
    w1=1/(1+p[9]*n);w2=1/(1+(p[13]*n)**p[14])
    terms=(low*w1,mid*(1-w1)*w2,high*(1-w2))
    e=ENERGY_GROUND_MEV+sum(terms)
    derivative=terms[0]*(dl-(1-w1))+terms[1]*(dm+w1-p[14]*(1-w2))+terms[2]*(dh+p[14]*w2)
    ddl=.25*u[4]/(1+u[4])**2-.25*sum(u[j]/(1+u[j])**2 for j in (2,3,5))
    ddm=p[8]*n/(1+p[8]*n)**2;ddh=-p[12]*n/(1+p[12]*n)**2
    l1=dl-(1-w1);l2=dm+w1-p[14]*(1-w2);l3=dh+p[14]*w2
    common1=w1*(1-w1);common2=p[14]**2*w2*(1-w2)
    second=terms[0]*(l1*l1+ddl-common1)+terms[1]*(l2*l2+ddm-common1-common2)+terms[2]*(l3*l3+ddh-common2)
    return e,derivative,second

def _logistic_negative(x):
    if x>=0:
        z=math.exp(-x);return z/(1+z)
    return 1/(1+math.exp(x))

def _pressure(rho):
    """Published P(rho) and exact derivative dlogP/dlogrho of that fit."""
    p=(None,)+PRESSURE_COEFFICIENTS;x=math.log10(rho)
    denom=1+p[4]*x
    a=(p[1]+p[2]*x+p[3]*x**3)/denom
    da=(p[2]+3*p[3]*x*x-p[4]*a)/denom
    f1=_logistic_negative(p[5]*(x-p[6]))
    f2=_logistic_negative(p[9]*(p[6]-x))
    f3=_logistic_negative(p[12]*(p[13]-x))
    f4=_logistic_negative(p[16]*(p[17]-x))
    b=p[7]+p[8]*x;c=p[10]+p[11]*x;d=p[14]+p[15]*x
    g=p[20]*(x-p[19]);h=p[23]*(x-p[22])
    tail1=p[18]/(1+g*g);tail2=p[21]/(1+h*h)
    logp=a*f1+b*f2+c*f3+d*f4+tail1+tail2-33.2047
    chi=(da*f1-a*p[5]*f1*(1-f1)+p[8]*f2+b*p[9]*f2*(1-f2)
        +p[11]*f3+c*p[12]*f3*(1-f3)+p[15]*f4+d*p[16]*f4*(1-f4)
        -2*tail1*g*p[20]/(1+g*g)-2*tail2*h*p[23]/(1+h*h))
    return 10**logp,chi

def state(n_b:float, *, thermodynamic_mode:str="source_fit")->dict:
    """Source-fit thermodynamics and full cell/cluster/free composition.

    Y_n_free and Y_p_free are free-particle baryon fractions; Y_h is the
    cluster baryon fraction, so Ye=Y_p_free+(Z_cluster/A_cluster)*Y_h.
    The source analytic fits smooth discontinuities and have finite errors.
    """
    if thermodynamic_mode not in ("source_fit","energy_derivative"):
        raise ValueError("thermodynamic_mode must be source_fit or energy_derivative")
    n=_density(n_b);e,de_dlogn,d2e_dlogn2=_energy(n)
    eps=n*(NEUTRON_REST_MEV+e);rho=eps*MEV_FM3_TO_G_CM3
    pressure_fit,chi=_pressure(rho);pressure_energy=n*de_dlogn
    pressure=pressure_fit if thermodynamic_mode=="source_fit" else pressure_energy
    sound_fit=chi*pressure_fit/eps
    sound_energy=n*(de_dlogn+d2e_dlogn2)/(eps+pressure_energy)
    sound=sound_fit if thermodynamic_mode=="source_fit" else sound_energy
    gamma=chi*(1+pressure_energy/eps) if thermodynamic_mode=="source_fit" else 1+d2e_dlogn2/de_dlogn
    p1,p2,p3,p4=PROTON_FRACTION_COEFFICIENTS
    ye=(p1*n**(-.75)-p2)*(1+p3*n)*(1+(p4*n)**4)
    x=(n-NEUTRON_DRIP_DENSITY)/(SOURCE_CORE_TRANSITION_DENSITY-NEUTRON_DRIP_DENSITY)
    a,b,c,d,k,f,g=FREE_NEUTRON_COEFFICIENTS
    yn_free=(a*x+d*x**k)/(1+b*x+c*x*x)+f*x*math.exp(g*(x-1))
    yp_free=0.  # Author 2023 correction: exactly zero below proton drip.
    yh=1-yn_free;zcell=40.;acell=zcell/ye;acluster=acell*yh
    if not (0<ye<1 and 0<=yn_free<1 and acluster>zcell):
        raise ArithmeticError("published composition fit left physical support")
    mu_e=.51099893*math.sqrt(1+(386.159268*(3*math.pi**2*n*ye)**(1/3))**2)
    nf=n*yn_free
    mu_n_composition=NEUTRON_REST_MEV+.211*NEUTRON_REST_MEV*nf**(2/3)/(1+9.41*math.sqrt(nf)+7.62*nf-259*nf*nf)
    mu_n_energy=NEUTRON_REST_MEV+e+de_dlogn
    return dict(model="BSk22",representation="published_analytic_ETFSI_fit",
        n_b=n,epsilon=eps,pressure=pressure,rho_g_cm3=rho,thermodynamic_mode=thermodynamic_mode,
        pressure_source_fit=pressure_fit,
        energy_per_baryon=NEUTRON_REST_MEV+e,energy_minus_neutron_rest_MeV=e,
        pressure_from_energy_derivative=pressure_energy,
        pressure_fit_thermodynamic_relative_residual=(pressure_fit-pressure_energy)/pressure_fit,
        dlogP_dlogrho=chi,gamma_from_function_derivatives=gamma,
        sound_speed_squared_over_c_squared=sound,
        sound_speed_squared_source_fit=sound_fit,sound_speed_squared_energy_derivative=sound_energy,
        mu_b_euler=(eps+pressure)/n,mu_b_energy_derivative=mu_n_energy,
        mu_n_composition_fit=mu_n_composition,mu_e=mu_e,
        Y_e=ye,Y_proton_total=ye,Y_neutron_total=1-ye,Y_mu=0.,
        Y_n_free=yn_free,Y_p_free=yp_free,Y_h=yh,
        A_cell=acell,Z_cell=zcell,A_cluster=acluster,Z_cluster=zcell,
        cluster_number_density=n/acell,
        composition_baryon_residual=yn_free+yp_free+yh-1,
        composition_charge_residual=ye-yp_free-zcell/acluster*yh,
        source=SOURCE_URL,author_code_sha256=AUTHOR_CODE_SHA256,
        limitation="Published finite-error representation; no independent ETF/Strutinsky minimization.")


def eos_data(densities, *, thermodynamic_mode="source_fit"):
    """Construct a composition-aware CompOSE carrier without joining endpoints."""
    import numpy as np
    from .eos_data import EOSData
    grid=np.asarray(densities,dtype=float)
    if grid.ndim!=1 or len(grid)<3 or np.any(np.diff(grid)<=0):
        raise ValueError("at least three strictly increasing densities required")
    rows=[state(n,thermodynamic_mode=thermodynamic_mode) for n in grid]
    col=lambda k:np.array([s[k] for s in rows])
    zero=np.zeros(len(grid));yes=np.ones(len(grid),dtype=bool)
    mn=NEUTRON_REST_MEV;mp=mn-1.2933
    data=EOSData(nb=grid,mn=mn,mp=mp,Q1=col("pressure")/grid,Q2=zero.copy(),
        Q3=col("mu_b_euler")/mn-1,Q4=-col("mu_e")/mn,Q5=zero.copy(),
        Q6=col("epsilon")/(grid*mn)-1,Q7=col("epsilon")/(grid*mn)-1,
        q_extra=zero.copy(),Iphase=np.full(len(grid),2,dtype=int),
        Yn=col("Y_n_free"),Yp=col("Y_p_free"),Ye=col("Y_e"),Ymu=zero.copy(),
        Aav=col("A_cluster"),Zav=col("Z_cluster"),Yh=col("Y_h"),
        composition_file_present=True,composition_available=True,beta_equilibrated=True,
        row_status=np.full(len(grid),2,dtype=int),
        Yn_present=yes.copy(),Yp_present=yes.copy(),Ye_present=yes.copy(),Ymu_present=yes.copy(),
        Abar_present=yes.copy(),Zbar_present=yes.copy(),Yh_present=yes.copy(),
        C3_applicable=yes.copy(),C4_applicable=yes.copy())
    return data,rows
