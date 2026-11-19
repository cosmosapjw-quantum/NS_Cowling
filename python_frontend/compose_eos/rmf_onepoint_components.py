"""Original-author vs local relativistic mean-field one-point component audit.

Read-only, no couplings fit to source. Fields use local convention:
phi positive and M*=M-g_sigma*phi, whereas Shen's sigma is negative.
Original GPPVA CompOSE packages do NOT publish these meson mean fields;
never label local fields as 'original provider' outputs.

Cold locally charge-neutral (or explicit fixed-charge) npe(mu) baryons,
n in fm^-3, P and epsilon in MeV/fm^3, mu and fields in MeV.
"""
from __future__ import annotations

import gzip
import hashlib
import math
from pathlib import Path

import numpy as np

from .cold_leptons import (HBARC,ELECTRON_MASS,MUON_MASS,free_fermion)
from . import gppva_rmf


SHEN_EOS4_T00_SHA256 = (
 "346b2701f7372817c516c3c945fdee0b3da5af9f08e2fe83f743951d1b5ebaef"
)
SHEN_ORIGINAL_M_NUCLEON_MEV=938.0
SHEN_ATOMIC_MASS_UNIT_MEV=931.494


def read_shen_eos4_t00_exact_row(path,log10_rho_gcm3=14.5,Yp=0.10):
    """Read a *single exact raw* 2020 author T=0 fixed-Yp TM1e row.

    The 16 scalar columns are mapped with published EOS4 User Note VI:
    log10(rho_B), n_B, Yp, F-938, Eint-931.494, S, A, Z,
    Mstar, Xn, Xp, Xalpha, XA, P, mun-938, mup-938.
    No interpolation, no source rebaselining, no lepton contributions.
    """
    path=Path(path)
    if not path.is_file():
        raise FileNotFoundError("original Shen EOS4 zero-temperature source missing")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SHEN_EOS4_T00_SHA256:
        raise ValueError("original Shen EOS4 source file hash drift")
    if not np.isfinite(log10_rho_gcm3) or not np.isfinite(Yp):
        raise ValueError("finite exact published grid values required")
    total=0
    selected=[]
    with gzip.open(path,"rt",encoding="utf-8",errors="replace") as stream:
        for line in stream:
            tokens=line.split()
            if len(tokens)!=16:
                continue
            try:row=[float(token.replace("D","E")) for token in tokens]
            except ValueError:
                continue
            if not np.all(np.isfinite(row)):
                raise ValueError("nonfinite author table row")
            total+=1
            if abs(row[0]-log10_rho_gcm3)<1e-6 and abs(row[2]-Yp)<1e-7:
                selected.append(row)
    if total!=7150:
        raise ValueError("original author table expected 65*110 numeric rows")
    if len(selected)!=1:
        raise ValueError("exact Shen EOS4 published grid point unavailable; no interpolation")
    row=selected[0]
    logRho,n,yp,F,Eint,entropy,A,Z,mstar,xn,xp,xalpha,xa,P,mun,mup=row
    if not 0<yp<1 or n<=0 or mstar<=0 or P<=0:
        raise ValueError("chosen point not physically tabulated as a uniform state")
    if abs(Eint-F-(SHEN_ORIGINAL_M_NUCLEON_MEV-SHEN_ATOMIC_MASS_UNIT_MEV))>.02:
        raise ValueError("author zero-temperature energy reference inconsistent")
    if abs(xa)>1e-8 or abs(xalpha)>1e-8 or abs(xp-yp)>1e-8 or abs(xn-(1-yp))>1e-8:
        raise ValueError("selected original-author point not uniform npe matter")
    return dict(
        status="EXACT_PRIMARY_AUTHOR_UNIFORM_TABULATION",source_sha256=SHEN_EOS4_T00_SHA256,
        author="Shen, Ji, Hu, Sumiyoshi (2020), EOS4/TM1e",raw_numeric_source_rows=total,
        reference_url="https://user.numazu-ct.ac.jp/~sumi/eos/table4/eos4.t00.gz",
        source_role="INDEPENDENT_ORIGINAL_TM1e_MODEL_AUTHOR_NOT_PAIS_GPPVA_GENERATOR",
        original_guide="https://user.numazu-ct.ac.jp/~sumi/eos/table4/guide_EOS4.pdf",
        temperature_mev=0.0,
        log10_rho_B_g_cm3=float(logRho),n_B_fm3=float(n),Yp=float(yp),
        free_energy_per_baryon_relative938_MeV=float(F),
        internal_energy_per_baryon_relative931p494_MeV=float(Eint),
        free_baryon_neutron_fraction=float(xn),free_baryon_proton_fraction=float(xp),
        alpha_fraction=float(xalpha),heavy_nucleus_fraction=float(xa),
        Mstar_MeV=float(mstar),pressure_baryons_MeV_fm3=float(P),
        mu_n_relative938_MeV=float(mun),mu_p_relative938_MeV=float(mup),
        source_baryon_rest_mass_MeV=SHEN_ORIGINAL_M_NUCLEON_MEV,
        no_source_interpolation=True,leptons_not_included_in_original_author_EOS4=True)


def nonlinear_component_budget(name,n_B,Yp,n_e,n_mu):
    """Analytic scalar/vector/fermion decomposition of local stationary state.

    Original source provides total P,e,Y only, never σ,ω,ρ or individual
    kinetic and meson pieces. This is a *local* decomposition at explicit
    baryon/proton/electron/muon densities, no external RMF assumptions.
    """
    vals=(float(n_B),float(Yp),float(n_e),float(n_mu))
    n,yp,ne,nm=vals
    if not all(map(math.isfinite,vals)) or n<=0 or not 0<=yp<=1 or min(ne,nm)<0:
        raise ValueError("valid positive density, bounded proton fraction, positive leptons")
    param=gppva_rmf.model(name)
    if param["kind"]!="nonlinear":
        raise ValueError("only original named nonlinear FSU2/TM1e with fields supported")
    b=gppva_rmf.baryons(name,n,yp)
    phi,w,r=float(b["phi"]),float(b["omega"]),float(b["rho"])
    density_n=n*(1-yp);density_p=n*yp;n3=density_p-density_n
    H=HBARC**3
    kn=free_fermion(density_n,b["effective_mass_n"])
    kp=free_fermion(density_p,b["effective_mass_p"])
    ke=free_fermion(ne,ELECTRON_MASS)
    km=free_fermion(nm,MUON_MASS)
    sigma2=.5*param["ms"]**2*phi**2/H
    sigma3=param["a3"]*phi**3/(3*H)
    sigma4=param["a4"]*phi**4/(4*H)
    omega2=.5*param["mw"]**2*w*w/H
    rho2=.5*param["mr"]**2*r*r/H
    omega4=param["b4"]*w**4/(4*H)
    cross=param["cross"]*w*w*r*r/H
    # On-shell field equations replace source terms in epsilon:
    # -W + omega*dW/domega + rho*dW/drho
    #  = .5mw²w² + .5mr²r² + 3b4w⁴/4 + 3cross w²r².
    E=dict(
        kinetic_neutrons=kn[0],kinetic_protons=kp[0],
        kinetic_electrons=ke[0],kinetic_muons=km[0],
        energy_sigma_mass=sigma2,energy_sigma_cubic=sigma3,
        energy_sigma_quartic=sigma4,
        energy_omega_mass=omega2,energy_rho_mass=rho2,
        energy_omega_quartic=3*omega4,
        energy_omega_rho_mixed=3*cross)
    P=dict(
        kinetic_neutrons=kn[1],kinetic_protons=kp[1],
        kinetic_electrons=ke[1],kinetic_muons=km[1],
        pressure_sigma_mass=-sigma2,pressure_sigma_cubic=-sigma3,
        pressure_sigma_quartic=-sigma4,
        pressure_omega_mass=omega2,pressure_rho_mass=rho2,
        pressure_omega_quartic=omega4,
        pressure_omega_rho_mixed=cross,
        pressure_rearrangement=n*float(b["rearrangement"]))
    esum=float(sum(E.values()))
    psum=float(sum(P.values()))
    epsilon_target=b["epsilon_baryons"]+ke[0]+km[0]
    pressure_target=b["pressure_baryons"]+ke[1]+km[1]
    mun=kn[2]+param["gw"]*w-param["gr"]*r/2+float(b["rearrangement"])
    mup=kp[2]+param["gw"]*w+param["gr"]*r/2+float(b["rearrangement"])
    epsilon_explicit=(
        kn[0]+kp[0]+ke[0]+km[0] +
        sigma2+sigma3+sigma4 -
        (omega2+rho2+omega4+cross) +
        (param["gw"]*w*n+param["gr"]*r*n3/2))
    p_explicit=(
        kn[1]+kp[1]+ke[1]+km[1] -
        (sigma2+sigma3+sigma4)+(omega2+rho2+omega4+cross) +
        n*float(b["rearrangement"]))
    return dict(
        model=param["name"],source_type="LOCAL_RMF_FIELDS_NOT_ORIGINAL_PROVIDER_FIELDS",
        n_B_fm3=n,Yp=yp,Ye=ne/n,Ymu=nm/n,
        mstar_MeV=float(b["effective_mass_n"]),
        local_sigma_author_sign_MeV=-phi,
        local_phi_positive_MeV=phi,local_omega0_MeV=w,
        local_rho30_MeV=r,
        rho_vertex_convention="g_rho*tau3/2",
        field_equation_residuals_fm3=dict(
            sigma=float(b["scalar_residual_fm3"]),
            omega=float(b["omega_residual_fm3"]),
            rho=float(b["rho_residual_fm3"])),
        epsilon_components_MeV_fm3={k:float(v) for k,v in E.items()},
        pressure_components_MeV_fm3={k:float(v) for k,v in P.items()},
        epsilon_baryons_MeV_fm3=float(b["epsilon_baryons"]),
        pressure_baryons_MeV_fm3=float(b["pressure_baryons"]),
        epsilon_total_MeV_fm3=float(esum),
        pressure_total_MeV_fm3=float(psum),
        epsilon_target_local_stationary_MeV_fm3=float(epsilon_target),
        pressure_target_local_stationary_MeV_fm3=float(pressure_target),
        epsilon_recomposition_residual_MeV_fm3=float(esum-epsilon_target),
        pressure_recomposition_residual_MeV_fm3=float(psum-pressure_target),
        epsilon_explicit_source_potential_residual_MeV_fm3=float(esum-epsilon_explicit),
        pressure_explicit_source_potential_residual_MeV_fm3=float(psum-p_explicit),
        mu_n_MeV=float(mun),mu_p_MeV=float(mup),
        mu_n_recomposition_residual_MeV=float(mun-b["mu_n"]),
        mu_p_recomposition_residual_MeV=float(mup-b["mu_p"]),
        beta_residual_MeV=float(mun-mup-ke[2]),
        charge_residual_fm3=float(density_p-ne-nm),
        chemical_potential_neutron_euler_check_local_MeV=float((esum+psum)/n-mun),
        c_hbar_kept_in_convention=True,
        field_terms_are_author_theory_but_numerical_fields_are_local=True)
