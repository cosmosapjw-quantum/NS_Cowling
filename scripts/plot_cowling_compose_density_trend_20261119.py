#!/usr/bin/env python3
"""Plot all exact original GPPVA phase-6 knots versus independently evaluated RMF.

FSU2/TM1e/DDME2 are independently sampled at *their own* original density
locations n_B in [0.08,0.45] fm^-3. No data interpolation, no EOS fitting.
Original Q3 and the algebraic Euler mu_B are separate curves.
"""
from __future__ import annotations
import csv
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"data/cowling_rmf_reconciliation_20261119/EXACT_DENSITY_AUDIT.json"
OUT = ROOT/"data/cowling_original_rmf_onepoint_20261119/figures"
MODELS = ("FSU2", "TM1e", "DDME2")
COLORS = {"FSU2":"#D05B41", "TM1e":"#2676AD", "DDME2":"#58885A"}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    original = json.loads(SOURCE.read_text())["models"]
    OUT.mkdir(parents=True, exist_ok=True)
    combined=[]
    fig,axes=plt.subplots(2,2,figsize=(13.8,9.1),constrained_layout=True)
    fig.suptitle("Original GPPVA CompOSE vs independently evaluated RMF: all exact phase-6 density knots",
                 weight="bold", fontsize=14)
    for name in MODELS:
        rows=original[name]["original_beta_equilibrium_fraction_comparison"]["homogeneous_rows"]
        if len(rows)<40:
            raise ValueError("missing original locked exact-density rows "+name)
        densities=np.asarray([q["density_fm3"] for q in rows],float)
        if not np.all(np.diff(densities)>0):
            raise ValueError("the source phase6 density grid is not strictly increasing")
        P=np.asarray([q["pressure_provider"] for q in rows],float)
        p=np.asarray([q["pressure_model"] for q in rows],float)
        E=np.asarray([q["epsilon_provider"] for q in rows],float)
        e=np.asarray([q["epsilon_model"] for q in rows],float)
        mu_Q3=np.asarray([q["mu_Q3_provider"] for q in rows],float)
        mu_Euler=np.asarray([q["mu_Euler_provider"] for q in rows],float)
        mu_local=np.asarray([q["mu_B_model"] for q in rows],float)
        if not all(np.all(np.isfinite(a)) for a in (densities,P,p,E,e,mu_Q3,mu_Euler,mu_local)):
            raise ValueError("nonfinite original/provider table")
        pressure_pct=100*(p/P-1)
        energy_per_baryon=(e-E)/densities
        q3gap=mu_local-mu_Q3
        eulergap=mu_local-mu_Euler
        original_euler_gap=mu_Q3-mu_Euler
        axes[0,0].plot(densities, pressure_pct, label=name, color=COLORS[name],lw=1.9)
        axes[0,1].plot(densities, energy_per_baryon,label=name,color=COLORS[name],lw=1.9)
        axes[1,0].plot(densities,q3gap,label=name,color=COLORS[name],lw=1.9)
        axes[1,1].plot(densities,original_euler_gap,label=name,color=COLORS[name],lw=1.9)
        for i,n in enumerate(densities):
            combined.append(dict(model=name,n_B_fm3=float(n),
               P_CompOSE_MeV_fm3=float(P[i]),P_local_MeV_fm3=float(p[i]),
               DeltaP_local_minus_CompOSE_MeV_fm3=float(p[i]-P[i]),
               Relative_deltaP_percent=float(pressure_pct[i]),
               epsilon_CompOSE_MeV_fm3=float(E[i]),epsilon_local_MeV_fm3=float(e[i]),
               Delta_epsilon_per_baryon_MeV=float(energy_per_baryon[i]),
               mu_Q3_CompOSE_MeV=float(mu_Q3[i]),
               mu_Euler_CompOSE_MeV=float(mu_Euler[i]),mu_local_MeV=float(mu_local[i]),
               Delta_mu_local_minus_Q3_MeV=float(q3gap[i]),
               Delta_mu_local_minus_Euler_MeV=float(eulergap[i]),
               Delta_mu_Q3_minus_Euler_CompOSE_MeV=float(original_euler_gap[i]),
               tier="original_source_knots_vs_local_independently_evaluated_RMF",
               original_meson_components_available=False))
    labels=(
     "Pressure: 100 x (local / original - 1), %",
     "Energy per baryon: (epsilon_local - epsilon_CompOSE) / n_B, MeV",
     "Baryon chemical potential: mu_local - mu_Q3, MeV",
     "CompOSE internal consistency: mu_Q3 - (epsilon+P)/n_B, MeV")
    for ax,label in zip(axes.flatten(),labels):
        ax.set_xlabel(r"Original $n_B$ (fm$^{-3}$)")
        ax.set_ylabel(label,fontsize=9)
        ax.grid(alpha=.22)
        ax.axhline(0,lw=.7,color="#424242")
        ax.legend(loc="best",fontsize=8)
    axes[0,0].set_title("A | Relative total-pressure mismatch")
    axes[0,1].set_title("B | Total energy-functional discrepancy")
    axes[1,0].set_title("C | Local minus original CompOSE Q3 chemical potential")
    axes[1,1].set_title("D | Q3 vs Euler relation within original CompOSE")
    for ax in axes.flatten():
        ax.axvline(.2,color="#8D8D8D",alpha=.5,lw=.9,ls="--")
    outfile=OUT/"COMPOSE_vs_LOCAL_RMF_original_density_trend.png"
    svg=OUT/"COMPOSE_vs_LOCAL_RMF_original_density_trend.svg"
    fig.savefig(outfile,dpi=190,bbox_inches="tight",facecolor="white")
    fig.savefig(svg,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    csvpath=OUT/"COMPOSE_vs_LOCAL_RMF_original_density_trend.csv"
    with csvpath.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(combined[0]))
        writer.writeheader()
        writer.writerows(combined)
    report=dict(schema="ns_cowling.original_rmf_all_exact_density_visualization.v1",
        source_relpath=str(SOURCE.relative_to(ROOT)),
        source_SHA256=sha(SOURCE),
        exact_original_knots=True,
        original_source_density_range_fm3=[.08,.45],
        no_source_interpolation_or_physics_fitting=True,
        total_rows=len(combined),
        model_count={name:sum(x["model"]==name for x in combined) for name in MODELS},
        observations={
          name:dict(n=len(original[name]["original_beta_equilibrium_fraction_comparison"]["homogeneous_rows"]),
             min_n=float(min(q["density_fm3"] for q in original[name]["original_beta_equilibrium_fraction_comparison"]["homogeneous_rows"])),
             max_n=float(max(q["density_fm3"] for q in original[name]["original_beta_equilibrium_fraction_comparison"]["homogeneous_rows"]))) for name in MODELS},
        figures={p.name:dict(sha256=sha(p),size_bytes=p.stat().st_size)
                 for p in (outfile,svg,csvpath)},
        caveat="original GPPVA has no released mean fields; data show total EOS discrepancies only",
        status="COMPUTED_FROM_FROZEN_STAGE05_ORIGINAL_EXACT_DENSITY_AUDIT")
    (OUT/"DENSITY_TREND_RECEIPT.json").write_text(json.dumps(report,indent=2)+"\n")
    print("COMPOSE_DENSITY_TREND_READY",len(combined),report["model_count"],
          outfile.stat().st_size,csvpath.stat().st_size,flush=True)
if __name__=="__main__":main()
