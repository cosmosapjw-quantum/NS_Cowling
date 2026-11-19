#!/usr/bin/env python3
"""Plot and CSV: immutable CompOSE GPPVA source against local RMF at exact knots.

All deltas are LOCAL - ORIGINAL and evaluated at the source's own unmodified
phase-6 baryon density. Original Q3 and the Euler relation are distinct.
The provider does NOT publish component-resolved meson fields or pressures.
"""
from __future__ import annotations
import csv,hashlib,json,os
from pathlib import Path

os.environ.setdefault("MPLBACKEND","Agg")
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/cowling_original_rmf_onepoint_20261119"
INPUT=DATA/"RESULT.json"
OUTPUT=DATA/"figures"
COLORS=dict(source="#25739A",local="#D37C37",Q3="#673AB7",Euler="#19A38F",nonlinear="#A15A5A")

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
    data=json.loads(INPUT.read_text())
    models=("FSU2","TM1e","DDME2")
    vals=[data["models"][m] for m in models]
    for m,v in zip(models,vals):
        if not v["selected_density_is_original_provider_exact_knot"] or not v["no_source_interpolation"]:
            raise ValueError("plotting cannot use interpolated or altered source "+m)
    n=np.asarray([v["source"]["original_source_n_B_fm3"] for v in vals])
    ps=np.asarray([v["source"]["original_pressure_MeV_fm3"] for v in vals])
    pl=np.asarray([v["local_beta_equilibrium"]["pressure"] for v in vals])
    es=np.asarray([v["source"]["original_epsilon_MeV_fm3"] for v in vals])
    el=np.asarray([v["local_beta_equilibrium"]["epsilon"] for v in vals])
    q3=np.asarray([v["source"]["original_mu_B_Q3_MeV"] for v in vals])
    euler=np.asarray([v["source"]["original_mu_B_Euler_MeV"] for v in vals])
    mul=np.asarray([v["local_beta_equilibrium"]["muB"] for v in vals])
    dp=(pl/ps-1)*100.
    de=el-es
    dmq=mul-q3
    dme=mul-euler
    x=np.arange(len(models))
    labels=[f"{m}\n$n_B$={v:.6f}" for m,v in zip(models,n)]
    OUTPUT.mkdir(parents=True,exist_ok=True)

    fig,ax=plt.subplots(2,2,figsize=(13.3,8.4),constrained_layout=True)
    fig.suptitle("GPPVA CompOSE vs independently evaluated RMF | exact phase-6 knots",
                 fontsize=16,fontweight="bold")
    axs=ax[0,0]
    axs.bar(x-.18,ps,.35,label="CompOSE original",color=COLORS["source"])
    axs.bar(x+.18,pl,.35,label="Local RMF",color=COLORS["local"])
    axs.set_xticks(x,labels)
    axs.set_ylabel(r"Total pressure $P$ (MeV fm$^{-3}$)")
    axs.set_title("A | Original and local pressure")
    axs.legend(loc="upper right",fontsize=9)
    for i,(a,b) in enumerate(zip(ps,pl)):
        axs.annotate(f"{a:.6f}",(x[i]-.18,a),xytext=(0,4),
                     textcoords="offset points",ha="center",fontsize=8)
        axs.annotate(f"{b:.6f}",(x[i]+.18,b),xytext=(0,4),
                     textcoords="offset points",ha="center",fontsize=8)
    axs.set_ylim(0,max(pl)*1.20)
    axs.grid(axis="y",alpha=.18)

    axs=ax[0,1]
    axs.bar(x,dp,color=[COLORS["local"],COLORS["local"],COLORS["source"]],width=.53)
    axs.axhline(0,color="0.25",lw=.7)
    axs.set_xticks(x,labels)
    axs.set_title(r"B | Relative pressure discrepancy")
    axs.set_ylabel(r"$(P_{local}/P_{CompOSE}-1)\times100$ (%)")
    axs.set_ylim(0,max(dp)*1.25)
    axs.grid(axis="y",alpha=.18)
    for i,v in enumerate(dp):
        axs.annotate(f"{v:+.4f}%",(i,v),xytext=(0,6),
                     textcoords="offset points",ha="center",fontsize=10)

    axs=ax[1,0]
    axs.bar(x,de,color=[COLORS["local"],COLORS["local"],COLORS["source"]],width=.53)
    axs.set_xticks(x,labels)
    axs.set_title(r"C | Total energy-density discrepancy")
    axs.set_ylabel(r"$\varepsilon_{local}-\varepsilon_{CompOSE}$ (MeV fm$^{-3}$)")
    axs.set_ylim(0,max(de)*1.25)
    axs.grid(axis="y",alpha=.18)
    for i,v in enumerate(de):
        axs.annotate(f"{v:+.6f}",(i,v),xytext=(0,6),
                     textcoords="offset points",ha="center",fontsize=10)

    axs=ax[1,1]
    axs.bar(x-.18,dmq,.35,color=COLORS["Q3"],
            label=r"$\mu_B^{local}-\mu_B^{Q_3,CompOSE}$")
    axs.bar(x+.18,dme,.35,color=COLORS["Euler"],
            label=r"$\mu_B^{local}-(\varepsilon+P)_{CompOSE}/n_B$")
    axs.set_xticks(x,labels)
    axs.set_title(r"D | Chemical-potential discrepancy")
    axs.set_ylabel(r"$\Delta\mu_B$ (MeV)")
    axs.set_ylim(0,max(dmq)*1.3)
    axs.legend(loc="upper right",fontsize=7)
    axs.grid(axis="y",alpha=.18)

    fig.text(.50,-.015,
        "Different models use their closest published n_B to 0.2 fm^-3; no regridding. "
        "Local nonlinear meson fields are not published in the original CompOSE EOS.",
        ha="center",va="top",fontsize=9,color="0.3")
    png=OUTPUT/"COMPOSE_vs_LOCAL_RMF_exact_onepoint.png"
    svg=OUTPUT/"COMPOSE_vs_LOCAL_RMF_exact_onepoint.svg"
    fig.savefig(png,dpi=200,bbox_inches="tight",facecolor="white")
    fig.savefig(svg,bbox_inches="tight",facecolor="white")
    plt.close(fig)

    csvname=OUTPUT/"COMPOSE_vs_LOCAL_RMF_exact_onepoint.csv"
    cols=[
        "model","n_B_original_fm3","P_CompOSE_MeV_fm3","P_local_MeV_fm3",
        "DeltaP_local_minus_CompOSE_MeV_fm3","DeltaP_relative_percent",
        "epsilon_CompOSE_MeV_fm3","epsilon_local_MeV_fm3",
        "Deltaepsilon_MeV_fm3","Deltaepsilon_relative_percent",
        "muB_CompOSE_Q3_MeV","muB_CompOSE_Euler_MeV","muB_local_MeV",
        "Deltamu_local_minus_Q3_MeV","Deltamu_local_minus_Euler_MeV",
        "CompOSE_Q3_minus_Euler_MeV","original_Yp","local_Yp",
        "source_has_original_meson_field_data"]
    with csvname.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
        for i,m in enumerate(models):
            z=vals[i];src=z["source"];loc=z["local_beta_equilibrium"]
            w.writerow(dict(model=m,n_B_original_fm3=n[i],
               P_CompOSE_MeV_fm3=ps[i],P_local_MeV_fm3=pl[i],
               DeltaP_local_minus_CompOSE_MeV_fm3=pl[i]-ps[i],
               DeltaP_relative_percent=dp[i],
               epsilon_CompOSE_MeV_fm3=es[i],epsilon_local_MeV_fm3=el[i],
               Deltaepsilon_MeV_fm3=de[i],
               Deltaepsilon_relative_percent=100*de[i]/es[i],
               muB_CompOSE_Q3_MeV=q3[i],muB_CompOSE_Euler_MeV=euler[i],
               muB_local_MeV=mul[i],
               Deltamu_local_minus_Q3_MeV=dmq[i],
               Deltamu_local_minus_Euler_MeV=dme[i],
               CompOSE_Q3_minus_Euler_MeV=q3[i]-euler[i],
               original_Yp=src["original_Yp"],local_Yp=loc["Yp"],
               source_has_original_meson_field_data=False))
    outcome=dict(schema="ns_cowling.compose_local_plot_exact_onepoint.v1",
       original_source_sha256=sha(INPUT),
       source_model_names=list(models),
       source_exact_density_fm3=n.tolist(),
       input_source_has_original_meson_fields=False,
       no_input_interpolation=True,
       scientific_interpretation="distinct source and local EOS, no source-level meson data; model-specific near 0.2 fm^-3",
       figures={
         "png":dict(relative_path=str(png.relative_to(ROOT)),sha256=sha(png),bytes=png.stat().st_size),
         "svg":dict(relative_path=str(svg.relative_to(ROOT)),sha256=sha(svg),bytes=svg.stat().st_size),
         "csv":dict(relative_path=str(csvname.relative_to(ROOT)),sha256=sha(csvname),bytes=csvname.stat().st_size)},
       max_abs_P_difference_MeV_fm3=float(max(abs(pl-ps))),
       max_relative_pressure_percent=float(max(abs(dp))),
       figure_created_from_actual_report=True)
    (OUTPUT/"FIGURE_RECEIPT.json").write_text(json.dumps(outcome,indent=2)+"\n")
    print("COMPose_PLOT_READY",png,"bytes",png.stat().st_size,
      "SVG",svg.stat().st_size,"CSV",csvname.stat().st_size,flush=True)
if __name__=="__main__":main()
