"""Publication plots: eight nonlinear EOS interpolation methods only.

Explicit user constraint: linear must never occur on displayed axes/legends.
The native physics failure of linear P(mu) is kept in machine diagnostic JSON.
"""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/native_dense_interpolation_20261118"
SUM=json.loads((DATA/"aggregated"/"SUMMARY.json").read_text())
OUT=ROOT/"docs/native_dense_interpolation_20261118/figures"
OUT.mkdir(parents=True,exist_ok=True)
METHODS=["Hermite","Exponential","Rational","Bernstein","Power","PCHIP","Akima","Cubic spline"]
KEYS=[k.lower().replace(" spline","") for k in METHODS]
assert all(k in SUM["interpolation_methods_plotted"] for k in KEYS)
assert "linear" not in KEYS
pretty={"rmf":"Native uniform RMF (5 models)","core":"Native core functionals (5 models)","tf":"Native spatial-TF N128 droplets (5 models)"}
figfiles=[]
for family in ("rmf","core","tf"):
    met=SUM["family_summary"][family]
    values=[met[k]["rms_baryon_relative"]*1e6 if met[k]["rms_baryon_relative"] is not None else np.nan for k in KEYS]
    coverage=[met[k]["coverage_fraction"]*100 for k in KEYS]
    physical=[met[k]["physical_fraction"]*100 for k in KEYS]
    fig,ax=plt.subplots(figsize=(11.4,5.8))
    bars=ax.bar(METHODS,values)
    for bar,eligible in zip(bars,[met[k]["qualifies_for_global_accuracy_rank"] for k in KEYS]):
        if not eligible:
            bar.set_hatch("//")
            bar.set_alpha(.45)
    for i,(num,cover,phy) in enumerate(zip(values,coverage,physical)):
        if np.isfinite(num):
            ax.annotate(f"{num:.3g} ppm\nC={cover:.1f}%, V={phy:.1f}%",
                xy=(i,num),xytext=(0,4),textcoords="offset points",ha="center",va="bottom",fontsize=7.9)
    ax.set_yscale("log")
    ax.set_ylim(bottom=max(min(v for v in values if np.isfinite(v))*.5,1e-3),
                top=max(v for v in values if np.isfinite(v))*3)
    ax.grid(axis="y",alpha=.22)
    ax.set_title(f"{pretty[family]}: withheld native baryon-density RMS error")
    ax.set_ylabel("RMS relative error (ppm), logarithmic scale")
    ax.tick_params(axis="x",rotation=23,labelsize=9)
    ax.text(.5,-.23,"C: coverage of withheld points, V: sampled causal/stable fraction of native-causal points.\n"
        "Hatched = insufficient full-family coverage/physical admission; do not rank by subset-only error.",
        ha="center",va="top",transform=ax.transAxes,fontsize=8)
    fig.tight_layout()
    fname=OUT/f"{family.upper()}_NATIVE_HOLDOUT_NO_LINEAR.png"
    fig.savefig(fname,dpi=190,bbox_inches="tight")
    plt.close(fig);figfiles.append(fname)
fig,ax=plt.subplots(figsize=(11.6,5.5))
xx=np.arange(len(KEYS))
for i,family in enumerate(("rmf","core","tf")):
    ratio=[SUM["family_summary"][family][k]["coverage_fraction"]*100 for k in KEYS]
    ax.bar(xx+(.27*(i-1)),ratio,width=.265,label=family.upper())
ax.set_ylim(0,106)
ax.set_ylabel("Native withheld points covered (%)")
ax.set_xticks(xx,METHODS,rotation=22,ha="right")
ax.set_title("Interpolation success coverage independent of reported RMS error")
ax.legend(loc="lower right")
ax.grid(axis="y",alpha=.2)
fig.tight_layout()
fname=OUT/"NATIVE_HOLDOUT_COVERAGE_NO_LINEAR.png";fig.savefig(fname,dpi=190)
plt.close(fig);figfiles.append(fname)

models=SUM["by_model"]
order=[(f,m) for f,mods in [("rmf",("NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL",
          "DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0")),
         ("core",("APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4")),
         ("tf",("DDME2","DDMEX","NL3","NL3wrL55","TM1e"))] for m in mods]
look={(r["family"],r["model"],r["method"]):r for r in models}
matrix=np.zeros((len(order),len(KEYS)))
mask=[]
for i,(f,m) in enumerate(order):
    cells=[]
    for j,method in enumerate(KEYS):
        row=look[(f,m,method)]
        v=row["rms_baryon_relative"]
        matrix[i,j]=np.log10(v*1e6) if v is not None and v>0 else np.nan
        cells.append(row["coverage_fraction"]>=.995 and row["physical_fraction"]>=.995)
    mask.append(cells)
fig,ax=plt.subplots(figsize=(12.5,7.4))
im=ax.imshow(matrix,aspect="auto",cmap="viridis",interpolation="nearest",
             vmin=np.nanmin(matrix),vmax=np.nanmax(matrix))
ax.set_yticks(range(len(order)),[f"{family.upper()}: {model}" for family,model in order],fontsize=8)
ax.set_xticks(range(len(KEYS)),METHODS,rotation=29,ha="right")
ax.set_title("Direct native holdout error by source and interpolation method")
for i in range(len(order)):
    for j in range(len(KEYS)):
        lab="*" if not mask[i][j] else ""
        val=10**matrix[i,j]
        if np.isfinite(val):
            ax.text(j,i,f"{val:.2g}{lab}",ha="center",va="center",fontsize=7,
                    color="white" if matrix[i,j]<(np.nanmin(matrix)+np.nanmax(matrix))/2 else "black")
fig.colorbar(im,ax=ax,label="log10 RMS baryon-density error (ppm)")
fig.text(.5,.01,"* Insufficient coverage or sampled physical admission: error is conditional on predicted subset",ha="center",fontsize=8)
fig.tight_layout(rect=(0,.03,1,1))
fname=OUT/"NATIVE_MODEL_METHOD_HEATMAP_NO_LINEAR.png"
fig.savefig(fname,dpi=190,bbox_inches="tight");plt.close(fig);figfiles.append(fname)
receipt=dict(status="CREATED_VERIFIED",exclude_linear_from_all_figures=True,
              methods=KEYS,figures=[dict(relative_path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
                        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in figfiles])
(DATA/"aggregated"/"FIGURES_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
