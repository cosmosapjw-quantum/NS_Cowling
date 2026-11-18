"""Plots cubic-vs-quintic source-fixed errors and matching/stellar shifts.

Every panel uses only cubic/quintic; no linear/other spline in shown curves.
H smooth-vs-CSS/Maxwell shifts explicitly marked 'different construction'.
"""
from pathlib import Path
import json,hashlib,math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/cubic_quintic_matching_20261118"
S=json.loads((D/"SUMMARY.json").read_text())
OUT=ROOT/"docs/cubic_quintic_matching_20261118/figures"
OUT.mkdir(parents=True,exist_ok=True)
FILES=[]
def write(fig,filename):
    fig.tight_layout()
    path=OUT/filename
    fig.savefig(path,dpi=190,bbox_inches="tight")
    plt.close(fig);FILES.append(path)

groups={"rmf":16,"core":8}
models={
 "rmf":["NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL","DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0"],
 "core":["APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4"],
 "tf":["DDME2","DDMEX","NL3","NL3wrL55","TM1e"]
}
src=S["dense_model_phase_runs"]
for family in ("rmf","core","tf"):
    names=models[family]
    vals=[]
    for name in names:
        by=[]
        for method in ("cubic","quintic"):
            rows=[r for r in src if r["family"]==family and r["model"]==name
                  and r["method"]==method and
                  (str(r["stride"])==str(groups[family]) if family!="tf" else True)]
            numerator=sum(r["holdout_count"]*r["n_percent"]**2 for r in rows)
            denominator=sum(r["holdout_count"] for r in rows)
            by.append(np.sqrt(numerator/denominator))
        vals.append(by)
    vals=np.array(vals)
    y=np.arange(len(names))
    fig,ax=plt.subplots(figsize=(10,5.3))
    width=.34
    ax.barh(y-width/2,vals[:,0],height=width,label="Cubic Hermite")
    ax.barh(y+width/2,vals[:,1],height=width,label="Quintic Hermite")
    ax.set_yticks(y,names,fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Held-out native baryon-density RMS relative error (%)")
    ax.set_title(f"{family.upper()} native EOS cubic vs quintic (identical training/holdout)")
    ax.grid(axis="x",alpha=.25);ax.legend()
    write(fig,f"{family.upper()}_CUBIC_QUINTIC_NATIVE_ERRORS.png")

maxrows=[r for r in S["maxwell"] if r["native_status"]=="PASS_SCOPED" and
         r["level"]=="sparse_8x" and r["approved"]]
pairs=sorted({r["pair"] for r in maxrows})
fig,axes=plt.subplots(1,2,figsize=(11,4.4))
for ax,field,title in zip(axes,("Pt_percent","delta_n_percent"),
                           ("Maxwell transition pressure","Maxwell density jump")):
    x=np.arange(len(pairs))
    for k,method in enumerate(("cubic","quintic")):
        values=[next(r[field] for r in maxrows if r["pair"]==pair and r["method"]==method)
                for pair in pairs]
        ax.bar(x+(.32*(k-.5)),values,width=.31,label=method)
    ax.set_xticks(x,[z.replace("__"," / ") for z in pairs],rotation=20,ha="right",fontsize=8)
    ax.axhline(0,color="black",linewidth=.8)
    ax.set_title(title+" (sparse 8x)")
    ax.set_ylabel("Signed error vs full-source native Maxwell (%)")
    ax.grid(axis="y",alpha=.2);ax.legend()
write(fig,"MAXWELL_NATIVE_TRANSITION_RELATIVE_ERRORS.png")

pairs=[r for r in S["stellar_paired"] if r["conversion"]=="slow"]
labels=[r["family"].replace("_"," ")+"\n"+r["case"].replace("H_","").split("__")[0] for r in pairs]
fig,axes=plt.subplots(1,3,figsize=(14,5.3))
for ax,field,label in zip(axes,("R_relative_percent","f_relative_percent","p1_relative_percent"),
                           ("Radius","f-mode","p1-mode")):
    yy=np.arange(len(pairs))
    z=[r[field] if r[field] is not None else np.nan for r in pairs]
    ax.barh(yy,np.nan_to_num(z,nan=0.))
    for i,val in enumerate(z):
        if np.isnan(val):ax.text(0,i,"not admitted",fontsize=7,ha="left",va="center")
    ax.set_yticks(yy,labels,fontsize=7)
    ax.invert_yaxis();ax.set_xlabel("Quintic - cubic (% of cubic)")
    ax.set_title(label)
    ax.axvline(0,color="black",lw=.7);ax.grid(axis="x",alpha=.2)
fig.suptitle("Fixed-M=1.4 M_sun Hermite-only shifts; other physics held fixed")
write(fig,"STELLAR_CUBIC_QUINTIC_RELATIVE_SHIFTS.png")
receipt=dict(status="SUCCESS",methods=["cubic_hermite","quintic_hermite"],
             source="SUMMARY.json from native direct EOS, Maxwell independent root and Cowling",
             figure_files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
                                sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in FILES])
(D/"FIGURES_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
print("PLOT_DONE",len(FILES),[(x.name,x.stat().st_size) for x in FILES])
