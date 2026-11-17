#!/usr/bin/env python3
"""Plot recorded finite-domain Maxwell successes without recomputing EOS physics."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np

CATEGORY = {
    "old_native_crust__new_native_core": (0, "#176B9C", "Native crust + native core"),
    "reference_assisted_hybrid": (1, "#80559B", "Public crust + native core"),
    "new_public_crust__old_public_GM1": (2, "#B56822", "Public crust + public core"),
}
NUMERIC_FIELDS = (
    "n_low", "n_high", "pressure_transition", "mu_transition",
    "density_jump", "latent_energy_density",
)
CSV_FIELDS = (
    "pair_id", "category", "crust", "core", "crust_source_id", "core_source_id",
    "candidate_index", "crust_component", "core_component",
    *NUMERIC_FIELDS, "source_json_path",
)

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def extract(campaign):
    rows=[]
    for ci, cell in enumerate(campaign["cells"]):
        for ri, candidate in enumerate(cell["candidate_results"]):
            if candidate.get("status") != "GENERATED_VERIFIED_FINITE_DOMAIN":
                continue
            native=candidate["v2_report"]["native_report"]
            row={"pair_id":cell["id"], **{k:cell[k] for k in (
                    "category","crust","core","crust_source_id","core_source_id")},
                 "candidate_index":ri,
                 "crust_component":candidate["crust_component"]["id"],
                 "core_component":candidate["core_component"]["id"],
                 "source_json_path":f"cells[{ci}].candidate_results[{ri}].v2_report.native_report",
                 **{k:native[k] for k in NUMERIC_FIELDS}}
            if not all(np.isfinite(row[k]) for k in NUMERIC_FIELDS):
                raise ValueError("Recorded transition contains nonfinite values")
            if not (0<row["n_low"]<row["n_high"] and row["pressure_transition"]>0):
                raise ValueError("Recorded transition is outside plotting domain")
            if row["category"] not in CATEGORY:
                raise ValueError("Unrecognized success category; update explicit legend")
            rows.append(row)
    return sorted(rows,key=lambda r:(CATEGORY[r["category"]][0],
                  1 if "GM_GM1"==r["core"] else 0,
                  0 if r["crust"].startswith("SPG_") else 1,
                  r["crust"],r["candidate_index"]))

def label(row):
    crust=row["crust"].replace("SPG_","SPG ").replace("GPPVA_","GPPVA ").replace("XMLSLZ_","XMLSLZ ")
    if row["core"]=="PK1_CANONICAL":
        core="PK1"
    elif row["core"]=="GM_GM1":
        core="GM1 [public]" if row["core_source_id"].startswith("public_original:") else "GM1 [native]"
    else:
        core=row["core"]
    return f"{crust}  "+r"$\rightarrow$"+f"  {core}"

def make_figure(rows,campaign):
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":10.5,
        "axes.spines.top":False,"axes.spines.right":False,
        "axes.edgecolor":"#9AA3AE","axes.labelcolor":"#243244",
        "xtick.color":"#475569","ytick.color":"#243244",
        "pdf.fonttype":42,"ps.fonttype":42,
        "svg.fonttype":"none","svg.hashsalt":"new-model-maxwell-20261117",
        "savefig.facecolor":"white",
    })
    fig=plt.figure(figsize=(15.6,8.8),facecolor="white")
    ax=fig.add_axes([.31,.245,.405,.565])
    ap=fig.add_axes([.765,.245,.195,.565],sharey=ax)
    y=np.arange(len(rows))
    labels=[label(r) for r in rows]
    for i,row in enumerate(rows):
        color=CATEGORY[row["category"]][1]
        for target in (ax,ap):
            if i%2==0:target.axhspan(i-.5,i+.5,color="#F4F6F8",zorder=0)
        low,high=row["n_low"],row["n_high"]
        ax.plot([low,high],[i,i],color=color,linewidth=3.1,solid_capstyle="round",zorder=3)
        ax.scatter([low],[i],s=49,facecolors="white",edgecolors=color,linewidths=1.8,zorder=4)
        ax.scatter([high],[i],s=49,color=color,edgecolors="white",linewidths=.7,zorder=4)
        ax.annotate(f"{low:.5f}",(low,i),xytext=(-9,0),textcoords="offset points",
                    ha="right",va="center",color="#334155",fontsize=9.4)
        ax.annotate(f"{high:.5f}",(high,i),xytext=(9,0),textcoords="offset points",
                    ha="left",va="center",color="#334155",fontsize=9.4)
        pressure=row["pressure_transition"]
        ap.scatter([pressure],[i],s=52,color=color,edgecolors="white",linewidths=.7,zorder=4)
        ap.annotate(f"{pressure:.5f}",(pressure,i),xytext=(9,0),textcoords="offset points",
                    ha="left",va="center",fontsize=9.4,color="#334155")
    boundaries=[i-.5 for i in range(1,len(rows)) if rows[i]["category"]!=rows[i-1]["category"]]
    for boundary in boundaries:
        for target in (ax,ap):target.axhline(boundary,color="#CBD2DB",linewidth=1)
    ax.set_yticks(y,labels=labels,fontsize=10.8)
    ax.tick_params(axis="y",length=0,pad=17)
    ap.tick_params(axis="y",left=False,labelleft=False)
    ax.set_ylim(len(rows)-.45,-.55)
    ax.set_xlim(0,.0665)
    ax.xaxis.set_major_locator(FixedLocator(np.arange(0,.061,.01)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f"{x:.2f}"))
    ap.set_xscale("log");ap.set_xlim(.01,.24)
    ap.xaxis.set_major_locator(FixedLocator([.01,.02,.05,.1,.2]))
    ap.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f"{x:.2f}"))
    ap.minorticks_off()
    for target in (ax,ap):
        target.grid(axis="x",color="#D9DFE5",linewidth=.7,zorder=1)
        target.spines["left"].set_visible(False)
        target.set_axisbelow(True)
        target.tick_params(axis="x",length=3.5,pad=7)
    ax.set_xlabel(r"Baryon density $n_B$ [fm$^{-3}$]",labelpad=11)
    ap.set_xlabel(r"$P_t$ [MeV fm$^{-3}$]  (log scale)",labelpad=11)
    ax.set_title("Density discontinuity",loc="left",pad=14,fontsize=12,weight="bold")
    ap.set_title("Transition pressure",loc="left",pad=14,fontsize=12,weight="bold")
    fig.text(.045,.952,"Nine finite-domain sharp-Maxwell candidates",fontsize=22,weight="bold",color="#142536")
    fig.text(.045,.911,
        f"{campaign['completed_cells']} scoped model/variant pairs  |  "
        "6 native pairings, 2 reference-assisted hybrids, 1 public pairing",
        fontsize=11.4,color="#475569")
    handles=[]
    for category,(order,color,title) in CATEGORY.items():
        count=sum(r["category"]==category for r in rows)
        handles.append(Line2D([0],[0],color=color,marker="o",lw=2.5,
                              label=f"{title} ({count})"))
    fig.legend(handles=handles,loc="upper left",bbox_to_anchor=(.039,.879),
               ncols=3,frameon=False,fontsize=10.4,columnspacing=2.9,handlelength=2.1)
    fig.text(.045,.817,"Crust source  "+r"$\rightarrow$"+"  core source",
             fontsize=11.1,weight="bold",color="#243244")
    fig.text(.045,.162,
        "Open endpoint: "+r"$n_{\rm low}$"+"; filled endpoint: "+r"$n_{\rm high}$"+
        ". Right-hand dots show the recorded "+r"$P_t$"+
        ". Exact extracted values are supplied in the CSV.",
        fontsize=10,color="#334155")
    fig.text(.045,.119,
        "Scope: admitted finite source components; thermodynamics only, with no abundance matching. "
        "Global ground state is not claimed.",
        fontsize=10,color="#475569")
    fig.text(.045,.082,
        f"Campaign native inputs: {campaign['physical_model_count_new_native']} new physical RMF models / "
        f"{campaign['native_variant_count']} frozen parameter variants. Nine outputs count source pairings.",
        fontsize=10,color="#475569")
    fig.text(.045,.045,
        "Public tables retain their printed source precision. All plotted numbers are copied from "
        "CAMPAIGN.json native_report; no EOS physics is recomputed.",
        fontsize=9.4,color="#64748B")
    return fig

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign",required=True,type=Path)
    parser.add_argument("--output-dir",required=True,type=Path)
    parser.add_argument("--expected-successes",default=9,type=int)
    args=parser.parse_args()
    source=args.campaign.resolve();out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    campaign=json.loads(source.read_text())
    rows=extract(campaign)
    if len(rows)!=args.expected_successes:
        raise ValueError(f"Expected {args.expected_successes} recorded successes, found {len(rows)}")
    csv_path=out/"new_model_maxwell_successes.csv"
    with csv_path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=CSV_FIELDS)
        writer.writeheader();writer.writerows(rows)
    fig=make_figure(rows,campaign)
    moment=datetime.fromisoformat(campaign["created_at_utc"])
    stem=out/"new_model_maxwell_successes"
    outputs=[]
    for extension in ("png","svg","pdf"):
        target=stem.with_suffix("."+extension)
        if extension=="pdf":
            metadata={"Title":"Nine finite-domain sharp-Maxwell candidates",
                "Author":"ns_cowling project","CreationDate":moment,"ModDate":moment,
                "Subject":"Descriptive figure from recorded native Maxwell reports"}
        else:
            metadata={"Title":"Nine finite-domain sharp-Maxwell candidates",
                      "Date":campaign["created_at_utc"]}
        fig.savefig(target,dpi=220,metadata=metadata)
        outputs.append(target)
    plt.close(fig)
    receipt={
        "schema":"ns_cowling.new_model_maxwell_figure.v1",
        "status":"FIGURE_GENERATED_FROM_RECORDED_VALUES",
        "campaign_path":str(source),"campaign_sha256":sha256(source),
        "script_path":str(Path(__file__).resolve()),"script_sha256":sha256(__file__),
        "command":[sys.executable,*sys.argv],
        "python":sys.version,"matplotlib":matplotlib.__version__,"numpy":np.__version__,
        "platform":platform.platform(),"successes":len(rows),
        "categories":{k:sum(r["category"]==k for r in rows) for k in CATEGORY},
        "numerical_source":"candidate_results[].v2_report.native_report",
        "units":{"n_low":"fm^-3","n_high":"fm^-3","pressure_transition":"MeV fm^-3",
            "mu_transition":"MeV","density_jump":"fm^-3","latent_energy_density":"MeV fm^-3"},
        "physics_recomputed":False,"source_interpolation":False,
        "display_rounding":"density and pressure: 5 decimal places; CSV stores round-trip float strings",
        "caption_scope":{"finite_components_only":True,"thermo_only":True,
            "new_physical_models":campaign["physical_model_count_new_native"],
            "frozen_native_parameter_variants":campaign["native_variant_count"],
            "global_ground_state":campaign["global_ground_state"]},
        "outputs":{p.name:{"sha256":sha256(p),"bytes":p.stat().st_size} for p in [csv_path,*outputs]},
    }
    (out/"FIGURE_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"successes":len(rows),"outputs":[str(p) for p in outputs],
                      "source_sha256":receipt["campaign_sha256"]},indent=2))
if __name__=="__main__":main()
