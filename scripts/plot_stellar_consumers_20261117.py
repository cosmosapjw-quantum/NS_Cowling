#!/usr/bin/env python3
"""Render stage-D measured sequence and eigenfunctions; no new physics solve."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args()
    root=args.root.resolve()
    data=root/"data/stellar_consumers_20261117"
    out=data/"figures"
    out.mkdir(parents=True,exist_ok=True)
    seq_path=data/"tov/sequence.csv"
    c_path=data/"tov/CANONICAL_1P4.json"
    seq=np.genfromtxt(seq_path,delimiter=",",names=True)
    canonical=json.loads(c_path.read_text())
    if "canonical" in canonical: canonical=canonical["canonical"]
    plt.rcParams.update({"font.size":10,"axes.grid":True,"grid.alpha":.22,
                         "savefig.dpi":180,"font.family":"DejaVu Sans"})
    fig,axes=plt.subplots(1,2,figsize=(11.2,4.5),layout="constrained")
    blue="#205c9c"; orange="#bd5a28"
    axes[0].plot(seq["radius_km"],seq["mass_msun"],"o-",ms=3,color=blue,lw=1.5)
    axes[0].scatter([canonical["radius_km"]],[canonical["mass_msun"]],s=65,
                    marker="*",color=orange,zorder=3,label=r"$1.4\,M_\odot$")
    axes[0].set(xlabel="Radius [km]",ylabel=r"Gravitational mass [$M_\odot$]",
                title="Exact B-contract TOV sequence")
    axes[0].legend(loc="upper right",frameon=False)
    axes[0].annotate(f'R = {canonical["radius_km"]:.5f} km',
                    (canonical["radius_km"],canonical["mass_msun"]),
                    xytext=(16.2,1.40),arrowprops={"arrowstyle":"-","color":orange},color=orange)
    axes[1].plot(seq["central_density_fm3"],seq["mass_msun"],"o-",ms=3,color=blue,lw=1.5)
    axes[1].scatter([1.],[seq["mass_msun"][-1]],color=orange,s=30,zorder=3)
    axes[1].annotate("Actual EOS support ends here\nMaximum mass is not bracketed",
                    (1.,seq["mass_msun"][-1]),xytext=(.31,1.68),
                    arrowprops={"arrowstyle":"->","color":orange},color=orange,fontsize=9)
    axes[1].set(xlabel=r"Central baryon density [fm$^{-3}$]",
                ylabel=r"Gravitational mass [$M_\odot$]",title="37 measured stellar models")
    fig.suptitle("Reference-assisted BSk22 + SPG + SLy4 hybrid | cold equilibrium fluid",fontsize=12)
    for ext in ("png","svg"):fig.savefig(out/f"TOV_SEQUENCE.{ext}")
    plt.close(fig)
    inputs=[seq_path,c_path]
    paths=[data/"cowling/ode_tight_f_profile.csv",data/"cowling/ode_tight_p1_profile.csv"]
    results=data/"cowling/ode_tight.json"
    if all(p.exists() for p in paths) and results.exists():
        fig,axes=plt.subplots(1,2,figsize=(11.2,4.3),layout="constrained")
        for p,label,color in zip(paths,["f: 0 radial nodes","p1: 1 radial node"],[blue,orange]):
            a=np.genfromtxt(p,delimiter=",",names=True)
            for ax,field in zip(axes,["W","V_km_inverse"]):
                ax.plot(a["r_km"],a[field],label=label,color=color,lw=1.5)
                ax.axhline(0,color="black",lw=.5)
                ax.set_xlabel("Radius [km]")
                ax.legend(frameon=False,fontsize=9)
        axes[0].set_ylabel("W / max(|W|)")
        axes[1].set_ylabel(r"V / max(|W|) [km$^{-1}$]")
        axes[0].set_title("Radial displacement variable")
        axes[1].set_title("Tangential displacement variable")
        fig.suptitle(r"$\ell=2$ equilibrium-barotropic Cowling | tight ODE control",fontsize=12)
        for ext in ("png","svg"):fig.savefig(out/f"COWLING_MODES.{ext}")
        plt.close(fig)
        inputs+=paths+[results]
    receipt={"status":"RENDERED_FROM_MEASURED_OUTPUTS",
             "inputs":{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
             "outputs":{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.glob("*")) if p.suffix in (".png",".svg")},
             "plot_scope":"TOV whole measured domain; Cowling ode_tight configuration profiles, not a frozen-composition calculation."}
    (out/"PROVENANCE.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":receipt["status"],"outputs":list(receipt["outputs"])}))
if __name__=="__main__":
    main()
