#!/usr/bin/env python3
"""Plot actually computed E stellar support and equilibrium derivative."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"data/stellar_extension_20261117"
def main():
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    seq=np.genfromtxt(OUT/"tov/sequence.csv",delimiter=",",names=True)
    old=np.genfromtxt(ROOT/"data/stellar_consumers_20261117/tov/sequence.csv",delimiter=",",names=True)
    result=json.loads((OUT/"tov/TOV_RESULT.json").read_text())
    refined=json.loads((OUT/"tov/MAXIMUM_REFINEMENT.json").read_text())["representative_maximum"]
    x=refined["central_density_fm3"]
    y=refined["mass_msun"]
    causal=result["native_causal_endpoint"]["native_causal_root_fm3"]
    fig,ax=plt.subplots(1,2,figsize=(11.8,4.4),layout="constrained")
    nkey="central_density_fm3"
    ax[0].plot(old[nkey],old["mass_msun"],color="#6b7280",lw=2,label="Preserved D support")
    ax[0].plot(seq[nkey],seq["mass_msun"],"o-",ms=3,color="#176B87",label="Native SLy4 extension")
    ax[0].scatter([x],[y],marker="*",s=125,color="#b44e18",zorder=5,label="Mass turning point")
    ax[0].axvline(causal,color="#a12653",ls=":",label=r"Native $\chi_{eq}=1$")
    ax[0].set(xlabel=r"$n_c$ [fm$^{-3}$]",ylabel=r"$M/M_\odot$",xlim=(.4,1.29),ylim=(1.0,2.10))
    ax[0].legend(frameon=False,fontsize=8,loc="lower right")
    ax[1].plot(old["radius_km"],old["mass_msun"],color="#6b7280",lw=2,label="Preserved D support")
    ax[1].plot(seq["radius_km"],seq["mass_msun"],"o-",ms=3,color="#176B87",label="E extension")
    ax[1].scatter([refined["radius_km"]],[y],marker="*",s=125,color="#b44e18",zorder=5)
    ax[1].set(xlabel=r"$R$ [km]",ylabel=r"$M/M_\odot$",xlim=(9.7,13.2),ylim=(1.0,2.10))
    fig.suptitle("Preserved hybrid crust + unchanged native SLy4 core")
    out=OUT/"figures";out.mkdir(exist_ok=True)
    for ext in ("png","svg"):
        fig.savefig(out/f"SLY4_EXTENSION.{ext}",dpi=190)
    plt.close(fig)
    print("SLY4_EXTENSION.png/svg written")
if __name__=="__main__":main()
