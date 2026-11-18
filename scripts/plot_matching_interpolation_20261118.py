"""Source-grounded illustrative plots of the H interpolation sensitivity experiment."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.whole_crust import WholeCrustBarotrope,CausalExponentialInterval
from matching_interpolation_shapes_20261118 import CausalShapeInterval,SlopeMatchedBernsteinInterval
D=ROOT/"data/matching_interpolation_20261118"
FIG=ROOT/"docs/matching_interpolation_20261118/figures"
FIG.mkdir(parents=True,exist_ok=True)
P=json.loads((ROOT/"data/low_density_crust_20261117/stellar/causal_extension/TOV_PLAN.json").read_text())
B={c["id"]:c for c in json.loads((D/"BERNSTEIN_RESULTS.json").read_text())["cases"]}
S={c["id"]:c for c in json.loads((D/"SUMMARY.json").read_text())["results"]}
primary=[x for x in P["cases"] if x["label"]=="primary"]
names=["DDME2","DDME-X","NL3","NL3wrL55"]
fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
for i,(entry,title) in enumerate(zip(primary,names)):
    c=json.loads((ROOT/entry["contract"]).read_text())
    bar=WholeCrustBarotrope.from_contract(c)
    l,h=c["lower_bridge"]["endpoint_states"]
    v=[l[k] for k in ("n","p","epsilon")]+[h[k] for k in ("n","p","epsilon")]
    outer_cs=float(bar.outer.gamma[-1]*bar.bridge.p_left/(bar.bridge.epsilon_left+bar.bridge.p_left))
    near=bar.upper_intervals[0]
    inner_cs=float(near.number(0)/(near.mu_left*near.susceptibility(0)))
    methods=dict(Exponential=CausalExponentialInterval(*v),
                 Rational=CausalShapeInterval(*v,family="rational"),
                 Power=CausalShapeInterval(*v,family="power"),
                 Bernstein=SlopeMatchedBernsteinInterval(*v,left_cs2=outer_cs,right_cs2=inner_cs))
    bk=methods["Bernstein"];record=B[entry["id"]]
    assert bk.degree==record["degree"]
    bk.bernstein_coefficients=np.array(record["bernstein_coefficients"])
    x=np.linspace(0,1,501)
    ax=axes.flat[i]
    for label in ("Exponential","Rational","Power","Bernstein"):
        b=methods[label]
        c2=b.evaluate_mu(b.mu_left+b.width*x)["cs2_over_c2"]
        ax.semilogy(x,c2,label=label,lw=1.7)
    ax.scatter([0,1],[outer_cs,inner_cs],marker="x",color="black",s=38,zorder=9,
               label="Adjacent EOS one-sided $c_s^2$" if i==0 else None)
    ax.set_title(f"{title}, $n_R={entry['upper_start_density']:.4g}$ fm$^{{-3}}$")
    ax.set_ylim(5e-5,2)
    ax.grid(alpha=0.25)
    if i>=2:ax.set_xlabel("Normalized chemical-potential coordinate $x$")
    if i%2==0:ax.set_ylabel("$c_s^2/c^2$")
handles,labels=axes.flat[0].get_legend_handles_labels()
fig.legend(handles,labels,loc="upper center",ncol=5,bbox_to_anchor=(.5,1.02))
fig.suptitle("Cold crust bridge: same endpoints, different causal interpolation",y=1.07)
fig.tight_layout()
fig.savefig(FIG/"CS2_METHOD_SENSITIVITY.png",dpi=180,bbox_inches="tight")
plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5))
idx=np.arange(4)
for j,(name,key) in enumerate([("Power","dRpow"),("Rational","dRrat"),("Bernstein","dRbern")]):
    ax.bar(idx+(j-1)*.23,[S[e["id"]][key] for e in primary],width=.22,label=name)
ax.axhline(0,color="black",linewidth=.8)
ax.set_xticks(idx,names)
ax.set_ylabel("$\\Delta R$ relative to exponential (mm), $M=1.4 M_\\odot$")
ax.set_title("Finite interpolation sensitivity of stellar radii")
ax.grid(axis="y",alpha=.2)
ax.legend()
fig.tight_layout()
fig.savefig(FIG/"TOV_RADIUS_METHOD_SENSITIVITY.png",dpi=180)
plt.close(fig)
files=[FIG/"CS2_METHOD_SENSITIVITY.png",FIG/"TOV_RADIUS_METHOD_SENSITIVITY.png"]
receipt=dict(source="H strict primary contracts plus SUMMARY/BERNSTEIN_RESULTS",
    engine=f"Matplotlib {matplotlib.__version__}",plots=[dict(path=str(x.relative_to(ROOT)),
    size=x.stat().st_size,sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in files],
    caveat="Plots are finite EOS representation comparisons, not microscopic uncertainties or statistical error bars.")
(D/"FIGURES_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
