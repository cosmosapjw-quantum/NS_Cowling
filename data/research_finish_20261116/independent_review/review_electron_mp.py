"""Independent high-precision electron-screening oracle for retained counterexamples."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.spatial_tf import RadialGrid,E2,HBARC,ELECTRON_MASS
from compose_eos.spatial_tf_equilibrium import screened_electrons
mp.mp.dps=70
M=mp.mpf(str(ELECTRON_MASS));H=mp.mpf(str(HBARC));K=4*mp.pi*mp.mpf(str(E2))
def density(q):
    if q<=M:return mp.mpf(0)
    return ((q-M)*(q+M))**mp.mpf("1.5")/(3*mp.pi**2*H**3)
def susceptibility(q):
    if q<=M:return mp.mpf(0)
    return q*mp.sqrt((q-M)*(q+M))/(mp.pi**2*H**3)
def solve_mp(case):
    N=case["points"];d=case["dimension"];R=mp.mpf(str(case["radius"]))
    p=[mp.mpf(str(x)) for x in case["protons"]]
    w=[(mp.mpf(i+1)/N)**d-(mp.mpf(i)/N)**d for i in range(N)]
    avg=mp.fsum(wi*pi for wi,pi in zip(w,p))
    def mean(v):return mp.fsum(wi*vi for wi,vi in zip(w,v))
    def poisson(charge):
        net=mean(charge);acc=mp.mpf(0);v=[mp.mpf(0)]
        for i in range(N-1):
            acc+=w[i]*(charge[i]-net)
            conductance=d*(mp.mpf(i+1)/N)**(d-1)*N
            v.append(v[-1]-K*R**2*acc/conductance)
        gauge=mean(v);return [vi-gauge for vi in v]
    def neutral(V,guess):
        lo=M-max(V);hi=mp.sqrt(M*M+H*H*(3*mp.pi**2*avg)**(mp.mpf(2)/3))-min(V)
        mu=min(max(guess,lo),hi)
        for _ in range(120):
            q=[mu+v for v in V];n=[density(x) for x in q]
            C=mean(n)-avg
            if abs(C)<avg*mp.mpf("1e-55"):return mu,n
            if C>0:hi=mu
            else:lo=mu
            derivative=mean([susceptibility(x) for x in q])
            trial=mu-C/derivative if derivative else (lo+hi)/2
            mu=trial if lo<trial<hi else (lo+hi)/2
        raise RuntimeError("independent scalar neutrality did not converge")
    mu=mp.sqrt(M*M+H*H*(3*mp.pi**2*avg)**(mp.mpf(2)/3))
    V=poisson([pi-avg for pi in p])
    for iteration in range(60):
        mu,n=neutral(V,mu)
        newV=poisson([pi-ni for pi,ni in zip(p,n)])
        delta=max(abs(a-b) for a,b in zip(V,newV))
        V=newV
        if delta<mp.mpf("1e-55"):break
    else:raise RuntimeError("independent Poisson fixed point did not converge")
    mu,n=neutral(V,mu)
    return np.array(list(map(float,n))),np.array(list(map(float,V))),float(mu),{
        "iterations":iteration+1,"mean_protons":float(avg),
        "relative_charge_residual":float(abs(mean(n)-avg)/avg),
        "potential_iteration_change":float(delta)}
results=[]
module=ROOT/"python_frontend/compose_eos/spatial_tf_equilibrium.py"
before=hashlib.sha256(module.read_bytes()).hexdigest()
for name in ("warm","cold"):
    path=ROOT/f"data/research_finish_20261116/electron_repair/{name}_failure_input.json"
    case=json.loads(path.read_text());nref,vref,muref,oracle=solve_mp(case)
    grid=RadialGrid(case["points"],case["dimension"]);p=np.array(case["protons"])
    n,v,mu,error=screened_electrons(grid,p,case["radius"],(np.array(case["initial_v"]),case["initial_mu"]))
    results.append({"case":name,"input_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
        "oracle":oracle,"ne_max_abs":float(np.max(abs(n-nref))),
        "ne_max_relative":float(np.max(abs(n-nref))/max(nref)),
        "V_max_abs_MeV":float(np.max(abs(v-vref))),"mu_abs_MeV":abs(mu-muref),
        "production_relative_neutrality":abs(grid.mean(n-p))/grid.mean(p),
        "production_PDE_max":error,"gauge_max":abs(grid.mean(v))})
after=hashlib.sha256(module.read_bytes()).hexdigest()
out={"method":"70-digit scalar neutrality plus independent finite-volume flux-integral Poisson fixed point; no production Newton or eigen-Poisson inverse reused",
     "functional_units":"MeV, fm; physical constants retained explicitly","module_sha256":before,"module_unchanged":before==after,
     "runner_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"cases":results}
Path(__file__).with_name("ELECTRON_ORACLE.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
print(json.dumps(out,indent=2))
