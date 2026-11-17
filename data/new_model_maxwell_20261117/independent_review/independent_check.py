#!/usr/bin/env python3
"""C9 bounded independent reference: source-sign fields and momentum quadrature.
No imports from compose_eos or native RMF kernels; no production state is changed.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

ROOT = Path("/home/cosmosapjw/ns_cowling_wt/new_model_maxwell_20261117")
DATA = ROOT / "data/new_model_maxwell_20261117"
OUT = DATA / "independent_review"
OUT.mkdir(exist_ok=True)
PARAM = DATA / "sources/PARAMETER_INPUTS.json"
RAW = json.loads(PARAM.read_text())
HC = 197.3269804
ME, MMU = .51099895, 105.6583755
PI = np.pi
T, W = leggauss(192)
T, W = (T+1)/2, W/2

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def parameters(name):
    p = dict(RAW["models"][name])
    return {**parameters(p["inherits"]), **p} if "inherits" in p else p

def couplings(p,n):
    if p["kind"] != "density_dependent":
        return np.array([p["gs"],p["gw"],p["Gamma_rho"]])
    x=n/p["n0"]
    vals=[]
    for pref,base in [("s","gs"),("w","gw")]:
        a,b,c,d = [p[pref+k] for k in ["a","b","c","d"]]
        vals.append(p[base]*a*(1+b*(x+d)**2)/(1+c*(x+d)**2))
    vals.append(p["Gamma_rho"]*np.exp(-p["ra"]*(x-1)))
    return np.array(vals)

def fermion(n,m):
    m=np.asarray(m)
    k=HC*(3*PI**2*n)**(1/3) if n>0 else 0.
    en=np.sqrt(k*k*T*T+m[...,None]**2)
    pref=k**3/(PI**2*HC**3)
    eps=pref*np.sum(W*T*T*en,axis=-1)
    scalar=pref*m*np.sum(W*T*T/en,axis=-1)
    pres=pref*k*k/3*np.sum(W*T**4/en,axis=-1)
    deriv=pref*k*k*np.sum(W*T**4/en**3,axis=-1)
    return eps,scalar,pres,deriv,np.hypot(k,m)

def source_state(p,n,y,scan=False):
    gs,gw,G=couplings(p,n)
    dgs,dgw,dG=np.imag(couplings(p,n+1e-20j))/1e-20
    nn,np_=n*(1-y),n*y
    g2=p.get("g2",0)*HC;g3=p.get("g3",0);c3=p.get("c3",0)
    def terms(sig):
        mn=p["mn"]+gs*sig;mp=p["mp"]+gs*sig
        fn=fermion(nn,mn);fp=fermion(np_,mp)
        u=.5*p["ms"]**2*sig**2+g2*sig**3/3+g3*sig**4/4
        f=(p["ms"]**2*sig+g2*sig**2+g3*sig**3)/HC**3+gs*(fn[1]+fp[1])
        cur=(p["ms"]**2+2*g2*sig+3*g3*sig**2)/HC**3+gs**2*(fn[3]+fp[3])
        return f,cur,fn[0]+fp[0]+u/HC**3
    cap=min(p["mn"],p["mp"])/gs*(1-1e-10)
    if scan:
        grid=np.linspace(-cap,0,4097);f=terms(grid)[0]
        ii=np.where((f[:-1]<=0)&(f[1:]>=0))[0]
        roots=[brentq(lambda s:terms(s)[0],grid[i],grid[i+1],xtol=2e-13) for i in ii]
        roots=[s for s in roots if terms(s)[1]>0]
    else:
        roots=[brentq(lambda s:terms(s)[0],-cap,0,xtol=2e-13)]
    if not roots:raise RuntimeError("independent scan found no stable scalar root")
    sig=min(roots,key=lambda s:terms(s)[2])
    w0=gw*n*HC**3/p["mw"]**2
    omega=brentq(lambda w:p["mw"]**2*w+c3*w**3-gw*n*HC**3,0,w0,xtol=2e-13) if c3 else w0
    rho=G*(nn-np_)*HC**3/p["mr"]**2
    fn=fermion(nn,p["mn"]+gs*sig);fp=fermion(np_,p["mp"]+gs*sig)
    u=.5*p["ms"]**2*sig**2+g2*sig**3/3+g3*sig**4/4
    sr=dgw*omega*n+dG*rho*(nn-np_)+dgs*sig*(fn[1]+fp[1])
    en=fn[0]+fp[0]+(u+.5*p["mw"]**2*omega**2+.75*c3*omega**4+.5*p["mr"]**2*rho**2)/HC**3
    pr=fn[2]+fp[2]+(-u+.5*p["mw"]**2*omega**2+.25*c3*omega**4+.5*p["mr"]**2*rho**2)/HC**3+n*sr
    mun=fn[4]+gw*omega+G*rho+sr
    mup=fp[4]+gw*omega-G*rho+sr
    def number(mu,m):return max(mu*mu-m*m,0.)**1.5/(3*PI**2*HC**3)
    upper=np.hypot(HC*(3*PI**2*np_)**(1/3),ME)
    mue=brentq(lambda mu:number(mu,ME)+number(mu,MMU)-np_,ME,upper,xtol=2e-13) if upper>MMU else upper
    ne,nmu=number(mue,ME),number(mue,MMU)
    fe=fermion(ne,ME);fm=fermion(nmu,MMU)
    return {"epsilon":float(en+fe[0]+fm[0]),"pressure":float(pr+fe[2]+fm[2]),
            "mu_n":float(mun),"mu_p":float(mup),"mu_e":float(mue),
            "phi":float(-sig),"omega":float(omega),"rho":float(-rho),"rearrangement":float(sr),
            "Y_e":float(ne/n),"Y_mu":float(nmu/n),"beta":float(mun-mup-mue),
            "scalar_curvature":float(terms(sig)[1]),"scan_root_count":len(roots)}

models=["NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL","DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0"]
rows=[];branch_scans=[];provenance=[];grid=[]
maxima={}
for name in models:
    p=parameters(name);folder=DATA/"generated"/name
    fine=np.genfromtxt(folder/"fine/states.csv",delimiter=",",names=True)
    coarse=np.genfromtxt(folder/"coarse/states.csv",delimiter=",",names=True)
    for key in fine.dtype.names:
        assert np.array_equal(coarse[key],fine[key][::2]),(name,key)
    grid.append({"model":name,"coarse_rows":len(coarse),"fine_rows":len(fine),"common_rows_exact":True})
    for level in ["coarse","fine"]:
        d=folder/level
        generation=json.loads((d/"GENERATION.json").read_text())
        run=json.loads((d/"RUN.json").read_text())
        observed={"source":digest(PARAM)==run["source_parameter_sha256"],
                  "parameter":digest(d/"parameters.nml")==run["parameters_nml_sha256"],
                  "density":digest(d/"density.txt")==run["density_sha256"],
                  "executable":digest(run["command"][0])==run["executable_sha256"],
                  "tables":all(digest(d/"compose"/f)==h for f,h in generation["table_sha256"].items()),
                  "run_exit0":run["returncode"]==0}
        assert all(observed.values()),(name,level,observed)
        provenance.append({"model":name,"level":level,"checks":observed,"states_sha256":digest(d/"states.csv")})
    for index in [0,160,512,736,768]:
        actual=fine[index]
        n,y=float(actual["n_b"]),float(actual["Y_p"])
        reference=source_state(p,n,y,scan=True)
        residual={}
        for key in ["epsilon","pressure","mu_n","mu_p","mu_e","phi","omega","rho","rearrangement","Y_e","Y_mu","scalar_curvature"]:
            v=abs(reference[key]-actual[key])/max(1,abs(actual[key]))
            residual[key]=float(v);maxima[key]=max(maxima.get(key,0),float(v))
        dy=1e-5
        derivative=(source_state(p,n,y+dy)["beta"]-source_state(p,n,y-dy)["beta"])/(2*dy)
        assert derivative<0,(name,n,derivative)
        assert reference["scan_root_count"]==int(actual["scalar_root_count"]),(name,n)
        rows.append({"model":name,"fine_index":index,"density":n,"yp":y,"independent_reference":reference,
                     "scaled_residuals":residual,"d_beta_d_yp_mev":derivative,
                     "local_composition_second_derivative_energy_mevfm3":-n*derivative})
    yy=np.linspace(0,1,17)
    bf=[source_state(p,.16,float(y))["beta"] for y in yy]
    assert bf[0]>0 and bf[-1]<0 and np.all(np.diff(bf)<0),(name,bf)
    branch_scans.append({"model":name,"density":.16,"yp_grid":yy.tolist(),"beta_values_mev":bf,
                        "sampled_monotone":True,"scope":"17 sampled compositions; not global convexity proof"})

build=json.loads((DATA/"uniform_validation/build/BUILD.json").read_text())
source_ids={f:digest(f) for f in build["sources"]}
comment_equivalence = {}
for filename, old_hash in build["sources"].items():
    if source_ids[filename] != old_hash:
        current = Path(filename).read_text()
        old = current.replace(
            " ! Resolve all detected upward roots on a 256-subinterval positive-mass scan.\n ! Select the lowest scalar+kinetic energy among detected stable roots only.\n ! Arbitrarily close unresolved root pairs/global minimality are not excluded.\n",
            " ! Isolate every upward root in the positive-mass interval and select the\n ! lowest scalar+kinetic energy among locally stable stationary points.\n")
        assert hashlib.sha256(old.encode()).hexdigest() == old_hash, filename
        comment_equivalence[filename] = {"production_sha256":old_hash,"current_sha256":source_ids[filename],
          "required_relation":"SEMANTIC_IDENTITY","evidence":"Replacing exactly the three explanatory comments with the originally reviewed two comments reconstructs the recorded production byte hash."}
source_ids.update({str(ROOT/p):digest(ROOT/p) for p in ["python_frontend/compose_eos/rmf_candidate.py","scripts/generate_rmf_candidates_20261117.py","python_frontend/tests/test_rmf_candidate.py"]})
ledger=json.loads((DATA/"sources/SOURCE_LEDGER.json").read_text())
source_matches={name: digest(DATA/"sources"/name)==entry["sha256"] for name,entry in ledger["source_byte_identities"].items()}
assert all(source_matches.values()),{k:v for k,v in source_matches.items() if not v}
assert max(maxima.values())<2e-8,maxima
result={"schema":"ns_cowling.independent_uniform_review.v1","status":"NUMERICALLY_CHECKED",
        "selection":"5 frozen inputs, 5 existing fine-grid rows each; one independent check subset",
        "method":"source sigma sign, full Gamma_rho vertex, 192-node momentum quadrature, complex-step DD derivatives, 4096-bin scalar scan; no imports from compose_eos",
        "precision_gate":"scaled field/state residual <2e-8; no tolerance tuning",
        "max_scaled_residuals":maxima,"max_abs_beta_mev":max(abs(x["independent_reference"]["beta"]) for x in rows),
        "checks":rows,"composition_scans":branch_scans,"common_grid_checks":grid,
        "provenance_checks":provenance,"source_ledger_file_count":len(source_matches),"source_ledger_bytes_match":all(source_matches.values()),
        "source_ids":source_ids,"comment_only_source_equivalence":comment_equivalence,"parameter_sha256":digest(PARAM),
        "runtime":{"python":__import__("sys").version,"numpy":np.__version__,"scipy":__import__("scipy").__version__},
        "scope_limit":"Finite selected states and scans; no global microscopic minimum or complete-continuum root uniqueness claim.",
        "conclusion":"Reference reconstruction supports native uniform thermodynamics and selected physical branches."}
(OUT/"UNIFORM_SUBSET.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
print(json.dumps({k:result[k] for k in ["status","selection","max_scaled_residuals","max_abs_beta_mev","source_ledger_file_count","source_ledger_bytes_match","conclusion"]},indent=2))
