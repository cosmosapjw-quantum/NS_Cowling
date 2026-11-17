"""Explicit-parameter native uniform RMF driver and lossless CompOSE adapter.

All field and thermodynamic evaluations are Fortran calls into the unchanged
TF model kernels. Python handles source binding, subprocesses and serialization.
No provider residual is used to modify a parameter or physical energy.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import os
import subprocess
import tempfile
import numpy as np
from .eos_data import EOSData
from .api import write_compose, read_compose

HBARC = 197.3269804
FIELDS = ("mn mp n0 ms mw mr gs gw gr sa sb sc sd wa wb wc wd ra a3 a4 b4 cross").split()
PRIMARY_MODELS = ("NL3_PROVIDER_PRECISION", "PK1_CANONICAL",
                  "TM1_CANONICAL", "DDMEX_PROVIDER_ROUNDED_N0")
DIAGNOSTIC_MODELS = ("DDMEX_ORIGINAL_TABLE_RHO0",)

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_parameters(source_json, model_id):
    source=Path(source_json).resolve()
    raw=json.loads(source.read_text())
    def resolve(name):
        p=dict(raw["models"][name])
        if "inherits" in p:
            p={**resolve(p["inherits"]),**p}
        return p
    p=resolve(model_id)
    if "g2" in p:
        p.update(a3=-float(p["g2"])*HBARC,a4=float(p["g3"]),b4=float(p["c3"]))
    p.update(id=model_id,source_parameter_file=str(source),source_parameter_sha256=sha256(source))
    p["density_dependent"]=p["kind"]=="density_dependent"
    return p

def write_parameters(path, parameters):
    p=parameters
    lines=["&candidate",f" p%id='{p['id']}'",
           " p%density_dependent="+(".true." if p["density_dependent"] else ".false.")]
    for key in FIELDS:
        value=float(p.get(key,0.))
        if not np.isfinite(value):
            raise ValueError("parameters must be finite")
        lines.append(f" p%{key}={value:.17e}")
    lines.append("/")
    Path(path).write_text("\n".join(lines)+"\n")

def build_driver(repo, build_dir):
    repo=Path(repo).resolve();build_dir=Path(build_dir).resolve()
    build_dir.mkdir(parents=True,exist_ok=True)
    sources=[repo/"fortran_eos/src/compose_eos_types.f90",
             repo/"fortran_eos/src/compose_eos_tf_models.f90",
             repo/"fortran_eos/tools/rmf_candidate_driver.f90"]
    exe=build_dir/"rmf_candidate_driver"
    command=["gfortran","-std=f2008","-O2","-fcheck=all","-fbacktrace",
             "-ffree-line-length-none","-J",str(build_dir),"-I",str(build_dir),
             *(str(s) for s in sources),"-o",str(exe)]
    proc=subprocess.run(command,cwd=build_dir,text=True,capture_output=True)
    (build_dir/"BUILD.stdout.log").write_text(proc.stdout)
    (build_dir/"BUILD.stderr.log").write_text(proc.stderr)
    receipt={"command":command,"returncode":proc.returncode,"cwd":str(build_dir),
             "compiler":subprocess.check_output(["gfortran","--version"],text=True).splitlines()[0],
             "sources":{str(s):sha256(s) for s in sources}}
    if proc.returncode==0:receipt["executable_sha256"]=sha256(exe)
    (build_dir/"BUILD.json").write_text(json.dumps(receipt,indent=2)+"\n")
    if proc.returncode:
        raise RuntimeError(f"native RMF build failed ({proc.returncode}): {proc.stderr}")
    return exe

def evaluate_states(executable, parameters, density, *, mode="beta", yp=.5, work_dir=None):
    n=np.atleast_1d(np.asarray(density,dtype=float))
    if n.ndim!=1 or n.size==0 or not np.all(np.isfinite(n)) or np.any(n<=0):
        raise ValueError("density must be a nonempty finite positive list")
    if mode not in ("beta","baryons"):raise ValueError("invalid evaluation mode")
    if not np.isfinite(yp) or not 0<=yp<=1:raise ValueError("Yp must be in [0,1]")
    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix="rmf-candidate-") as temp:
            return evaluate_states(executable,parameters,n,mode=mode,yp=yp,work_dir=temp)
    out=Path(work_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    parameter_file=out/"parameters.nml";density_file=out/"density.txt";csv_file=out/"states.csv"
    write_parameters(parameter_file,parameters)
    np.savetxt(density_file,n,fmt="%.17e")
    exe=Path(executable).resolve()
    command=[str(exe),str(parameter_file),str(density_file),str(csv_file),mode,format(float(yp),".17g")]
    env=dict(os.environ);env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    proc=subprocess.run(command,env=env,capture_output=True,text=True)
    (out/"RUN.stdout.log").write_text(proc.stdout);(out/"RUN.stderr.log").write_text(proc.stderr)
    receipt={"command":command,"returncode":proc.returncode,
             "executable_sha256":sha256(exe),"parameters_nml_sha256":sha256(parameter_file),
             "density_sha256":sha256(density_file),"model":parameters["id"],
             "thread_environment":{k:env[k] for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")},
             "source_parameter_sha256":parameters.get("source_parameter_sha256")}
    (out/"RUN.json").write_text(json.dumps(receipt,indent=2)+"\n")
    if proc.returncode:
        raise RuntimeError(f"native RMF evaluation failed ({proc.returncode}): {proc.stderr}")
    a=np.genfromtxt(csv_file,delimiter=",",names=True)
    result={name:np.atleast_1d(a[name]) for name in a.dtype.names}
    if len(result["n_b"])!=len(n) or any(not np.all(np.isfinite(v)) for v in result.values()):
        raise RuntimeError("native RMF output is incomplete or nonfinite")
    return result

def residual_report(states):
    s=states
    report={key:float(np.max(np.abs(s[key]))) for key in (
        "scalar_residual_fm3","omega_residual_fm3","rho_residual_fm3",
        "beta_residual_mev","charge_residual_fm3","euler_residual_mevfm3")}
    report.update(min_effective_mass_n=float(np.min(s["effective_mass_n"])),
                  min_effective_mass_p=float(np.min(s["effective_mass_p"])),
                  min_scalar_curvature=float(np.min(s["scalar_curvature"])),
                  root_counts=np.unique(s["scalar_root_count"]).astype(int).tolist(),
                  composition_branches=np.unique(s["composition_branch"]).astype(int).tolist())
    report["physical_state_checks"] = bool(
        report["scalar_residual_fm3"]<2e-11 and report["omega_residual_fm3"]<2e-11
        and report["rho_residual_fm3"]<2e-11 and report["beta_residual_mev"]<2e-8
        and report["charge_residual_fm3"]<2e-13 and report["euler_residual_mevfm3"]<2e-8
        and min(report["min_effective_mass_n"],report["min_effective_mass_p"],report["min_scalar_curvature"])>0)
    return report

def export_compose(directory, states, parameters):
    directory=Path(directory).resolve();directory.mkdir(parents=True,exist_ok=True)
    s=states;n=s["n_b"];zero=np.zeros(len(n));one=np.ones(len(n))
    if len(n)<2 or np.any(np.diff(n)<=0):raise ValueError("export requires increasing density")
    mn=float(parameters["mn"]);mp=float(parameters["mp"])
    data=EOSData(nb=n,mn=mn,mp=mp,Q1=s["pressure"]/n,Q2=zero,
        Q3=s["mu_B"]/mn-1,Q4=-s["mu_e"]/mn,Q5=zero,
        Q6=s["epsilon"]/(n*mn)-1,Q7=s["epsilon"]/(n*mn)-1,
        q_extra=zero,Iphase=np.zeros(len(n),dtype=int),Yn=one-s["Y_p"],
        Yp=s["Y_p"],Ye=s["Y_e"],Ymu=s["Y_mu"],Aav=zero,Zav=zero,Yh=zero,
        composition_file_present=True,composition_available=True,beta_equilibrated=True)
    target=directory/"compose";write_compose(target,data);readback=read_compose(target)
    comparisons={"density":(readback.nb,n),"pressure":(readback.p,s["pressure"]),
                 "epsilon":(readback.eps,s["epsilon"]),"mu_B":(readback.mub,s["mu_B"]),
                 "Yp":(readback.Yp,s["Y_p"]),"Ye":(readback.Ye,s["Y_e"]),
                 "Ymu":(readback.Ymu,s["Y_mu"])}
    residuals={key:float(np.max(np.abs(a-b)/np.maximum(1.,np.abs(b)))) for key,(a,b) in comparisons.items()}
    if max(residuals.values())>5e-13:raise RuntimeError(f"native readback mismatch: {residuals}")
    report=residual_report(s)
    if not report["physical_state_checks"]:raise RuntimeError(f"physical state residual failed: {report}")
    receipt={"schema":"ns_cowling.rmf_candidate_generation.v1","model":parameters["id"],
       "status":"IMPLEMENTATION_VERIFIED","source_identity":"UNRESOLVED_UNTIL_PROVIDER_COMPARISON",
       "physical_energy_adjustment":"NONE","rows":len(n),"density_min_fm3":float(n[0]),
       "density_max_fm3":float(n[-1]),"compose_directory":str(target),
       "source_parameters":parameters,"residuals":report,"native_readback_scaled_residuals":residuals,
       "units":{"n":"fm^-3","pressure":"MeV fm^-3","epsilon":"MeV fm^-3","mu":"MeV"},
       "table_sha256":{f.name:sha256(f) for f in sorted(target.glob("eos.*"))}}
    (directory/"GENERATION.json").write_text(json.dumps(receipt,indent=2)+"\n")
    return receipt
