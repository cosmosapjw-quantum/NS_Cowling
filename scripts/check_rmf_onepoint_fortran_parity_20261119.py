#!/usr/bin/env python3
"""Independent source-locked Fortran term-by-term parity at local field solutions.

Fortran gets primitive physical state variables, not Python energy or pressure
components; both independently calculate kinetic integrals/meson potential
pieces and satisfy equations of motion. Original GPPVA provider has no fields.
"""
from __future__ import annotations
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import gppva_rmf
from compose_eos.rmf_onepoint_components import nonlinear_component_budget
DATA=ROOT/"data/cowling_original_rmf_onepoint_20261119"
BIN=ROOT/"build/rmf_onepoint_native/rmf_onepoint_components"
SRC=[
"fortran_eos/src/compose_eos_types.f90",
"fortran_eos/src/compose_eos_tf_models.f90",
"fortran_eos/app/rmf_original_onepoint_components.f90",
]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def build():
    d=BIN.parent
    d.mkdir(parents=True,exist_ok=True)
    cmd=["gfortran","-std=f2018","-O0","-Wall","-Wextra","-fcheck=all",
         "-ffree-line-length-none","-J",str(d),"-I",str(d),
         *SRC,"-o",str(BIN)]
    x=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True)
    (DATA/"FORTRAN_COMPONENT_BUILD.log").write_text(x.stdout+x.stderr)
    if x.returncode:raise RuntimeError("Fortran component compile failed "+x.stderr[-1300:])
    return dict(command=cmd,compiler=subprocess.check_output(["gfortran","--version"],text=True).splitlines()[0],
       compiler_log_sha256=sha(DATA/"FORTRAN_COMPONENT_BUILD.log"),
       source_sha256={name:sha(ROOT/name) for name in SRC},
       compiled_with_bounds_check=True)

def compare(model,b,tag):
    n=b["n_B_fm3"]
    line=(model+"\n"+
        " ".join(f"{v:.17g}" for v in
        (n,b["Yp"],n*b["Ye"],n*b["Ymu"],b["local_phi_positive_MeV"],
         b["local_omega0_MeV"],b["local_rho30_MeV"]))+"\n")
    r=subprocess.run([str(BIN)],input=line,text=True,capture_output=True)
    if r.returncode:raise RuntimeError("Fortran component execution failed "+r.stderr[-1300:])
    f={}
    for record in r.stdout.splitlines():
        tokens=record.split()
        if len(tokens)==2 and tokens[0] not in {"STOP"}:
            f[tokens[0]]=float(tokens[1])
    if len(f)<29:raise RuntimeError(f"insufficient Fortran component output {len(f)}")
    fields={
        "epsilon_total":"epsilon_total_MeV_fm3",
        "pressure_total":"pressure_total_MeV_fm3",
        "Mstar_n":"mstar_MeV",
        "local_phi":"local_phi_positive_MeV",
        "local_omega0":"local_omega0_MeV",
        "local_rho30":"local_rho30_MeV",
    }
    fortran_minus_python={key:f[key]-b[original] for key,original in fields.items()}
    paired={}
    suffix={"neutrons":"neutron","protons":"proton",
            "electrons":"electron","muons":"muon"}
    for typ,source,pre in [
        ("epsilon",b["epsilon_components_MeV_fm3"],"energy"),
        ("pressure",b["pressure_components_MeV_fm3"],"pressure")]:
        for name,val in source.items():
            if name=="pressure_rearrangement":continue
            if name.startswith("kinetic_"):
                fkey=f"{typ}_kinetic_{suffix[name.removeprefix('kinetic_')]}"
            else:
                fkey=f"{typ}_{name.removeprefix(pre+'_')}"
            if fkey not in f:raise RuntimeError("Missing F90 term "+fkey)
            paired[fkey]=f[fkey]-val
    EOM={t:f[t+"_EOM_fm3"] for t in ("scalar","omega","rho")}
    if max(abs(v) for v in [*fortran_minus_python.values(),*paired.values()])>1e-8:
        raise ValueError("Fortran/Python direct component parity failed "+tag)
    if max(abs(v) for v in EOM.values())>1e-7:
        raise ValueError("Fortran EOM residual failed "+tag)
    return dict(status="INDEPENDENT_FORTRAN_PYTHON_TERM_PARITY_PASSED",
       model=model,case=tag,n_B_fm3=n,
       maximum_component_delta_MeV_fm3=float(max(abs(v) for v in paired.values())),
       maximum_total_and_field_delta=float(max(abs(v) for v in fortran_minus_python.values())),
       fortran_minus_python=fortran_minus_python,
       termwise_fortran_minus_python=paired,Fortran_EOM_fm3=EOM,
       local_fields_not_original_provider=True,fortran_raw_values=f)

def main():
    build_info=build()
    z=json.loads((DATA/"RESULT.json").read_text())
    rows={}
    for name in ("FSU2","TM1e"):
        entry=z["models"][name]
        for tag,key in (
            ("local_beta_equilibrium","local_mean_field_components_at_local_beta_equilibrium"),
            ("original_printed_composition","local_mean_field_components_at_original_printed_Yp_Ye_Ymu")):
            rows[name+"_"+tag]=compare(name,entry[key],tag)
    rows["TM1e_shen_original_author_fixedYp_baryons_only"]=compare(
       "TM1e",z["independent_author_TM1e"]["local_fixed_Yp_baryons_only"],
        "Shen2020_original_fixed_Yp_0.10_T0_uniform_baryons_only")
    overall=dict(schema="ns_cowling.original_rmf_20261119_fortran_python_component_parity.v1",
        source_provider_comparison_sha256=sha(DATA/"RESULT.json"),
        build=build_info,independent_cases=rows,
        no_provider_meson_field_data=True,
        status="ALL_FIVE_LOCAL_STATES_INDEPENDENT_NATIVE_F90_COMPONENT_PARITY_VERIFIED")
    (DATA/"FORTRAN_COMPONENT_PARITY.json").write_text(json.dumps(overall,indent=2)+"\n")
    print("ALL_FORTRAN_COMPONENTS_PARITY_PASSED",len(rows),
        "max_abs",max(x["maximum_component_delta_MeV_fm3"] for x in rows.values()))
if __name__=="__main__":main()
