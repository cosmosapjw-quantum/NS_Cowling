#!/usr/bin/env python3
"""E native SLy4 causal-domain extension; no functional or energy convention change.

C1-C8: D identity, original native formulas, finite physical tail, source-grid/
ODE/center convergence. Every native process preserves its actual exit/logs.
The causal cutoff is a domain restriction, never a causal completion.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, platform, subprocess, sys, time
from copy import deepcopy
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import brentq, minimize_scalar

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose, write_compose
from compose_eos.core_extension import extend_whole_crust
from compose_eos.crust_assembly import ROW_FIELDS, select_rows
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_background import ContractEnthalpyEOS, solve_background
from build_whole_crust_20261117 import serialized_check

BASE="300e5a0532dd137cb09738568fca0f35cfcc9282"
SOURCE=ROOT/"data/core_support_20261117/primary_v1"
OLDCORE=ROOT/"data/maxwell_three_families_20261117/native_core_extended_v2/RG_SLY4"
DEFAULT_BUILD=ROOT.parent/"stellar_consumers_20261117/build/compose_eos_fortran"
DEFAULT_OUT=ROOT/"data/stellar_extension_20261117"
CONFIG=dict(r0_km=.001,rtol=2e-10,atol=1e-12,max_step_h=.002,surface_h=0.)
DRIVER="""program native_endpoint
 use compose_eos_types,only:wp
 use compose_eos_core_native,only:core_state_t,solve_core_state
 use compose_eos_core_domain,only:core_sound_speed_analytic,core_sound_speed_fd
 implicit none
 type(core_state_t)::s
 real(wp)::n,rel,cs,frozen,hyy,beta,fd,cm,fl
 character(256)::arg,status
 call get_command_argument(1,arg);read(arg,*)n
 call get_command_argument(2,arg);read(arg,*)rel
 call solve_core_state('RG_SLY4',n,s,status)
 if(status/='ok')then;write(*,*)trim(status);stop 1;endif
 call core_sound_speed_analytic('RG_SLY4',n,cs,frozen,hyy,beta,status)
 if(status/='ok')then;write(*,*)trim(status);stop 2;endif
 call core_sound_speed_fd('RG_SLY4',n,rel,fd,cm,fl,status)
 if(status/='ok')then;write(*,*)trim(status);stop 3;endif
 write(*,'(*(es25.16e3,:,","))')n,s%epsilon,s%pressure,s%mu_n,cs,fd,cm,fl,frozen,hyy,beta
end program
"""
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def rows(path,values):
    with Path(path).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
def relative(path):
    try:return str(Path(path).relative_to(ROOT))
    except ValueError:return str(path)
def command_log(command,path):
    p=subprocess.run([str(x) for x in command],capture_output=True,text=True,
        env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1"))
    Path(str(path)+".stdout.log").write_text(p.stdout)
    Path(str(path)+".stderr.log").write_text(p.stderr)
    rec=dict(command=[str(x) for x in command],exit_code=p.returncode,
        stdout_sha256=sha(str(path)+".stdout.log"),stderr_sha256=sha(str(path)+".stderr.log"))
    dump(str(path)+".json",rec)
    if p.returncode:raise RuntimeError(f"native command failed ({p.returncode}): {command}")
    return p,rec

def binding(build):
    paths=[SOURCE/x for x in ("REPRESENTATION.json","ROW_PROVENANCE.json","RESULT.json")]
    paths += [OLDCORE/"DISPATCH.json",OLDCORE/"RECEIPT.json",
        ROOT/"data/stellar_consumers_20261117/NEXT_STEP_PROPOSED_E.json",
        ROOT/"python_frontend/compose_eos/stellar_background.py",
        ROOT/"python_frontend/compose_eos/core_extension.py",
        ROOT/"python_frontend/compose_eos/whole_crust.py"]
    paths += [p for d in (SOURCE/"compose",OLDCORE/"generated") for p in sorted(d.iterdir()) if p.is_file()]
    paths += sorted((ROOT/"fortran_eos/src").glob("*.f90"))
    paths += sorted((ROOT/"fortran_eos/app").glob("*.f90"))
    bindings={}
    for p in paths:
        rel=relative(p)
        original=subprocess.check_output(["git","show",BASE+":"+rel],cwd=ROOT)
        if sha(p)!=hashlib.sha256(original).hexdigest():raise ValueError("D source changed: "+rel)
        bindings[rel]=sha(p)
    # Build objects come from D and carry the same source code; bind the actual
    # linked modules and generator executable separately from source identities.
    binary={str(p):sha(p) for p in [build/"compose_eos_generate"]+
        [build/"obj"/(name+".o") for name in ("compose_eos_types","compose_eos_core_native","compose_eos_core_domain")]}
    return dict(base_commit=BASE,sha256=bindings,binary_sha256=binary,
        native_definition=json.loads((OLDCORE/"DISPATCH.json").read_text())["request"]["contract"],
        energy_reference_shift="NONE",tail_physics="UNCHANGED_NATIVE_SLY4_BETA_EQUILIBRIUM")

def make_eos(out,build):
    eosout=out/"eos"
    if (eosout/"EOS_RESULT.json").exists():return json.loads((eosout/"EOS_RESULT.json").read_text())
    if eosout.exists():raise FileExistsError("partial EOS attempt exists; preserve and use a new output")
    eosout.mkdir(parents=True)
    dump(eosout/"CHECKPOINT.json",dict(stage="E_C1_C8",status="NATIVE_EXTENSION_RUNNING",
        objective="first native causal boundary or causal TOV turning point",base_commit=BASE))
    dump(eosout/"SOURCE_BINDING.json",binding(build))
    driver=eosout/"native_endpoint_driver.f90";driver.write_text(DRIVER)
    binary=eosout/"native_endpoint_driver"
    objs=[build/"obj"/(name+".o") for name in ("compose_eos_types","compose_eos_core_native","compose_eos_core_domain")]
    command_log(["gfortran","-std=f2018","-O2","-ffree-line-length-none","-fcheck=all",
        "-I"+str(build/"mod"),driver,*objs,"-o",binary],eosout/"driver_compile")
    evaluations=[]
    def evaluate(n,rel=2e-5):
        p=subprocess.run([str(binary),format(n,".17g"),format(rel,".17g")],capture_output=True,text=True)
        rec=dict(n_fm3=float(n),relative_fd_step=rel,exit_code=p.returncode,
            stdout=p.stdout.strip(),stderr=p.stderr)
        evaluations.append(rec)
        if p.returncode:raise RuntimeError("native EOS point failed: "+p.stdout+p.stderr)
        vals=[float(v) for v in p.stdout.strip().split(",")]
        keys=("n","epsilon","p","mu","cs2","cs2_fd","cs2_mu_fd","firstlaw_relative","cs2_frozen","composition_hessian","beta_residual")
        return dict(zip(keys,vals))
    scan=[evaluate(n) for n in np.linspace(1.,1.3999,81)]
    rows(eosout/"native_domain_scan.csv",scan)
    for s in scan:
        if min(s["p"],s["epsilon"],s["mu"],s["cs2"],s["composition_hessian"])<=0:
            raise ValueError("mechanical/composition physical limit precedes causal analysis")
    first=next(i for i,s in enumerate(scan) if s["cs2"]>=1.)
    left,right=scan[first-1]["n"],scan[first]["n"]
    if any(s["cs2"]>=1 for s in scan[:first]):raise AssertionError("first-crossing scan broken")
    root=brentq(lambda n:evaluate(n)["cs2"]-1,left,right,xtol=2e-13,rtol=2e-14)
    root_state=evaluate(root)
    derivative_checks=[]
    for step in (1e-3,3e-4,1e-4,2e-5,1e-5):
        s=evaluate(root,step)
        rr=brentq(lambda n:evaluate(n,step)["cs2_fd"]-1,left,right,xtol=2e-11,rtol=2e-13)
        derivative_checks.append(dict(relative_step=step,root_fd_fm3=rr,
            root_delta_fm3=rr-root,cs2_fd=s["cs2_fd"],cs2_analytic=s["cs2"],
            cs2_delta=s["cs2_fd"]-s["cs2"],firstlaw_relative=s["firstlaw_relative"]))
    # Safety offset is declared before interpolation checks. Its resulting
    # endpoint error is measured in the TOV pass, not hidden in a tolerance.
    safe=root-2e-5
    endpoint=dict(native_causal_root_fm3=root,native_state=root_state,
        native_bracket_fm3=[left,right],native_bracket_chi=[scan[first-1]["cs2"],scan[first]["cs2"]],
        root_analytic_tolerance_fm3=2e-13,fd_convergence=derivative_checks,
        monotonic_chi_on_scan=bool(np.all(np.diff([s["cs2"] for s in scan])>0)),
        admissible_export_endpoint_fm3=safe,cutoff_offset_fm3=root-safe,
        interpretation="first measured native chi=c_s^2/c^2=1 crossing; above it diagnostic only",
        caveat="finite scan plus local root, not a global analytic theorem")
    dump(eosout/"NATIVE_CAUSAL_ENDPOINT.json",endpoint)
    dump(eosout/"native_endpoint_calls.json",dict(binary_sha256=sha(binary),calls=evaluations))
    old=read_compose(SOURCE/"compose")
    oldcontract=json.loads((SOURCE/"REPRESENTATION.json").read_text())
    oldmodel=WholeCrustBarotrope.from_contract(oldcontract)
    anchor=select_rows(read_compose(OLDCORE/"generated"),[-1])
    result=dict(status="CAUSAL_PREFIX_NATIVE_EXTENSION_COMPUTED",causal_endpoint=endpoint,contracts={})
    for label,count in (("coarse",65),("medium",129),("fine",257)):
        dest=eosout/label;dest.mkdir()
        request=dest/"request.nml"
        request.write_text("&generation\n model='RG_SLY4'\n family='core'\n"
            +f" n_points={count}\n n_min=1.0\n n_max={safe:.17e}\n/\n")
        _,receipt=command_log([build/"compose_eos_generate",request,dest/"generated"],dest/"generation")
        core=read_compose(dest/"generated")
        # The native n=1 anchor is regenerated independently as an identity check.
        for k in ("nb","p","eps","mub"):
            a=float(getattr(anchor,k)[0]);b=float(getattr(core,k)[0])
            if abs(a-b)>128*np.finfo(float).eps*max(abs(a),abs(b)):raise ValueError("regenerated native anchor mismatch")
        combined,contract,report=extend_whole_crust(old,oldcontract,core,anchor)
        contract["stellar_extension"]=dict(stage="E",variant=label,
            parent_contract=relative(SOURCE/"REPRESENTATION.json"),parent_sha256=sha(SOURCE/"REPRESENTATION.json"),
            native_root_fm3=root,physical_support_restriction_fm3=safe,
            high_density_extrapolation="NONE",artificial_causal_completion="NONE",energy_reference_shift="NONE")
        write_compose(dest/"compose",combined)
        back=read_compose(dest/"compose")
        readback=serialized_check(combined,back)
        prefix=serialized_check(old,select_rows(back,slice(0,old.n)))
        model=WholeCrustBarotrope.from_contract(contract)
        pp=np.unique(np.r_[0.,old.p,np.geomspace(old.p[0],old.p[-1],401)])
        before=oldmodel.state_at_pressure(pp);after=model.state_at_pressure(pp)
        unchanged={k:bool(np.array_equal(before[k],after[k])) for k in ("n","p","epsilon","mu","enthalpy")}
        if not all(unchanged.values()):raise ValueError("old continuous state changed")
        if not np.array_equal(before["cs2_over_c2"][:-1],after["cs2_over_c2"][:-1]):
            raise ValueError("old interior derivative changed")
        # Pointwise analytic/FD diagnostics are separate from full continuous
        # Hermite interval validation already required by extend_whole_crust.
        sound=np.genfromtxt(dest/"generated/sound_speed.csv",delimiter=",",names=True,dtype=None,encoding=None)
        fd_delta=float(np.max(abs(sound["cs2_over_c2_fd"]-sound["cs2_over_c2_analytic"])))
        if np.any(sound["cs2_over_c2_analytic"]>=1):raise ValueError("export includes native acausal knot")
        provenance=deepcopy(json.loads((SOURCE/"ROW_PROVENANCE.json").read_text()))
        for j in report["native_source_row_indices"]:
            provenance["provenance"].append(dict(origin="native_generated",piece="upper",
                source_row=j,physical_role="native_high_density_extension",
                source_directory=relative(dest/"generated"),source_generation_sha256=sha(dest/"generation.json")))
        provenance["row_count"]=back.n
        provenance["pieces"][-1]["last_row"]=back.n-1
        provenance["energy_reference_shift"]="NONE"
        dump(dest/"ROW_PROVENANCE.json",provenance)
        dump(dest/"REPRESENTATION.json",contract)
        record=dict(**report,generation=receipt,readback=readback,old_prefix=prefix,
            old_continuous_states_exact=unchanged,old_interior_cs2_exact=True,
            old_upper_knot_one_sided_cs2=dict(left=before["cs2_over_c2"][-1],right=after["cs2_over_c2"][-1]),
            max_native_analytic_fd_chi_delta=fd_delta,
            representation_sha256=sha(dest/"REPRESENTATION.json"),
            generated_sha256={p.name:sha(p) for p in sorted((dest/"generated").iterdir()) if p.is_file()})
        dump(dest/"RESULT.json",record)
        result["contracts"][label]=dict(path=relative(dest/"REPRESENTATION.json"),sha256=record["representation_sha256"],
            tail_intervals=count-1,support_n_fm3=list(model.density_domain),continuous_chi_max=report["extension_cs2_max"])
        print("EOS",label,json.dumps(result["contracts"][label]),flush=True)
    dump(eosout/"EOS_RESULT.json",result)
    print("NATIVE_ENDPOINT",json.dumps(endpoint),flush=True)
    return result

def background_pass(out,eosresult):
    tov=out/"tov";tov.mkdir(exist_ok=True)
    if (tov/"TOV_RESULT.json").exists():raise FileExistsError("completed TOV pass exists")
    contracts={k:ContractEnthalpyEOS.from_path(ROOT/v["path"]) for k,v in eosresult["contracts"].items()}
    eos=contracts["fine"]
    endpoint=eosresult["causal_endpoint"]
    safe=endpoint["admissible_export_endpoint_fm3"]
    seq=[]
    for n in np.r_[np.linspace(1.,safe,21),safe]:
        if seq and n==seq[-1]["central_density_fm3"]:continue
        star=solve_background(eos,n,**CONFIG)
        seq.append(star.summary())
        print("TOV",n,star.mass_msun,star.radius_km,flush=True)
        dump(tov/"sequence_checkpoint.json",seq)
    masses=np.array([r["mass_msun"] for r in seq]);densities=np.array([r["central_density_fm3"] for r in seq])
    slopes=np.diff(masses)/np.diff(densities)
    if np.any(slopes<0):
        imax=int(np.argmax(masses))
        if imax==0 or imax==len(seq)-1:raise ValueError("no interior maximum bracket")
        bracket=[densities[imax-1],densities[imax+1]]
        opt=minimize_scalar(lambda n:-solve_background(eos,n,**CONFIG).mass_msun,
            bounds=bracket,method="bounded",options=dict(xatol=2e-7))
        if not opt.success:raise RuntimeError("mass maximum optimizer failed")
        nm=float(opt.x);star=solve_background(eos,nm,**CONFIG)
        classification="CAUSAL_MAXIMUM_MASS_REFINED"
        target_n=nm
        maximum=dict(status=classification,central_density_fm3=nm,bracket_fm3=bracket,summary=star.summary())
    else:
        if not np.all(slopes>0):raise ValueError("mass monotonicity unresolved")
        classification="CAUSAL_BOUNDARY_BEFORE_OBSERVED_MASS_TURNING_POINT"
        target_n=safe
        maximum=dict(status="NOT_DETERMINED_BEFORE_CAUSAL_BOUNDARY",turning_point_claim=False)
    targets=[dict(label="near_endpoint",central_density_fm3=.98*target_n,reason="within common causal support"),
             dict(label="endpoint",central_density_fm3=target_n,reason=classification)]
    selection=dict(primary_contract=eosresult["contracts"]["fine"]["path"],
        comparison_contract=eosresult["contracts"]["coarse"]["path"],models=targets,
        native_causal_density_fm3=endpoint["native_causal_root_fm3"],
        domain_label=classification,endpoint_is_artificial_completion=False)
    dump(tov/"MODE_SELECTION.json",selection)
    print("MODE_SELECTION_READY",json.dumps(selection),flush=True)
    summaries={}
    for target in targets:
        b=solve_background(eos,target["central_density_fm3"],**CONFIG)
        summaries[target["label"]]=b.summary()
        dump(tov/(target["label"].upper()+".json"),b.summary())
        hs=np.unique(np.r_[np.linspace(0.,b.h_start,801),eos.break_h[(eos.break_h>0)&(eos.break_h<b.h_start)]])[::-1]
        rows(tov/(target["label"]+"_profile.csv"),[b.state_at_enthalpy(h) for h in hs])
    ref=summaries["endpoint"];convergence=[]
    variants=[("ode_loose",dict(CONFIG,rtol=1e-8,atol=1e-11,max_step_h=.005)),
        ("reference",CONFIG),("ode_tight",dict(CONFIG,rtol=3e-11,atol=2e-13,max_step_h=.001)),
        ("center_10m",dict(CONFIG,r0_km=.01)),("center_0p3m",dict(CONFIG,r0_km=.0003))]
    for label,cfg in variants:
        b=solve_background(eos,target_n,**cfg)
        convergence.append(dict(label=label,**cfg,mass_msun=b.mass_msun,radius_km=b.radius_km,
            delta_mass_msun=b.mass_msun-ref["mass_msun"],delta_radius_m=1000*(b.radius_km-ref["radius_km"])))
    for label,other in contracts.items():
        for target in targets:
            b=solve_background(other,target["central_density_fm3"],**CONFIG)
            reference=summaries[target["label"]]
            convergence.append(dict(label="grid_"+label+"_"+target["label"],**CONFIG,
                mass_msun=b.mass_msun,radius_km=b.radius_km,delta_mass_msun=b.mass_msun-reference["mass_msun"],
                delta_radius_m=1000*(b.radius_km-reference["radius_km"])))
    rows(tov/"convergence.csv",convergence)
    # One-sided derivative near the actual support endpoint distinguishes a
    # boundary-limited mass from a maximum. Decrease density step explicitly.
    derivatives=[]
    for delta in (.004,.002,.001,.0005):
        b0=solve_background(eos,target_n,**CONFIG)
        b1=solve_background(eos,target_n-delta,**CONFIG)
        b2=solve_background(eos,target_n-2*delta,**CONFIG)
        deriv=(3*b0.mass_msun-4*b1.mass_msun+b2.mass_msun)/(2*delta)
        derivatives.append(dict(step_fm3=delta,dM_dn_Msun_per_fm3=deriv,
            mass_endpoint=b0.mass_msun,mass_minus_h=b1.mass_msun,mass_minus_2h=b2.mass_msun))
    dump(tov/"ENDPOINT_DERIVATIVE.json",derivatives)
    result=dict(status=classification,maximum_mass=maximum,endpoint=ref,
        nearby=summaries["near_endpoint"],native_causal_endpoint=endpoint,
        sequence_points=len(seq),all_sampled_mass_slopes_positive=bool(np.all(slopes>0)),
        last_secant_Msun_per_fm3=float(slopes[-1]),endpoint_derivatives=derivatives,
        convergence=convergence,domain_restriction_error_estimate_Msun=
            abs(derivatives[-1]["dM_dn_Msun_per_fm3"]*(endpoint["native_causal_root_fm3"]-safe)),
        stability_scope="equilibrium sequence and turning-point evidence only; no complete radial spectrum",
        unsupported_acausal_models="NOT_INTEGRATED",metric_signature="(-,+,+,+)")
    dump(tov/"TOV_RESULT.json",result)
    rows(tov/"sequence.csv",[{k:v for k,v in s.items() if k not in ("config","backend")} for s in seq])
    dump(tov/"PROVENANCE.json",dict(base_commit=BASE,runner_sha256=sha(__file__),
        source_contracts=eosresult["contracts"],python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
        platform=platform.platform(),threads={k:os.environ.get(k) for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")},
        invocation=sys.argv,backend="unchanged D exact-contract DOP853 TOV"))
    print("TOV_FINISHED",json.dumps({k:v for k,v in result.items() if k not in ("convergence","native_causal_endpoint")}),flush=True)
    return result

def refine_peak(out,eosresult):
    """Resolve observed adaptive-ODE derivative noise with measured controls.

    Local quadratic interpolation estimates the TOV turning point only; it is
    never used to extend the EOS or replace a native state. Direct symmetric
    derivatives at fixed bracket endpoints certify the sign change.
    """
    tov=out/"tov"
    initial=json.loads((tov/"TOV_RESULT.json").read_text())
    if initial["status"]!="CAUSAL_MAXIMUM_MASS_REFINED":
        raise ValueError("this bounded refinement requires an observed interior maximum")
    refined=tov/"MAXIMUM_REFINEMENT.json"
    if refined.exists():raise FileExistsError("preserve completed maximum refinement")
    contracts={k:ContractEnthalpyEOS.from_path(ROOT/v["path"]) for k,v in eosresult["contracts"].items()}
    control=dict(CONFIG,rtol=1e-12,atol=1e-14,max_step_h=.0005)
    strict=dict(control,max_step_h=.0001)
    center=1.20603
    cache={}
    calls=[]
    def star(label,n,cfg):
        key=(label,float(n),tuple(sorted(cfg.items())))
        if key not in cache:
            b=solve_background(contracts[label],n,**cfg)
            cache[key]=b
            calls.append(dict(grid=label,**b.summary()))
            print("REFINE",label,n,b.mass_msun,b.radius_km,cfg,flush=True)
        return cache[key]
    fits=[]
    for label in ("coarse","medium","fine"):
        for step in (.001,.0005):
            b0=star(label,center,control)
            bm=star(label,center-step,control);bp=star(label,center+step,control)
            derivative=(bp.mass_msun-bm.mass_msun)/(2*step)
            curvature=(bp.mass_msun-2*b0.mass_msun+bm.mass_msun)/step**2
            if curvature>=0:raise ValueError("local turning point is not a maximum")
            nfit=center-derivative/curvature
            mfit=b0.mass_msun-derivative**2/(2*curvature)
            fits.append(dict(grid=label,step_fm3=step,center_n_fm3=center,
                derivative=derivative,curvature=curvature,
                quadratic_turning_density_fm3=nfit,quadratic_maximum_mass_msun=mfit))
    derivative_brackets=[]
    for label in ("coarse","medium","fine"):
        for cfgname,cfg in ([("tight",control),("steplimit",strict)] if label=="fine" else [("tight",control)]):
            vals=[]
            for n in (1.2060,1.2061):
                step=.0005
                bm=star(label,n-step,cfg);bp=star(label,n+step,cfg)
                vals.append((bp.mass_msun-bm.mass_msun)/(2*step))
            rec=dict(grid=label,controls=cfgname,bracket_fm3=[1.2060,1.2061],
                     central_difference_step_fm3=.0005,derivative_at_bracket=vals)
            if not vals[0]>0>vals[1]:raise ValueError("measured turning-point bracket not confirmed")
            derivative_brackets.append(rec)
    final_n=float(fits[-1]["quadratic_turning_density_fm3"])
    final=star("fine",final_n,strict)
    refinements=[]
    for label in ("coarse","medium","fine"):
        b=star(label,final_n,strict)
        refinements.append(dict(label="source_"+label,**b.summary(),
            delta_mass_msun=b.mass_msun-final.mass_msun,delta_radius_m=1000*(b.radius_km-final.radius_km)))
    for label,cfg in (("ode_tight",control),("center_0p3m",dict(strict,r0_km=.0003)),
                      ("center_10m",dict(strict,r0_km=.01))):
        b=star("fine",final_n,cfg)
        refinements.append(dict(label=label,**b.summary(),
            delta_mass_msun=b.mass_msun-final.mass_msun,delta_radius_m=1000*(b.radius_km-final.radius_km)))
    half=(1.2061-1.2060)/2
    mass_bracket_envelope=max(abs(s["curvature"])*half**2/2 for s in fits)
    result=dict(status="EQUILIBRIUM_BAROTROPIC_MAXIMUM_BRACKET_CONVERGED",
        support="all maximum samples below native equilibrium-barotropic causal limit",
        baseline_pass=relative(tov/"TOV_RESULT.json"),initial_optimizer_n_fm3=initial["endpoint"]["central_density_fm3"],
        measured_bracket_fm3=[1.2060,1.2061],quadratic_fits=fits,
        measured_derivative_brackets=derivative_brackets,representative_maximum=final.summary(),
        convergence=refinements,
        mass_bracket_curvature_envelope_msun=mass_bracket_envelope,
        observed_refinement_mass_spread_msun=max(abs(r["delta_mass_msun"]) for r in refinements),
        density_precision_statement="report n_c=1.2060--1.2061 fm^-3; extra digits identify reproducible selected model only",
        initial_derivative_problem="unsegmented adaptive DOP853 gives an isolated approximately 1e-6 Msun excursion at the original control; tighter tolerance and maximum-step refinement restore a consistent sign bracket",
        full_radial_spectrum="NOT_COMPUTED",
        microscopic_causality="NOT_ESTABLISHED_BY_EQUILIBRIUM_BAROTROPE",
        calls=calls)
    dump(refined,result)
    # Resolve the response convention directly against the actual native API.
    binary=out/"eos/native_endpoint_driver"
    native_calls=[]
    def native(n):
        p=subprocess.run([str(binary),format(n,".17g"),"2e-5"],capture_output=True,text=True)
        native_calls.append(dict(n_fm3=n,exit_code=p.returncode,stdout=p.stdout,stderr=p.stderr))
        if p.returncode:raise RuntimeError("native response probe failed")
        values=[float(x) for x in p.stdout.strip().split(",")]
        return dict(n=values[0],epsilon=values[1],p=values[2],mu=values[3],
                    chi_eq=values[4],chi_eq_fd=values[5],chi_fixed_Yp_lepton_equilibrated=values[8])
    scan=np.genfromtxt(out/"eos/native_domain_scan.csv",delimiter=",",names=True)
    j=int(np.flatnonzero(scan["cs2_frozen"]>=1)[0])
    frozen_root=brentq(lambda n:native(n)["chi_fixed_Yp_lepton_equilibrated"]-1,
                      float(scan["n"][j-1]),float(scan["n"][j]),xtol=2e-13)
    scope=dict(equilibrium_barotropic_causal_density_fm3=eosresult["causal_endpoint"]["native_causal_root_fm3"],
        fixed_Yp_lepton_equilibrated_causal_density_fm3=frozen_root,
        at_representative_maximum=native(final_n),
        convention="native frozen output fixes proton fraction but retains common electron/muon chemical potential; it is not fixed individual lepton fractions",
        chi_eq_formula="n*(H_nn-H_ny**2/H_yy)/mu_B",
        chi_fixed_Yp_formula="n*H_nn/mu_B",
        full_fixed_species_response="NOT_COMPUTED",
        maximum_inside_both_examined_response_bounds=bool(final_n<frozen_root<eosresult["causal_endpoint"]["native_causal_root_fm3"]),
        native_calls=native_calls)
    dump(out/"eos/RESPONSE_SCOPE.json",scope)
    print("REFINEMENT_FINISHED",json.dumps({k:v for k,v in result.items() if k not in ("calls","convergence","quadratic_fits")}),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=DEFAULT_OUT)
    parser.add_argument("--build",type=Path,default=DEFAULT_BUILD)
    parser.add_argument("--eos-only",action="store_true")
    parser.add_argument("--refine-only",action="store_true")
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    start=time.time()
    result=make_eos(out,args.build)
    if args.refine_only:refine_peak(out,result)
    elif not args.eos_only:
        background_pass(out,result)
        refine_peak(out,result)
    print("FINISHED elapsed_seconds",time.time()-start,flush=True)
if __name__=="__main__":main()
