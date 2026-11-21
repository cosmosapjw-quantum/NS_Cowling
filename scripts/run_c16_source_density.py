#!/usr/bin/env python3
"""C16 source-locked local f/p1 node-identity branch continuation.

All physical EOS sound speeds, sharp Maxwell junctions, DVODE
shooting, determinant frequency scan, Brent polishing, and nodal
counts/unique-successor choice execute in Fortran. Python only
unpacks original source EOS, marshals star-density inputs, and
serializes typed records. OS processes isolate DVODE COMMON state.

This is LOCAL node-labelled continuation across ±2% central
density, not a global avoided-crossing or full-GR proof.
"""
from __future__ import annotations
import argparse,ctypes as ct,hashlib,io,json,sys,time,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import stellar_cowling_native as native
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (
   NativeSharpCowlingDVODE,solve_sharp_background_native)
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
MODELS={
 "M2":dict(model="SPG_M2",kind="maxwell",nc=.30385128926615257,
      conversion="rapid",original_p1_reference_hz=5143.719721563306,
      anchor={"f":1842.3116593617328,"p1":5143.71967334719}),
 "M3":dict(model="SPG_M3",kind="css",nc=.192510767556007,
      conversion="slow",original_p1_reference_hz=3736.2088951581254,
      anchor={"f":1564.2448940763072,"p1":3736.2106883570214})
}
FACTORS=(.98,.99,1.,1.01,1.02)
FREQS=np.ascontiguousarray(np.r_[
    np.linspace(700.,3400.,55),
    np.linspace(3450.,5500.,42)],dtype=np.float64)
MAX_BRANCH_STEP_FRACTION=.12

def candidates_for_density(eos,n,config):
    t=time.perf_counter()
    star=solve_sharp_background_native(eos,n,surface_h=3e-8,
       rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    time_bg=time.perf_counter()-t
    opts=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,
       max_step_h=.002,root_xtol_hz=1e-5)
    t=time.perf_counter()
    solver=NativeSharpCowlingDVODE(star,opts,samples=4097,
       coefficient_backend="native_source_cs2",
       source_restart_relative_jump=.01)
    time_coef=time.perf_counter()-t
    t=time.perf_counter()
    brackets,scan=solver.scan_source_brackets_native(FREQS,conversion=config["conversion"])
    time_scan=time.perf_counter()-t
    candidates=[]
    for low,high in brackets:
        item=dict(bracket_hz=[low,high])
        try:
            mode=solver.solve_mode((low,high),conversion=config["conversion"],
                profile_points=360,root_backend="fortran_source_exact")
            item.update(status="PHYSICAL_ROOT_ADMITTED",
                frequency_hz=mode["frequency_hz"],
                source_eos_rhs_callback_calls=mode["source_callback_rhs_calls"],
                phasewise_nodal_backend=mode["node_count_backend"],
                W_nodes=mode["W_nodes_per_phase"],
                V_nodes=mode["V_nodes_per_phase"],
                relative_join_mismatch=mode["relative_join_mismatch"],
                normalized_signed_determinant=mode["matching_residual"],
                original_strong_knot_restarts=mode["root_fortran_restart_count"],
                native_root_iterations=mode["root_iterations"],
                native_root_rhs_evaluations=mode["root_native_rhs_calls"])
            assert item["source_eos_rhs_callback_calls"]==0
            assert item["phasewise_nodal_backend"]=="fortran_phasewise"
            assert item["relative_join_mismatch"]<1e-5
            # Temporary frontend handoff for independent Fortran
            # profile-shape comparison; removed before science JSON write.
            item["_raw_physical_profile"]=mode["profile"]
        except (ValueError,RuntimeError) as exc:
            item.update(status="ROOT_OR_INTERFACE_REJECTED",reason=str(exc))
        candidates.append(item)
    return dict(n_c_fm3=n,mass_msun=star.mass_msun,radius_km=star.radius_km,
        compactness=2*star.mass_geom_km/star.radius_km,
        maxwell_interface_r_km=star.interfaces[0]["r_km"],
        source_restart_count=solver.source_restart_count,
        background_backend=star.config["backend"],
        scan_backend="fortran_source_exact",
        scan_grid_hz=[float(FREQS[0]),float(FREQS[-1]),len(FREQS)],
        scan_point_failures=sum(row["status"]!=0 for row in scan),
        original_scan_samples=scan,
        scan_sign_change_candidate_count=len(brackets),
        candidate_eigenvalues=candidates),dict(background_s=time_bg,
               coefficient_s=time_coef,scan_s=time_scan)

def select_branch_native(candidates,previous_hz,target_nodes):
    if previous_hz is None:
        return dict(status="BLOCKED_PREVIOUS_UNCERTIFIED",frequency_hz=None)
    if not candidates:
        return dict(status="NO_SCAN_BRACKET",frequency_hz=None)
    n=len(candidates)
    if n>128:raise RuntimeError("candidate budget exceeded")
    freq=np.ascontiguousarray([
      c.get("frequency_hz",.5*sum(c["bracket_hz"])) for c in candidates],dtype=np.float64)
    W=np.ascontiguousarray([c.get("W_nodes",-1) for c in candidates],dtype=np.int32)
    V=np.ascontiguousarray([c.get("V_nodes",-1) for c in candidates],dtype=np.int32)
    st=np.ascontiguousarray([
      0 if c["status"]=="PHYSICAL_ROOT_ADMITTED" else 1 for c in candidates],dtype=np.int32)
    idx,status=ct.c_int(-1),ct.c_int(-1)
    with native._LOCK:
        lib=native._native_library()
        lib.cowling_native_select_mode_branch(n,freq,W,V,st,target_nodes,
           float(previous_hz),MAX_BRANCH_STEP_FRACTION,
           ct.byref(idx),ct.byref(status))
    if status.value:
        return dict(status={3:"AMBIGUOUS_ELIGIBLE_MODES",4:"NO_ELIGIBLE_NODE_IDENTIFIED"}
                    .get(status.value,"INVALID_NATIVE_BRANCH_SELECTION"),
                    selector_status=status.value,frequency_hz=None,
                    predecessor_hz=float(previous_hz))
    candidate=candidates[idx.value]
    return dict(status="LOCAL_NODE_CONTINUATION_ADMITTED",selector_status=0,
        candidate_index=idx.value,frequency_hz=candidate["frequency_hz"],
        predecessor_hz=float(previous_hz),
        relative_frequency_step=abs(candidate["frequency_hz"]-previous_hz)/previous_hz,
        W_nodes=candidate["W_nodes"],V_nodes=candidate["V_nodes"],
        relative_join_mismatch=candidate["relative_join_mismatch"])

def compare_mode_shapes_native(previous_profile,current_profile):
    """One-sided W profile correlation in Fortran, not Python.

    Each phase is mapped separately to u in [0,1] so no artificial
    interpolation smears the density jump. This diagnostic does not
    define a conserved or gauge-invariant mode inner product.
    """
    def packed(profile):
        phase=np.ascontiguousarray([r["phase_index"] for r in profile],np.int32)
        radius=np.ascontiguousarray([r["r_km"] for r in profile],np.float64)
        w=np.ascontiguousarray([r["W"] for r in profile],np.float64)
        return phase,radius,w
    phase_a,ra,wa=packed(previous_profile)
    phase_b,rb,wb=packed(current_profile)
    total,inner,outer,status=(ct.c_double(0),ct.c_double(0),
                              ct.c_double(0),ct.c_int(-1))
    with native._LOCK:
        native._native_library().cowling_native_mode_shape_overlap(
          len(ra),phase_a,ra,wa,len(rb),phase_b,rb,wb,
          ct.byref(total),ct.byref(inner),ct.byref(outer),
          ct.byref(status))
    if status.value:
        return dict(status="FORTRAN_SHAPE_OVERLAP_REJECTED",
                    fortran_status=int(status.value))
    return dict(status="FORTRAN_SHAPE_OVERLAP_EVALUATED",
                absolute_radial_W_overlap=float(total.value),
                high_phase_overlap=float(inner.value),
                low_phase_overlap=float(outer.value),
                method="phasewise_normalized_radius_W_only_trapezoid",
                covariant_mode_inner_product=False)


def run(key,out):
    config=MODELS[key]
    raw=SOURCE.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:
        raise RuntimeError("immutable original source EOS archive mismatch")
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
      with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inner:
        contract=inner.read(f"FRESH_5x2_WHOLE_STARS/{config['model']}/{config['kind']}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(contract)))
    results=[]
    for factor in FACTORS:
        nc=config["nc"]*factor
        row=dict(density_factor=factor,n_c_fm3=nc)
        try:
            output,timing=candidates_for_density(eos,nc,config)
            row.update(output)
            # Performance counters deliberately excluded from scientific
            # JSON so the exact same results can be hashed across processes.
            print("C16_WORKER_TIME",key,factor,timing,flush=True)
            print("C16_DENSITY",key,factor,"M",row["mass_msun"],
               "ROOT_CANDIDATES",len(row["candidate_eigenvalues"]),
               "STATUS",[(c["status"],round(c.get("frequency_hz",0),5),c.get("W_nodes")) for c in row["candidate_eigenvalues"]],flush=True)
        except (ValueError,RuntimeError) as exc:
            row.update(status="TOV_OR_GEOMETRY_REJECTED",reason=str(exc),
                candidate_eigenvalues=[])
        results.append(row)
    assert len(results)==5
    byfactor={v["density_factor"]:v for v in results}
    for label,target in (("f",0),("p1",1)):
        # One anchored center plus two one-step density continuations
        # in each direction. The selector will not bridge a missing node.
        center=byfactor[1.0]
        center_choice=select_branch_native(center.get("candidate_eigenvalues",[]),
                            config["anchor"][label],target)
        center.setdefault("local_branches",{})[label]=center_choice
        for direction in ((.99,.98),(1.01,1.02)):
            last=(center_choice["frequency_hz"] if
                      center_choice["status"]=="LOCAL_NODE_CONTINUATION_ADMITTED" else None)
            previous_profile=(center["candidate_eigenvalues"][center_choice["candidate_index"]]
               .get("_raw_physical_profile") if last is not None else None)
            for factor in direction:
                row=byfactor[factor]
                selected=select_branch_native(row.get("candidate_eigenvalues",[]),last,target)
                row.setdefault("local_branches",{})[label]=selected
                if selected["status"]=="LOCAL_NODE_CONTINUATION_ADMITTED":
                    current_profile=row["candidate_eigenvalues"][selected["candidate_index"]].get(
                        "_raw_physical_profile")
                    if previous_profile is not None and current_profile is not None:
                        selected["cross_density_shape_overlap"]=compare_mode_shapes_native(
                            previous_profile,current_profile)
                    previous_profile=current_profile
                    last=selected["frequency_hz"]
                else:
                    previous_profile=None
                    last=None
    # No real mode assigned without phasewise native node identity and the
    # previously locked Maxwell root + profile mismatch gate.
    for row in results:
        f=row["local_branches"]["f"]
        p=row["local_branches"]["p1"]
        if f["status"]=="LOCAL_NODE_CONTINUATION_ADMITTED" and p["status"]=="LOCAL_NODE_CONTINUATION_ADMITTED":
            row["mode_order_valid"]=bool(f["frequency_hz"]<p["frequency_hz"])
            if not row["mode_order_valid"]:
                row["local_branch_classification"]="MODE_ORDER_CONFLICT_UNRESOLVED"
            else:
                row["local_branch_classification"]="NODE_LABELED_MODE_ORDER_PASS"
        else:
            row["local_branch_classification"]="PARTIAL_OR_MISSING_BRANCH"
    # Keep science data and numerical outputs; temporary profile lists
    # are intentionally omitted to reduce published evidence payload.
    for row in results:
        for c in row["candidate_eigenvalues"]:
            c.pop("_raw_physical_profile",None)
    output=dict(schema="ns_cowling.c16.original_source_node_labeled_local_branches.v1",
       parent_commit="8c0f70eed5595b10933763b3917ad4e045b6cf94",
       original_source_archive_sha256=SHA,
       source_contract_sha256=hashlib.sha256(contract).hexdigest(),
       model=config["model"],kind=config["kind"],
       conversion=config["conversion"],nc_anchor_fm3=config["nc"],
       density_factors=list(FACTORS),
       branch_anchors_hz=config["anchor"],
       scan_grid_hz=[float(x) for x in FREQS],
       branch_rule="native Fortran unique node-count plus <12% step, disallow ambiguity",
       cross_density_shape_diagnostic="Fortran phasewise radial W normalized shape overlap; not physical conserved modal inner product",
       physical_scope="cold Gamma1=Gamma_eq, exactly one Maxwell, original contract cs2, no full GR",
       status="SCOPED_LOCAL_CONTINUATION_WITH_REJECTIONS_PRESERVED",
       densities=results)
    out.parent.mkdir(exist_ok=True,parents=True)
    out.write_text(json.dumps(output,sort_keys=True,indent=2,allow_nan=False)+"\n")
    for row in results:
        for label in ("f","p1"):
            record=row["local_branches"][label]
            print("C16_BRANCH",key,row["density_factor"],label,record["status"],
                  record["frequency_hz"],flush=True)
    print("C16_WORKER_DONE",key,out,flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--case",choices=MODELS,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    run(args.case,args.output)
