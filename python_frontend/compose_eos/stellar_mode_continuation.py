"""Sign-invariant phase-aware eigenfunction continuation for sharp Cowling modes.

Mode labels f/p1 must be backed by eigenfunction shape as the background varies.
Never interpolate a profile across a Maxwell jump or infer identity solely from
an ordered sequence of positive-frequency roots.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import linear_sum_assignment


def _signature(mode,points_per_phase=96):
    profile=mode.get("profile")
    if not isinstance(profile,list) or not profile:
        raise ValueError("nonempty Cowling profile required")
    phases=sorted(set(int(r["phase_index"]) for r in profile))
    if phases != list(range(len(phases))):
        raise ValueError("phase indices must be contiguous from center to surface")
    radii=[]
    components=[]
    for phase in phases:
        subset=[r for r in profile if int(r["phase_index"])==phase]
        if len(subset)<4:
            raise ValueError("too few points in a mode phase")
        data=np.asarray([[float(r["r_km"]),float(r["W"]),
                          float(r["V_km_inverse"])] for r in subset],float)
        if not np.all(np.isfinite(data)) or np.any(data[:,0]<=0):
            raise ValueError("nonfinite or nonphysical eigenfunction profile")
        if np.any(np.diff(data[:,0])<-5e-10):
            raise ValueError("radius must not decrease within a phase")
        radius,indices=np.unique(data[:,0],return_index=True)
        if len(radius)<4 or radius[-1]<=radius[0]:
            raise ValueError("degenerate phase radial support")
        data=data[indices]
        u=(radius-radius[0])/(radius[-1]-radius[0])
        xs=np.linspace(0,1,points_per_phase)
        W=np.interp(xs,u,data[:,1])
        radialV=np.interp(xs,u,radius*data[:,2])
        radii.append([float(radius[0]),float(radius[-1])])
        components.append((W,radialV))
    rtotal=max(x[1] for x in radii)
    # Weight each physical phase by its own radial extent, not its
    # sampled knot count; interpolate *separately* inside each phase.
    length=np.array([x[1]-x[0] for x in radii])
    if any(length<=0) or rtotal<=0:
        raise ValueError("nonphysical phase radius domain")
    weights=np.sqrt(length/rtotal)
    w=np.concatenate([components[k][0]*weights[k] for k in range(len(phases))])
    rv=np.concatenate([components[k][1]*weights[k] for k in range(len(phases))])
    nw,nv=np.linalg.norm(w),np.linalg.norm(rv)
    if min(nw,nv)<1e-14 or not np.all(np.isfinite([nw,nv])):
        raise ValueError("zero or nonfinite Cowling profile norm")
    return np.r_[w/nw/np.sqrt(2),rv/nv/np.sqrt(2)],phases


def eigenfunction_overlap(reference,candidate,*,points_per_phase=96):
    """Phase-separated L2 quadrature proxy with global sign invariance.

    Compare (W,r V) by normalized radius on each phase, independently
    normalizing the two components. This does not compute a canonical
    Sturm-Liouville energy inner product; it is a documented numerical
    shape diagnostic, and cannot certify global QNM mode-family identity.
    """
    if isinstance(points_per_phase,bool) or points_per_phase<12:
        raise ValueError(">=12 collocation points required per phase")
    a,p=_signature(reference,points_per_phase)
    b,q=_signature(candidate,points_per_phase)
    if p!=q:
        raise ValueError("phase layout differs, cannot interpolate across interface")
    dot=float(np.dot(a,b))
    if not np.isfinite(dot):
        raise ValueError("nonfinite eigenfunction overlap")
    return dict(absolute_overlap=float(np.clip(abs(dot),0,1)),
                signed_overlap=dot,n_phases=len(p),
                global_sign_invariant=True,
                inner_product="phase-wise radius-resampled (W,rV) normalized L2 proxy",
                phase_jump_interpolation="NONE")


def track_mode_family(previous,candidates,*,min_overlap=.88,
                     max_fractional_step=.35,min_ambiguity_margin=.03,
                     v_node_policy="strict"):
    """One-to-one adjacent-state assignment with node and overlap gates.

    Frequencies are a *soft* tie breaker, not an identity criterion.
    Strict mode requires both node counts unchanged. The explicitly selected
    phase-aware shape-continuity mode allows one V-node transition only at
    overlap >= .995 and reports the ambiguity as conditional, never as an
    unconditional Sturm node classification.
    """
    if not isinstance(previous,list) or not previous or not isinstance(candidates,list):
        raise ValueError("nonempty previous and candidate mode lists required")
    if v_node_policy not in ("strict","shape_continuity"):
        raise ValueError("v_node_policy must be strict or shape_continuity")
    if (not 0<min_overlap<=1 or not 0<max_fractional_step<2
        or not 0<=min_ambiguity_margin<1):
        raise ValueError("mode-continuation thresholds outside valid ranges")
    if not candidates:
        return dict(status="UNRESOLVED",pairs=[],unresolved_indices=list(range(len(previous))))
    scores=np.full((len(previous),len(candidates)),float("-inf"))
    overlaps=np.zeros_like(scores)
    for i,p in enumerate(previous):
        fp=float(p["frequency_hz"])
        if not np.isfinite(fp) or fp<=0:raise ValueError("invalid reference frequency")
        for j,c in enumerate(candidates):
            fc=float(c["frequency_hz"])
            if not np.isfinite(fc) or fc<=0:continue
            if int(p["W_nodes_per_phase"])!=int(c["W_nodes_per_phase"]):
                continue
            v_change=abs(int(p["V_nodes_per_phase"])-int(c["V_nodes_per_phase"]))
            if v_change and (v_node_policy=="strict" or v_change>1):
                continue
            fraction=abs(fc-fp)/fp
            if fraction>max_fractional_step:continue
            try:
                shape=eigenfunction_overlap(p,c)["absolute_overlap"]
            except ValueError:
                continue
            # A changed V-node count can result from a near-interface zero.
            # Accept *only* the immediately adjacent count with near-identical
            # full W/rV profile. Record the change as conditional, not as a
            # theorem of mode family identity.
            if v_change and shape<.995:
                continue
            overlaps[i,j]=shape
            # Small frequency regularizer changes only ties between similar shapes.
            scores[i,j]=shape-.01*fraction
    cost=np.where(np.isfinite(scores),-scores,1e7)
    assigned_i,assigned_j=linear_sum_assignment(cost)
    accepted=[];unresolved=set(range(len(previous)))
    for i,j in zip(assigned_i,assigned_j):
        top=scores[i,j]
        if not np.isfinite(top) or overlaps[i,j]<min_overlap:continue
        alternatives=[scores[i,k] for k in range(len(candidates))
                      if k!=j and np.isfinite(scores[i,k])]
        if alternatives and top-max(alternatives)<min_ambiguity_margin:continue
        record=dict(previous_index=int(i),current_index=int(j),
                   absolute_overlap=float(overlaps[i,j]),
                   frequency_change_hz=float(candidates[j]["frequency_hz"]-
                                             previous[i]["frequency_hz"]))
        old_v=int(previous[i]["V_nodes_per_phase"])
        new_v=int(candidates[j]["V_nodes_per_phase"])
        if old_v!=new_v:
            record["V_node_transition"]=[old_v,new_v]
            record["conditional_on_shape_overlap_at_least"]=.995
        accepted.append(record)
        unresolved.discard(int(i))
    accepted.sort(key=lambda x:x["previous_index"])
    changed=any("V_node_transition" in p for p in accepted)
    status=("UNRESOLVED" if unresolved else
        "CONTINUED_WITH_V_NODE_CHANGE" if changed else "CONTINUED")
    policy=(
        "conditional phase-separated overlap: W-node unchanged, V-node step <=1 "
        "only for overlap >=0.995; not a Sturm mode-identity theorem"
        if v_node_policy=="shape_continuity" else
        "strict same W/V-node phase-separated eigenfunction overlap; "
        "ambiguous branches remain unresolved")
    return dict(status=status,
                pairs=accepted,unresolved_indices=sorted(unresolved),
                policy=policy)
