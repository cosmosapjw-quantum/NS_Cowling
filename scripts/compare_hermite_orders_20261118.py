"""Two thermodynamic Hermite families for cold P(mu) native EOS.

Cubic: P and dP/dmu=n are matched at frozen training knots.
Quintic: additionally d2P/dmu2=dn/dmu, estimated *only* from the
same training n(mu) via shape-preserving PCHIP. No withheld-data leakage.

No artificial causality clipping and no interpolation across sharp phase
boundaries. Correct energy is recovered from epsilon=mu*n-P.
"""
from __future__ import annotations
from pathlib import Path
import hashlib, json, sys
import numpy as np
from scipy.interpolate import BPoly, CubicHermiteSpline, PchipInterpolator

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from benchmark_native_dense_eos_20261118 import (
    load_native,train_holdout_indices,phase_runs,native_truth_cs2,summarize,
)


def quintic_hermite(mu,p,n,chi):
    mu,p,n,chi=(np.asarray(z,float) for z in (mu,p,n,chi))
    if (any(z.ndim!=1 or z.shape!=mu.shape or np.any(~np.isfinite(z))
            for z in (mu,p,n,chi)) or len(mu)<2 or np.any(np.diff(mu)<=0)
        or np.any(n<=0) or np.any(chi<=0)):
        raise ValueError("positive finite strictly increasing chemical-potential native training support required")
    return BPoly.from_derivatives(mu,
        [[float(p[i]),float(n[i]),float(chi[i])] for i in range(len(mu))])


def derivative_training(mu,n):
    mu=np.asarray(mu,float);n=np.asarray(n,float)
    if np.any(np.diff(mu)<=0) or np.any(np.diff(n)<=0):
        raise ValueError("training native density or chemical potential not increasing")
    return PchipInterpolator(mu,n).derivative()(mu)


def predict(method,source,train,hold):
    if method not in ("cubic","quintic"):
        raise ValueError("only cubic and quintic Hermite methods")
    mu=np.asarray(source["mu"],float)
    p=np.asarray(source["p"],float)
    n=np.asarray(source["n"],float)
    idx=np.asarray(train,int);pick=np.asarray(hold,int)
    if len(idx)<2 or np.any(np.diff(idx)<=0) or np.intersect1d(idx,pick).size:
        raise ValueError("strict independent training and holdout indices required")
    if method=="cubic":
        spline=CubicHermiteSpline(mu[idx],p[idx],n[idx])
    else:
        dndmu=derivative_training(mu[idx],n[idx])
        spline=quintic_hermite(mu[idx],p[idx],n[idx],dndmu)
    ms=mu[pick]
    pp=np.asarray(spline(ms),float)
    nn=np.asarray(spline(ms,nu=1),float)
    chi=np.asarray(spline(ms,nu=2),float)
    eps=ms*nn-pp
    with np.errstate(divide="ignore",invalid="ignore"):
        cs=nn/(ms*chi)
    return dict(mu=ms,p=pp,n=nn,epsilon=eps,cs2=cs,dn_dmu=chi,
                training_second_derivative_policy="PCHIP from sparse training native n(mu) only" if method=="quintic" else "not_used")


def interval_hermite(method,low,high,dn_dmu_left=None,dn_dmu_right=None):
    """A single physics-admitted P(mu) cubic/quintic for TOV and Maxwell segments.

    Does not change physical endpoint states. Causality is checked separately
    against native physics, rather than changed by this constructor.
    """
    mu=np.array([(low["epsilon"]+low["p"])/low["n"],
                 (high["epsilon"]+high["p"])/high["n"]])
    p=np.array([low["p"],high["p"]]);n=np.array([low["n"],high["n"]])
    if method=="cubic":
        return CubicHermiteSpline(mu,p,n)
    if method=="quintic":
        if dn_dmu_left is None or dn_dmu_right is None:
            raise ValueError("quintic needs independent declared one-sided susceptibilities")
        return quintic_hermite(mu,p,n,[dn_dmu_left,dn_dmu_right])
    raise ValueError("unsupported Hermite interpolation method")
