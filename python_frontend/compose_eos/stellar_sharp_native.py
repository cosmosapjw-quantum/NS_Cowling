"""Native Fortran DVODE TOV/Cowling for equilibrium-barotropic sharp Maxwell stars.

Fortran integrates every TOV and Cowling ODE segment on its own physical EOS
branch. The original Python contract supplies one-sided thermodynamic samples
and phase identity. No interpolation takes place across a density jump.

Slow conversion: continuous radial displacement and Lagrangian pressure.
Rapid conversion: continuous isobar displacement and Lagrangian pressure.
Both share the Tonetto-Lugones convention; this is NOT a general g-mode
or full-GR gravitational radiation solver.
"""
from __future__ import annotations
import ctypes as ct
from dataclasses import dataclass
from math import pi

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from .stellar_background import MEV_FM3_TO_GEOM_KM2
from .stellar_cowling import CowlingSettings, C_KM_S
from . import stellar_cowling_native as native
from .stellar_sharp import SharpContractEnthalpyEOS, SharpStellarBackground, _PiecewiseSolution


def _one_phase_grid(hstart, surface_h, start, end, break_h, samples, *, extras=False):
    """One-sided s samples on descending h interval, including exact endpoints."""
    assert start > end
    slo=np.sqrt(max(hstart-start,0.))
    shi=np.sqrt(max(hstart-end,0.))
    global_ds=np.sqrt(hstart-surface_h)/(samples-1)
    n=max(33,int((shi-slo)/global_ds)+1)
    grid=np.linspace(slo,shi,n)
    cuts=np.asarray(break_h,dtype=float)
    cuts=cuts[(cuts>end)&(cuts<start)]
    if len(cuts):
        marks=np.sqrt(np.maximum(hstart-cuts,0.))
        widths=np.minimum(0.01*global_ds,
             np.maximum(1e-10,2e-12/np.maximum(marks,1e-4)))
        close=np.r_[marks-widths,marks,marks+widths]
        grid=np.r_[grid,close[(close>slo)&(close<shi)]]
    if extras and surface_h>0:
        hs=np.geomspace(surface_h,start,max(1025,n//2))
        additional=np.sqrt(np.maximum(hstart-hs,0.))
        grid=np.r_[grid,additional[(additional>slo)&(additional<shi)]]
    grid=np.unique(np.r_[grid,slo,shi])
    if len(grid)<4:
        raise ValueError("not enough one-sided support samples")
    return grid


def native_shoot_waypoints(targets, source_break_s):
    """Merge physical one-sided EOS breakpoints into a monotone DVODE path.

    Return (path, selected output indices) without leaving the same phase or
    changing the caller's requested output points. No midpoint across a
    Maxwell jump is ever introduced; phase segmentation is upstream.
    """
    wanted=np.asarray(targets,dtype=float)
    cuts=np.asarray(source_break_s,dtype=float)
    if (wanted.ndim!=1 or cuts.ndim!=1 or len(wanted)<2 or
        not np.all(np.isfinite(wanted)) or not np.all(np.isfinite(cuts))):
        raise ValueError("finite 1-D requested targets and source cuts required")
    delta=np.diff(wanted)
    asc=bool(np.all(delta>0))
    desc=bool(np.all(delta<0))
    if not (asc or desc):
        raise ValueError("requested source targets must be strictly monotone")
    lo=float(min(wanted[0],wanted[-1]));hi=float(max(wanted[0],wanted[-1]))
    interior=cuts[(cuts>lo)&(cuts<hi)]
    increasing=np.unique(np.r_[wanted,interior])
    if asc:
        ix=np.searchsorted(increasing,wanted)
        out=increasing
    else:
        ix=len(increasing)-1-np.searchsorted(increasing,wanted)
        out=increasing[::-1]
    if not np.array_equal(out[ix],wanted):
        raise RuntimeError("original shooting output points not preserved")
    return np.ascontiguousarray(out),np.asarray(ix,dtype=int)


class _NativePhase:
    """Lightweight exact-endpoint native phase compatible with _PiecewiseSolution."""
    def __init__(self,hstart, h0,h1, sgrid, profile, nfev):
        self.t=np.asarray([h0,h1],float)
        self.y=np.asarray(profile,dtype=float).T
        self.nfev=int(nfev)
        self.success=True
        self._hstart=hstart
        self._interp=PchipInterpolator(sgrid,np.asarray(profile),axis=0,
                                       extrapolate=False)
    def sol(self,h):
        z=np.asarray(h,dtype=float)
        s=np.sqrt(np.maximum(self._hstart-z,0))
        out=self._interp(s)
        return np.moveaxis(out,-1,0)


def solve_sharp_background_native(eos,central_density_fm3, *,
         r0_km=.001,rtol=1e-9,atol=1e-11,surface_h=0.,
         max_step_h=.005,native_samples=8193,method="adams"):
    if not isinstance(eos,SharpContractEnthalpyEOS):
        raise TypeError("native sharp TOV requires SharpContractEnthalpyEOS")
    if not eos.transitions:
        raise ValueError("sharp native TOV requires at least one validated transition")
    if method not in ("adams","bdf"):
        raise ValueError("method must be adams or bdf")
    if not isinstance(native_samples,int) or native_samples<129:
        raise ValueError("native_samples must be >=129")
    opts=[central_density_fm3,r0_km,rtol,atol,surface_h,max_step_h]
    if not np.all(np.isfinite(opts)) or central_density_fm3<=0 or r0_km<=0 or min(rtol,atol,max_step_h)<=0 or surface_h<0:
        raise ValueError("invalid native TOV configuration")
    hc=float(eos.enthalpy_at_density(central_density_fm3))
    if eos._transition(hc) is not None:
        raise ValueError("central state must lie outside transition bands")
    center=eos.state_at_enthalpy(hc,"high")
    pc=center["p"]*MEV_FM3_TO_GEOM_KM2
    ec=center["epsilon"]*MEV_FM3_TO_GEOM_KM2
    h0=hc-2*pi*(ec/3+pc)*r0_km**2
    if not surface_h<h0<hc:
        raise ValueError("center expansion leaves physical support")
    trans=sorted([t["h"] for t in eos.transitions if surface_h<t["h"]<h0],
                 reverse=True)
    if not trans:
        raise ValueError("central pressure does not admit a resolved sharp interface")
    smax=np.sqrt(h0-surface_h)
    prev=[float(r0_km),4*pi*ec*r0_km**3/3]
    pieces=[];interfaces=[];bounds=[h0,*trans,float(surface_h)]
    mf={"adams":10,"bdf":22}[method]
    with native._LOCK:
        lib=native._native_library()
        for h_begin,h_end in zip(bounds[:-1],bounds[1:]):
            grid=_one_phase_grid(h0,surface_h,h_begin,h_end,eos.break_h,
                                 native_samples,extras=h_end==surface_h)
            midpoint=(h_begin+h_end)/2.
            therm=np.empty((len(grid),2),float)
            for i,s in enumerate(grid):
                h=h0-s*s
                # Preserve the exact one-sided state at phase edges.
                h = h_begin if i==0 else h_end if i==len(grid)-1 else h
                side=eos.side_for_segment(h,midpoint)
                st=eos.state_at_enthalpy(h,side=side)
                therm[i,:]=[st["p"],st["epsilon"]]
            if not np.all(np.isfinite(therm)) or np.any(therm<=0):
                raise ValueError("one-sided sharp EOS contains invalid thermodynamics")
            phys=therm*MEV_FM3_TO_GEOM_KM2
            pchip=PchipInterpolator(grid,phys,axis=0,extrapolate=False)
            slope=pchip(grid,nu=1)
            yout=np.empty((2,len(grid)),order="F",dtype=float)
            nfev,status=ct.c_int(),ct.c_int()
            max_ds=min(smax/200.,max_step_h/(2*smax))
            lib.tov_native_profile(len(grid),np.ascontiguousarray(grid),
                np.ascontiguousarray(phys[:,0]),np.ascontiguousarray(phys[:,1]),
                np.ascontiguousarray(slope[:,0]),np.ascontiguousarray(slope[:,1]),
                float(prev[0]),float(prev[1]),rtol,atol,
                min(1.e-15,atol*1.e-3),max_ds,mf,
                yout,ct.byref(nfev),ct.byref(status))
            if status.value:
                raise RuntimeError(f"native sharp TOV phase [{h_begin},{h_end}] status={status.value}")
            prev=list(map(float,yout[:,-1]))
            phase=_NativePhase(h0,h_begin,h_end,grid,yout.T.copy(),nfev.value)
            pieces.append(phase)
            if h_end in trans:
                low=eos.state_at_enthalpy(h_end,"low")
                high=eos.state_at_enthalpy(h_end,"high")
                interfaces.append(dict(h=float(h_end),r_km=prev[0],
                    m_geom_km=prev[1],r_jump_km=0.,m_jump_km=0.,
                    p_MeV_fm3=float(low["p"]),
                    n_low_fm3=float(low["n"]),n_high_fm3=float(high["n"]),
                    epsilon_low_MeV_fm3=float(low["epsilon"]),
                    epsilon_high_MeV_fm3=float(high["epsilon"])))
    r,m=prev
    if not (r>2*m>0 and r>r0_km):
        raise RuntimeError("native sharp TOV endpoint is nonphysical")
    solution=_PiecewiseSolution(pieces)
    cfg=dict(r0_km=r0_km,rtol=rtol,atol=atol,surface_h=surface_h,
             max_step_h=max_step_h,native_samples=native_samples,
             native_fortran_tov=True,original_python_eos_source=True,method=method,
             no_phase_interpolation=True)
    star=SharpStellarBackground(eos,float(central_density_fm3),hc,h0,
                     float(surface_h),float(r),float(m),solution,cfg,interfaces)
    # Unlike SharpStellarBackground.summary, the native path is explicit.
    star.config["backend"]="Fortran DVODE phasewise TOV"
    return star


def native_interface_transfer(y,source,target,omega,conversion):
    """Native Fortran transfer with independent reconstructed physical residuals."""
    if conversion not in ("slow","rapid"):
        raise ValueError("conversion must be explicitly slow or rapid")
    y=np.ascontiguousarray(np.asarray(y,dtype=float))
    if y.shape!=(2,) or not np.all(np.isfinite(y)):
        raise ValueError("two finite mode components required")
    if not np.isfinite(omega) or omega<=0:
        raise ValueError("positive frequency required")
    for k in ("r_km","Phi","lambda_metric","dPhi_dr_km","p"):
        if not np.isclose(source[k],target[k],rtol=4e-12,atol=1e-14):
            raise ValueError(f"interface must share geometry/pressure: {k}")
    after=np.empty(2,float)
    traction,kinematic,status=ct.c_double(),ct.c_double(),ct.c_int()
    with native._LOCK:
        lib=native._native_library()
        lib.cowling_native_interface_transfer(float(omega),
             float(source["r_km"]),float(source["Phi"]),
             float(source["lambda_metric"]),float(source["dPhi_dr_km"]),
             float(source["p"]),float(source["epsilon"]),float(target["epsilon"]),
             {"slow":1,"rapid":2}[conversion],y,after,
             ct.byref(traction),ct.byref(kinematic),ct.byref(status))
    if status.value:
        raise RuntimeError(f"native Maxwell interface transfer status {status.value}")
    r=source["r_km"]
    a=omega**2*np.exp(source["lambda_metric"]-2*source["Phi"])*r
    b=source["dPhi_dr_km"]
    qi=a*y[1]+b*y[0]
    qo=a*after[1]+b*after[0]
    hi=source["epsilon"]+source["p"]
    ho=target["epsilon"]+target["p"]
    tri=hi*qi;tro=ho*qo
    resid=abs(tri-tro)/max(abs(tri),abs(tro),1e-300)
    mode_index=0 if conversion=="slow" else 1
    kinematic_resid=abs(y[mode_index]-after[mode_index])/max(
        abs(y[mode_index]),abs(after[mode_index]),1e-300)
    diag=dict(traction_relative_residual=float(resid),
        kinematic_relative_residual=float(kinematic_resid),
        native_traction_relative_residual=float(traction.value),
        native_kinematic_relative_residual=float(kinematic.value),
        before=y.tolist(),after=after.tolist(),conversion=conversion,
        h_ratio=float(hi/ho))
    if resid>1e-9 or kinematic_resid>1e-9:
        raise RuntimeError("independent Maxwell junction invariant failed")
    return after,diag


def _nodes(values):
    a=np.asarray(values,float)
    peak=np.max(abs(a))
    if peak==0:
        return 0
    a=a[abs(a)>1e-9*peak]
    return int(np.count_nonzero(a[:-1]*a[1:]<0))


class NativeSharpCowlingDVODE:
    """One-sided native forward shooting through phase interfaces.

    An input profile is constructed only at initialization. Every frequency
    shot uses Fortran DVODE on each phase, then native interface transfer.
    No interpolation or node count spans the density jump.
    """
    def __init__(self,background,settings=CowlingSettings(),*,samples=8193,
                 method="adams",knot_waypoints=False,coefficient_backend="pchip",
                 source_restart_relative_jump=.01):
        if not isinstance(background,SharpStellarBackground) or not background.interfaces:
            raise ValueError("native sharp Cowling needs a sharp stellar background")
        settings.validate(background)
        if settings.surface_h<background.surface_h:
            raise ValueError("Cowling surface outside integrated TOV domain")
        if method not in ("adams","bdf"):
            raise ValueError("native Cowling method must be adams or bdf")
        if samples<129:
            raise ValueError("at least 129 one-sided samples required")
        self.background=background
        self.settings=settings
        self.method=method
        self.samples=samples
        if coefficient_backend not in ("pchip","direct_source","hybrid_cs2","native_source_cs2"):
            raise ValueError("coefficient_backend must be pchip, direct_source, hybrid_cs2 or native_source_cs2")
        if coefficient_backend=="native_source_cs2" and len(background.interfaces)!=1:
            raise ValueError("native_source_cs2 currently requires exactly one Maxwell interface")
        self.coefficient_backend=coefficient_backend
        if source_restart_relative_jump is not None:
            if (not isinstance(source_restart_relative_jump,
                               (int,float,np.integer,np.floating)) or
                isinstance(source_restart_relative_jump,bool) or
                not np.isfinite(source_restart_relative_jump) or
                not 0<source_restart_relative_jump<1):
                raise ValueError("source_restart_relative_jump must lie in (0,1) or None")
        self.source_restart_relative_jump=source_restart_relative_jump
        self._callback_rhs_calls=0
        if not isinstance(knot_waypoints,bool):
            raise ValueError("knot_waypoints must be bool")
        self.knot_waypoints=knot_waypoints
        self.max_interface_residual=0.
        h0=background.h_start; hsurface=settings.surface_h
        bounds=[h0,*sorted([x["h"] for x in background.interfaces
                           if hsurface<x["h"]<h0],reverse=True),hsurface]
        self._phases=[]
        for h_begin,h_end in zip(bounds[:-1],bounds[1:]):
            mid=(h_begin+h_end)/2.
            if self.coefficient_backend=="direct_source":
                # Native Fortran receives background data directly at RHS calls.
                # A minimal source table is needed only for endpoint validation
                # and for the immutable default PCHIP interface contract.
                x=np.linspace(np.sqrt(max(h0-h_begin,0.)),
                              np.sqrt(max(h0-h_end,0.)),33)
            else:
                x=_one_phase_grid(h0,hsurface,h_begin,h_end,background.eos.break_h,
                                  samples,extras=h_end==hsurface)
            values=np.empty((len(x),5),float)
            for i,s in enumerate(x):
                h=h0-s*s
                h=h_begin if i==0 else h_end if i==len(x)-1 else h
                st=background.state_at_enthalpy(h,background.eos.side_for_segment(h,mid))
                cs=st["cs2_over_c2"]
                if not np.isfinite(cs) or cs<=0:
                    raise ValueError("sharp EOS must have positive finite cs²")
                values[i]=[st["r_km"],st["Phi"],st["lambda_metric"],
                           st["dPhi_dr_km"],cs]
            if self.coefficient_backend=="native_source_cs2":
                # The Fortran kernel owns the PCHIP slope construction.
                # No SciPy PCHIP object is created or used in this backend.
                pchip=None
                vals_f=np.asfortranarray(values.T)
                slopes_f=np.empty_like(vals_f,order="F")
                status=ct.c_int(-1)
                with native._LOCK:
                    native._native_library().cowling_native_pchip_slopes(
                        len(x),np.ascontiguousarray(x),vals_f,
                        slopes_f,ct.byref(status))
                if status.value:
                    raise RuntimeError("native one-sided geometry PCHIP rejected: "+
                                       str(status.value))
            else:
                pchip=PchipInterpolator(x,values,axis=0,extrapolate=False)
                slopes_f=np.asfortranarray(pchip(x,nu=1).T)
                vals_f=np.asfortranarray(values.T)
            source_breaks=np.asarray(background.eos.break_h,dtype=float)
            source_breaks=source_breaks[(source_breaks>h_end)&(source_breaks<h_begin)]
            source_break_s=np.sqrt(np.maximum(h0-source_breaks,0.))
            source_restart_s=[]
            if self.coefficient_backend in ("direct_source","hybrid_cs2") and \
                    self.source_restart_relative_jump is not None:
                # Only restart at demonstrated strong EOS sound-speed jumps,
                # not every source table knot (often >1000 per star).
                # Evaluate strictly inside the same physical phase on either
                # side: a Maxwell interface is already a separate phase.
                for bpoint in source_break_s:
                    h=h0-float(bpoint)**2
                    dh=max(2e-11,1e-9*max(h,1e-5))
                    if h-dh<=h_end or h+dh>=h_begin:
                        continue
                    side_lo=background.eos.side_for_segment(h-dh,mid)
                    side_hi=background.eos.side_for_segment(h+dh,mid)
                    lo=background.state_at_enthalpy(h-dh,side_lo)["cs2_over_c2"]
                    hi=background.state_at_enthalpy(h+dh,side_hi)["cs2_over_c2"]
                    if (not np.isfinite(lo) or not np.isfinite(hi) or
                        min(lo,hi)<=0):
                        raise ValueError("source EOS near internal knot is invalid")
                    relative_jump=abs(hi-lo)/max(abs(hi),abs(lo))
                    if relative_jump>=self.source_restart_relative_jump:
                        source_restart_s.append(float(bpoint))
            self._phases.append(dict(h0=h_begin,h1=h_end,mid=mid,
                s=np.ascontiguousarray(x),vals=vals_f,
                slopes=slopes_f,pchip=pchip,
                source_break_s=np.unique(source_break_s),
                source_restart_s=np.unique(np.asarray(source_restart_s,float))))
        self._smax=float(np.sqrt(h0-hsurface))
        self.geometry_interpolation_backend=(
            "fortran_pchip" if self.coefficient_backend=="native_source_cs2"
            else "python_pchip")
        self.source_break_count=sum(len(q["source_break_s"]) for q in self._phases)
        self.source_restart_count=sum(len(q["source_restart_s"]) for q in self._phases)
        self._native_source_eos=None
        if self.coefficient_backend=="native_source_cs2":
            # Frontend packs original source-contract parameters only.
            # All c_s² law selection, inversion and evaluation are Fortran.
            b=background.eos.barotrope
            outer=np.asfortranarray(np.vstack([
                b.outer.nb,b.outer.p,b.outer.F,
                np.r_[b.outer.gamma,b.outer.gamma[-1]]]))
            upper=np.asfortranarray(np.vstack([
                b.upper_n,b.upper_p,b.upper_epsilon]))
            surface=np.ascontiguousarray([
                b.surface.rest_energy_MeV,b.surface.n0,
                b.surface.p0,b.surface.gamma],dtype=np.float64)
            bridge=np.ascontiguousarray([
                b.bridge.n_left,b.bridge.p_left,b.bridge.epsilon_left,
                b.bridge.n_right,b.bridge.p_right,b.bridge.epsilon_right,
                b.bridge.kappa],dtype=np.float64)
            jump=int(next(iter(b.sharp_interfaces)))
            self._native_source_eos=(outer,upper,surface,bridge,jump)
            if self.source_restart_relative_jump is not None:
                # Strong source derivative jumps are classified by the
                # original EOS evaluator inside Fortran, not in Python.
                with native._LOCK:
                    dll=native._native_library()
                    self._native_exact_eos_setup(dll)
                    for phase in self._phases:
                        capacity=outer.shape[1]+upper.shape[1]+5
                        out=np.zeros(capacity,dtype=np.float64)
                        count,status=ct.c_int(0),ct.c_int(-1)
                        dll.cowling_native_source_select_restarts(
                            capacity,h0,phase["h0"],phase["h1"],
                            self.source_restart_relative_jump,
                            2 if phase is self._phases[0] else 1,
                            out,ct.byref(count),ct.byref(status))
                        if status.value:
                            raise RuntimeError("Fortran strong EOS knot selection failed: "+
                                               str(status.value))
                        phase["source_restart_s"]=np.ascontiguousarray(
                            out[:count.value],dtype=np.float64)
                self.source_restart_count=sum(
                    len(q["source_restart_s"]) for q in self._phases)

    def _native_phase_profile_values(self,phase,coordinates):
        """Evaluate phase-local PCHIP fields in Fortran, no SciPy path."""
        if self.coefficient_backend!="native_source_cs2":
            raise ValueError("native profile interpolation requires native_source_cs2")
        x=np.ascontiguousarray(coordinates,dtype=np.float64)
        out=np.empty((5,len(x)),dtype=np.float64,order="F")
        st=ct.c_int(-1)
        with native._LOCK:
            native._native_library().cowling_native_pchip_values(
                len(phase["s"]),phase["s"],phase["vals"],phase["slopes"],
                len(x),x,out,ct.byref(st))
        if st.value:
            raise RuntimeError("native PCHIP profile query failed: "+str(st.value))
        return out

    def _native_exact_eos_setup(self,lib):
        if self._native_source_eos is None:
            raise RuntimeError("native source EOS contract not configured")
        outer,upper,surface,bridge,jump=self._native_source_eos
        status=ct.c_int(-1)
        lib.cowling_native_source_eos_setup(
            outer.shape[1],outer,upper.shape[1],upper,surface,bridge,
            jump,ct.byref(status))
        if status.value:
            raise RuntimeError("Fortran original EOS law setup rejected: "+str(status.value))

    def _source_coefficients(self,phase,svalue):
        """One original-EOS c_s² with phase-local PCHIP geometry only.

        Never interpolate thermodynamic derivatives across EOS derivative
        knots; the native integrator restarts at sufficiently strong jumps.
        Original physical Maxwell interfaces remain separate phases.
        """
        if self.coefficient_backend!="hybrid_cs2":
            raise ValueError("source coefficients require hybrid_cs2 backend")
        s=float(np.clip(svalue,phase["s"][0],phase["s"][-1]))
        h=float(np.clip(self.background.h_start-s*s,
                        phase["h1"],phase["h0"]))
        geom=phase["pchip"](s)
        cs=self.background.eos.state_at_enthalpy(
            h,self.background.eos.side_for_segment(h,phase["mid"]))["cs2_over_c2"]
        return (float(geom[0]),float(geom[1]),float(geom[2]),
                float(geom[3]),float(cs))

    def _integrate_phase(self, phase, targets, y, omega, lib):
        x=phase["s"]
        status=ct.c_int()
        lib.cowling_native_setup(len(x),x,phase["vals"],phase["slopes"],
                                  ct.byref(status))
        if status.value:
            raise RuntimeError("native one-sided phase setup rejected: "+str(status.value))
        native._ACTIVE=None
        targets=np.ascontiguousarray(targets,dtype=float)
        # Physical source discontinuity stops are stronger than ordinary
        # output waypoints: DVODE is restarted at each selected cs² jump.
        if self.coefficient_backend in ("direct_source","hybrid_cs2","native_source_cs2"):
            cut_s=(phase["source_break_s"] if self.knot_waypoints
                   else phase["source_restart_s"])
        elif self.knot_waypoints:
            cut_s=phase["source_break_s"]
        else:
            cut_s=np.asarray([],float)
        if len(cut_s):
            all_targets,wanted_indices=native_shoot_waypoints(targets,cut_s)
        else:
            all_targets=targets;wanted_indices=np.arange(len(targets))
        restart_indices=np.flatnonzero(np.isin(all_targets,phase["source_restart_s"]))
        restart_indices=restart_indices[(restart_indices>0)&
                                         (restart_indices<len(all_targets)-1)]
        restart_indices=np.ascontiguousarray(restart_indices+1,dtype=np.int32)
        out=np.empty((2,len(all_targets)),dtype=float,order="F")
        calls=ct.c_int()
        max_ds=min(self._smax/200,
                   self.settings.max_step_h/(2*self._smax))
        callback=None
        native_source_registered=False
        cb_errors=[]
        cb_status=ct.c_int()
        if self.coefficient_backend=="native_source_cs2":
            source_side=2 if phase is self._phases[0] else 1
            lib.cowling_native_source_phase_setup(
                len(x),x,phase["vals"],phase["slopes"],
                self.background.h_start,source_side,ct.byref(cb_status))
            if cb_status.value:
                raise RuntimeError("Fortran source geometry setup rejected: "+
                                   str(cb_status.value))
            lib.cowling_native_source_register(ct.byref(cb_status))
            if cb_status.value:
                raise RuntimeError("Fortran source EOS callback registration rejected: "+
                                   str(cb_status.value))
            native_source_registered=True
        if self.coefficient_backend in ("direct_source","hybrid_cs2"):
            # The function pointer is held strongly until Fortran returns.
            # Fortran and source EOS callbacks are process-local and always
            # serialized by native._LOCK in the calling matched-shot path.
            CB=ct.CFUNCTYPE(None,ct.c_double,ct.POINTER(ct.c_double),
                            ct.POINTER(ct.c_int))
            def exact_background(svalue,coeff_ptr,status_ptr):
                status_ptr[0]=1
                try:
                    svalue=float(np.clip(svalue,x[0],x[-1]))
                    h=self.background.h_start-svalue*svalue
                    h=float(np.clip(h,phase["h1"],phase["h0"]))
                    if self.coefficient_backend=="hybrid_cs2":
                        values=self._source_coefficients(phase,svalue)
                    else:
                        st=self.background.state_at_enthalpy(
                            h,self.background.eos.side_for_segment(h,phase["mid"]))
                        values=(st["r_km"],st["Phi"],st["lambda_metric"],
                                st["dPhi_dr_km"],st["cs2_over_c2"])
                    if (not np.all(np.isfinite(values)) or
                         values[0]<=0 or values[3]<=0 or values[4]<=0):
                        raise ValueError("source background has nonphysical coefficients")
                    for i,v in enumerate(values):coeff_ptr[i]=float(v)
                    self._callback_rhs_calls+=1
                    status_ptr[0]=0
                except Exception as error:
                    cb_errors.append(f"{type(error).__name__}: {error}")
            callback=CB(exact_background)
            lib.cowling_native_set_source_callback(
                  ct.cast(callback,ct.c_void_p),ct.byref(cb_status))
            if cb_status.value:
                raise RuntimeError("native EOS callback registration failed")
        try:
            if self.coefficient_backend in ("direct_source","hybrid_cs2","native_source_cs2") and len(restart_indices):
                lib.cowling_native_segment_restarts(
                    omega,self.settings.ell,self.settings.rtol,
                    self.settings.atol,max_ds,
                    {"adams":10,"bdf":22}[self.method],
                    len(all_targets),all_targets,
                    len(restart_indices),restart_indices,
                    np.ascontiguousarray(y),out,ct.byref(calls),ct.byref(status))
            else:
                lib.cowling_native_segment(omega,self.settings.ell,
                    self.settings.rtol,self.settings.atol,max_ds,
                    {"adams":10,"bdf":22}[self.method],
                    len(all_targets),all_targets,np.ascontiguousarray(y),
                    out,ct.byref(calls),ct.byref(status))
        finally:
            if callback is not None:
                lib.cowling_native_set_source_callback(None,ct.byref(cb_status))
                if cb_status.value:
                    raise RuntimeError("native EOS callback reset failed")
            if native_source_registered:
                lib.cowling_native_source_unregister(ct.byref(cb_status))
                if cb_status.value:
                    raise RuntimeError("Fortran original EOS callback reset failed")
        if cb_errors:
            raise RuntimeError("native direct source callback failed: "+cb_errors[0])
        if status.value:
            raise RuntimeError("native one-sided phase DVODE rejected: "
                               +str(status.value))
        return np.asfortranarray(out[:,wanted_indices]),int(calls.value)

    def _matched_shot(self,frequency_hz,conversion,profile_points=0):
        if conversion not in ("slow","rapid"):
            raise ValueError("conversion must be slow or rapid")
        if not np.isfinite(frequency_hz) or frequency_hz<=0:
            raise ValueError("finite positive frequency_hz required")
        if not self._phases:
            raise ValueError("empty native sharp background")
        omega=2*pi*float(frequency_hz)/C_KM_S
        # Keep match strictly inside the high-density phase. The mode
        # eigenvalue is independent of this arbitrary numerical separation.
        smatch=0.5*self._phases[0]["s"][-1]
        count=max(2,profile_points//(len(self._phases)+1))
        traces=[];junctions=[];nfev=0
        with native._LOCK:
            lib=native._native_library()
            if self.coefficient_backend=="native_source_cs2":
                self._native_exact_eos_setup(lib)
            phase=self._phases[0]
            x_inner=np.linspace(phase["s"][0],smatch,count)
            inner,calls=self._integrate_phase(phase,x_inner,[1.,-1./self.settings.ell],omega,lib)
            nfev+=calls
            ui=inner[:,-1].copy()
            # The surface starting state imposes DeltaP=0 exactly.
            sdata=self.background.state_at_enthalpy(self.settings.surface_h)
            asurf=omega**2*np.exp(sdata["lambda_metric"]-2*sdata["Phi"])*sdata["r_km"]
            bsurf=sdata["dPhi_dr_km"]
            if asurf<=0 or bsurf<=0:
                raise RuntimeError("nonpositive native surface coefficients")
            y=np.ascontiguousarray([1.,-bsurf/asurf],dtype=float)
            surface_bc=(asurf*y[1]+bsurf*y[0])/(
                np.hypot(asurf,bsurf)*np.linalg.norm(y))
            # Outward -> inward: integrate each one-sided physical phase
            # BACKWARD from the surface, crossing a junction exactly.
            for j in range(len(self._phases)-1,-1,-1):
                phase=self._phases[j]
                start=phase["s"][-1]
                stop=smatch if j==0 else phase["s"][0]
                xs_back=np.linspace(start,stop,count)
                yy,calls=self._integrate_phase(phase,xs_back,y,omega,lib)
                nfev+=calls
                traces.append((j,xs_back,yy.copy(),phase["pchip"]))
                y=np.ascontiguousarray(yy[:,-1])
                if j>0:
                    ht=phase["h0"]
                    low=self.background.state_at_enthalpy(ht,"low")
                    high=self.background.state_at_enthalpy(ht,"high")
                    y,diag=native_interface_transfer(y,low,high,omega,conversion)
                    diag["h"]=float(ht)
                    diag["r_km"]=float(low["r_km"])
                    diag["direction"]="outer_to_inner"
                    junctions.append(diag)
                    self.max_interface_residual=max(
                        self.max_interface_residual,
                        diag["traction_relative_residual"],
                        diag["kinematic_relative_residual"])
            uo=y.copy()
        den=np.linalg.norm(ui)*np.linalg.norm(uo)
        if den==0 or not np.isfinite(den):
            raise RuntimeError("singular matched native Cowling state")
        determinant=(ui[0]*uo[1]-ui[1]*uo[0])/den
        scale=float(np.dot(ui,uo)/np.dot(uo,uo))
        reljoin=float(np.linalg.norm(ui-scale*uo)/np.linalg.norm(ui))
        return dict(residual=float(determinant),ui=ui,uo=uo,
            scale=scale,join_relative=reljoin,inside=(x_inner,inner,self._phases[0]["pchip"]),
            backwards=traces,junctions=junctions,nfev=nfev,
            surface_bc_relative=abs(float(surface_bc)))

    def matching_residual(self,frequency_hz,*,conversion):
        return self._matched_shot(frequency_hz,conversion,0)["residual"]

    def scan_brackets(self,frequencies_hz,*,conversion):
        fs=np.asarray(frequencies_hz,float)
        if fs.ndim!=1 or len(fs)<2 or np.any(~np.isfinite(fs)) or np.any(np.diff(fs)<=0) or fs[0]<=0:
            raise ValueError("strict positive monotone frequency scan required")
        vals=[self.matching_residual(f,conversion=conversion) for f in fs]
        return ([(float(fs[i]),float(fs[i+1])) for i in range(len(fs)-1)
                 if vals[i]*vals[i+1]<0],
                [dict(frequency_hz=float(f),residual=float(v))
                 for f,v in zip(fs,vals)])

    def _native_root_phase_pchip(self,bracket_hz,*,conversion):
        """Fortran-only numerical root of phasewise PCHIP Maxwell determinant.

        A deliberately opt-in, source-approximate backend: the original
        hybrid_cs2/direct_source sound-speed callbacks must never be silently
        substituted by interpolated PCHIP sound speed. Exactly one Maxwell
        interface is admitted, with the original slow/rapid traction transfer.
        """
        if conversion not in ("slow","rapid"):
            raise ValueError("conversion must be slow or rapid")
        if self.coefficient_backend!="pchip":
            raise ValueError("Fortran phase root requires pchip; original source cs2 needs its own backend")
        if len(self._phases)!=2 or len(self.background.interfaces)!=1:
            raise ValueError("Fortran phase root requires exactly one Maxwell interface")
        if self.knot_waypoints:
            raise ValueError("Fortran phase root requires knot_waypoints=False for matched path")
        lo,hi=map(float,bracket_hz)
        if not np.isfinite(lo) or not np.isfinite(hi) or not 0<lo<hi:
            raise ValueError("positive ordered native root bracket required")
        high,low=self._phases
        ht=self.background.interfaces[0]["h"]
        state_low=self.background.state_at_enthalpy(ht,"low")
        state_high=self.background.state_at_enthalpy(ht,"high")
        for field in ("r_km","Phi","lambda_metric","dPhi_dr_km","p"):
            if not np.isclose(state_low[field],state_high[field],rtol=4e-12,atol=1e-14):
                raise ValueError("Maxwell interface source geometry is discontinuous: "+field)
        interface=np.ascontiguousarray([
            state_low["r_km"],state_low["Phi"],state_low["lambda_metric"],
            state_low["dPhi_dr_km"],state_low["p"],
            state_low["epsilon"],state_high["epsilon"]],dtype=np.float64)
        smatch=0.5*float(high["s"][-1])
        max_ds=min(self._smax/200,self.settings.max_step_h/(2*self._smax))
        frequency,residual=ct.c_double(0),ct.c_double(0)
        steps,calls,status=ct.c_int(0),ct.c_int(0),ct.c_int(-1)
        with native._LOCK:
            lib=native._native_library()
            try:
                lib.cowling_native_sharp_onejump_brent(
                    lo,hi,self.settings.ell,self.settings.rtol,self.settings.atol,
                    max_ds,{"adams":10,"bdf":22}[self.method],
                    len(high["s"]),high["s"],high["vals"],high["slopes"],
                    len(low["s"]),low["s"],low["vals"],low["slopes"],
                    smatch,interface,{"slow":1,"rapid":2}[conversion],
                    self.settings.root_xtol_hz,ct.byref(frequency),
                    ct.byref(residual),ct.byref(steps),ct.byref(calls),
                    ct.byref(status))
            finally:
                # The Fortran phase-root implementation replaces the module's
                # active interpolant on every physical phase. Invalidate the
                # smooth-star frontend's cached setup identity.
                native._ACTIVE=None
        if status.value==3:
            raise ValueError("native sharp phase root bracket lacks a sign change")
        if status.value==33:
            raise RuntimeError("native sharp phase root sign flip not physically admitted")
        if status.value:
            raise RuntimeError("Fortran sharp phase root rejected: "+str(status.value))
        if abs(residual.value)>1e-5 or not np.isfinite(frequency.value):
            raise RuntimeError("Fortran sharp root determinant outside acceptance threshold")
        return float(frequency.value),int(steps.value),int(calls.value)

    def _native_source_exact_root(self,bracket_hz,*,conversion):
        """Thin frontend for original-source Fortran Brent and knot restarts.

        All determinant shots, analytic c_s² evaluations, selected derivative
        jumps and DVODE history restarts are native. Existing SciPy Brent
        stays available as an explicit scientific reference.
        """
        if conversion not in ("slow","rapid"):
            raise ValueError("conversion must be slow or rapid")
        if self.coefficient_backend!="native_source_cs2":
            raise ValueError("Fortran source-exact root requires native_source_cs2 backend")
        if len(self._phases)!=2 or len(self.background.interfaces)!=1:
            raise ValueError("Fortran native sharp root requires one Maxwell interface")
        if self.knot_waypoints:
            raise ValueError("Fortran root requires knot_waypoints=False")
        if self.source_restart_relative_jump is None:
            raise ValueError("Fortran root requires explicit restart jump threshold")
        lo,hi=map(float,bracket_hz)
        if not (np.isfinite(lo) and np.isfinite(hi) and 0<lo<hi):
            raise ValueError("ordered positive frequency bracket required")
        high,low=self._phases
        ht=self.background.interfaces[0]["h"]
        state_low=self.background.state_at_enthalpy(ht,"low")
        state_high=self.background.state_at_enthalpy(ht,"high")
        for key in ("r_km","Phi","lambda_metric","dPhi_dr_km","p"):
            if not np.isclose(state_low[key],state_high[key],
                              rtol=4e-12,atol=1e-14):
                raise ValueError("Maxwell interface source geometry mismatch: "+key)
        iface=np.ascontiguousarray([
           state_low["r_km"],state_low["Phi"],state_low["lambda_metric"],
           state_low["dPhi_dr_km"],state_low["p"],
           state_low["epsilon"],state_high["epsilon"]],dtype=np.float64)
        root,det=ct.c_double(0),ct.c_double(0)
        iterations,calls,cuts,status=(ct.c_int(0),ct.c_int(0),
                                    ct.c_int(0),ct.c_int(-1))
        max_ds=min(self._smax/200,
                   self.settings.max_step_h/(2*self._smax))
        with native._LOCK:
            lib=native._native_library()
            self._native_exact_eos_setup(lib)
            try:
                lib.cowling_native_sharp_source_brent(
                    lo,hi,self.settings.ell,self.settings.rtol,self.settings.atol,
                    max_ds,{"adams":10,"bdf":22}[self.method],
                    len(high["s"]),high["s"],high["vals"],high["slopes"],
                    len(low["s"]),low["s"],low["vals"],low["slopes"],
                    self.background.h_start,
                    0.5*float(high["s"][-1]),iface,
                    {"slow":1,"rapid":2}[conversion],
                    self.settings.root_xtol_hz,self.source_restart_relative_jump,
                    ct.byref(root),ct.byref(det),ct.byref(iterations),
                    ct.byref(calls),ct.byref(cuts),ct.byref(status))
            finally:
                native._ACTIVE=None
        if status.value==3:
            raise ValueError("Fortran source-exact root bracket lacks sign change")
        if status.value==33:
            raise RuntimeError("Fortran source-exact root bracket discontinuous or underresolved")
        if status.value:
            raise RuntimeError("Fortran source-exact root failed status="+str(status.value))
        if not np.isfinite(root.value) or abs(det.value)>1e-5:
            raise RuntimeError("Fortran source-exact root outside determinant gate")
        return float(root.value),int(iterations.value),int(calls.value),int(cuts.value)

    def sweep_native_numerical_controls(self,bracket_hz,*,conversion,controls):
        """Frontend-only marshaling for a Fortran source-locked error budget.

        Each control is a mapping with rtol, max_step_s and method
        ('adams'|'bdf'). This routine calls ONE Fortran batch; neither
        the original EOS c_s² nor any root residual is evaluated here.
        A zero/failed numeric status is preserved, never promoted into
        an admissible mode or fabricated eigenfrequency.

        NOTE: Roots are diagnostics. Scientific mode/node admissibility
        still requires the independent eigenfunction reconstruction.
        """
        if self.coefficient_backend!="native_source_cs2":
            raise ValueError("sweep requires the native_source_cs2 backend")
        if conversion not in ("slow","rapid"):
            raise ValueError("conversion must be slow or rapid")
        if len(self._phases)!=2 or len(self.background.interfaces)!=1:
            raise ValueError("sweep requires exactly one Maxwell interface")
        if self.knot_waypoints or self.source_restart_relative_jump is None:
            raise ValueError("sweep requires native strong-knot configuration")
        a,b=map(float,bracket_hz)
        if not (np.isfinite(a) and np.isfinite(b) and 0<a<b):
            raise ValueError("ordered positive sweep bracket required")
        controls=list(controls)
        if not controls or len(controls)>10000:
            raise ValueError("sweep controls must have 1..10000 entries")
        method_map={"adams":10,"bdf":22}
        for item in controls:
            if item["method"] not in method_map:
                raise ValueError("unsupported DVODE method in controls")
        tol=np.ascontiguousarray([float(c["rtol"]) for c in controls],dtype=float)
        step=np.ascontiguousarray([float(c["max_step_s"]) for c in controls],dtype=float)
        mf=np.ascontiguousarray([method_map[c["method"]] for c in controls],dtype=np.int32)
        count=len(controls)
        frequency=np.zeros(count,dtype=float)
        det=np.zeros(count,dtype=float)
        iterations=np.zeros(count,dtype=np.int32)
        calls=np.zeros(count,dtype=np.int32)
        cuts=np.zeros(count,dtype=np.int32)
        status=np.full(count,-1,dtype=np.int32)
        high,low=self._phases
        ht=self.background.interfaces[0]["h"]
        stlow=self.background.state_at_enthalpy(ht,"low")
        sthigh=self.background.state_at_enthalpy(ht,"high")
        for key in ("r_km","Phi","lambda_metric","dPhi_dr_km","p"):
            if not np.isclose(stlow[key],sthigh[key],rtol=4e-12,atol=1e-14):
                raise ValueError("Maxwell interface geometry mismatch in sweep: "+key)
        geom=np.ascontiguousarray([
            stlow["r_km"],stlow["Phi"],stlow["lambda_metric"],
            stlow["dPhi_dr_km"],stlow["p"],stlow["epsilon"],
            sthigh["epsilon"]],dtype=np.float64)
        with native._LOCK:
            lib=native._native_library()
            self._native_exact_eos_setup(lib)
            try:
                lib.cowling_native_sharp_source_sweep(
                    count,a,b,self.settings.ell,self.settings.atol,
                    self.background.h_start,0.5*float(high["s"][-1]),
                    geom,{"slow":1,"rapid":2}[conversion],
                    self.settings.root_xtol_hz,self.source_restart_relative_jump,
                    len(high["s"]),high["s"],high["vals"],high["slopes"],
                    len(low["s"]),low["s"],low["vals"],low["slopes"],
                    tol,step,mf,frequency,det,iterations,calls,cuts,status)
            finally:
                native._ACTIVE=None
        results=[]
        for i,c in enumerate(controls):
            ok=int(status[i])==0 and abs(det[i])<=1e-5
            code=int(status[i]) if ok or status[i]!=0 else 91
            results.append(dict(
                rtol=float(tol[i]),max_step_s=float(step[i]),
                method=c["method"],status=code,
                frequency_hz=float(frequency[i]) if ok else None,
                matching_residual=float(det[i]) if ok else None,
                root_iterations=int(iterations[i]) if ok else 0,
                native_rhs_evaluations=int(calls[i]) if ok else 0,
                source_restart_count=int(cuts[i]) if ok else 0))
        return results

    def solve_mode(self,bracket_hz,*,conversion,profile_points=600,
                   root_backend="python_brentq"):
        lo,hi=map(float,bracket_hz)
        if not 0<lo<hi or profile_points<40:
            raise ValueError("positive bracket and >=40 profile points required")
        root_ncuts=None
        if root_backend=="python_brentq":
            root,info=brentq(lambda f:self.matching_residual(f,conversion=conversion),
                            lo,hi,xtol=self.settings.root_xtol_hz,
                            rtol=4*np.finfo(float).eps,full_output=True)
            root_converged=bool(info.converged)
            root_iterations=int(info.iterations)
            root_nfev=None
        elif root_backend=="fortran_phase_pchip":
            root,root_iterations,root_nfev=self._native_root_phase_pchip(
                bracket_hz,conversion=conversion)
            root_converged=True
        elif root_backend=="fortran_source_exact":
            root,root_iterations,root_nfev,root_ncuts=self._native_source_exact_root(
                bracket_hz,conversion=conversion)
            root_converged=True
        else:
            raise ValueError("unsupported sharp Cowling root_backend")
        shot=self._matched_shot(root,conversion,profile_points)
        residual=float(shot["residual"])
        if abs(residual)>1e-5 or shot["join_relative"]>1e-5:
            raise RuntimeError("mode root not physically admitted: "+
                f"matching={residual} join={shot['join_relative']}")
        ell=self.settings.ell
        x0,y0,i0=shot["inside"]
        segments=[(0,x0,y0,i0)]
        # Reverse every backpropagated segment to normal center->surface
        for j,xs,ys,interp in reversed(shot["backwards"]):
            segments.append((j,xs[::-1],(ys*shot["scale"])[:,::-1],interp))
        # Merge segments sharing the same physical phase, not across a jump.
        byphase={}
        for j,xs,ys,interp in segments:
            byphase.setdefault(j,[]).append((xs,ys,interp))
        profile=[];wcount=vcount=0;norm=0.
        for j in sorted(byphase):
            data=byphase[j]
            xs=np.concatenate([x for x,y,interp in data])
            ys=np.hstack([y for x,y,interp in data])
            interp=data[0][2]
            if self.coefficient_backend=="direct_source":
                # Profile radii are source evaluated too, not from the
                # deliberately sparse PCHIP sentinel table.
                phase=self._phases[j]
                r=np.asarray([self.background.state_at_enthalpy(
                    float(np.clip(self.background.h_start-s*s,
                                  phase["h1"],phase["h0"])),
                    self.background.eos.side_for_segment(
                        float(np.clip(self.background.h_start-s*s,
                                      phase["h1"],phase["h0"])),
                        phase["mid"]))["r_km"] for s in xs],float)
            elif self.coefficient_backend=="native_source_cs2":
                r=self._native_phase_profile_values(self._phases[j],xs)[0]
            else:
                r=interp(xs)[:,0]
            W=r**(ell+1)*ys[0]
            V=r**ell*ys[1]
            wcount+=_nodes(W);vcount+=_nodes(V)
            norm=max(norm,float(np.max(abs(W))))
            for s,radius,w,v in zip(xs,r,W,V):
                profile.append(dict(h=float(self.background.h_start-s*s),
                    r_km=float(radius),W=float(w),V_km_inverse=float(v),
                    phase_index=j))
        if norm<=0 or not np.isfinite(norm):
            raise RuntimeError("invalid native matched eigenfunction amplitude")
        for row in profile:
            row["W"]/=norm;row["V_km_inverse"]/=norm
        return dict(frequency_hz=float(root),bracket_hz=[lo,hi],
                    conversion=conversion,ell=ell,
                    perturbation_closure="GAMMA1_EQUALS_GAMMA_EQ",
                    frozen_composition_gamma1="NOT_DERIVABLE_FROM_EQUILIBRIUM_BAROTROPE",
                    bulk_composition_g_modes_supported=False,
                    root_converged=root_converged,root_iterations=root_iterations,
                    root_backend=root_backend,
                    root_native_rhs_calls=root_nfev,
                    root_fortran_restart_count=root_ncuts,
                    matching_residual=float(residual),
                    relative_join_mismatch=shot["join_relative"],
                    surface_BC_relative=shot["surface_bc_relative"],
                    interface_junctions=shot["junctions"],
                    W_nodes_per_phase=wcount,V_nodes_per_phase=vcount,
                    profile=profile,nfev=shot["nfev"],
                    backend="Fortran DVODE sharp "+self.method.capitalize(),
                    coefficient_backend=self.coefficient_backend,
                    geometry_interpolation_backend=self.geometry_interpolation_backend,
                    strong_knot_selector_backend=("fortran_source_exact"
                        if self.coefficient_backend=="native_source_cs2" else
                        "python_original_source"),
                    source_callback_rhs_calls=int(self._callback_rhs_calls),
                    source_restart_count=int(self.source_restart_count),
                    source_restart_relative_jump=self.source_restart_relative_jump,
                    knot_waypoints_used=self.knot_waypoints,
                    source_break_waypoint_count=sum(len(p["source_break_s"])
                          for p in self._phases) if self.knot_waypoints else 0,
                    background_backend=self.background.config.get("backend",
                       "Python sharp reference"))
