"""Opt-in Fortran DVODE backend for smooth relativistic Cowling modes.

Python owns EOS/TOV source provenance and frontend marshaling. Shooting ODE,
frequency bracketing/root polishing and profiles use native Fortran DVODE.
The original SciPy brentq reference remains explicitly selectable.
A sharp Maxwell interface is *not* supported by this smooth backend.
DVODE's global COMMON workspace is serialized by one process-wide lock.
"""
from __future__ import annotations

import ctypes as ct
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
import subprocess
import tempfile
import threading

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from .stellar_cowling import C_KM_S, CowlingSettings

_ROOT = Path(__file__).resolve().parents[2]
_FORTRAN = _ROOT / "fortran_eos"
_SOURCES = [
    _FORTRAN / "src/compose_eos_cowling_dvode.f90",
    _FORTRAN / "src/compose_eos_cowling_pchip.f90",
    _FORTRAN / "src/compose_eos_cowling_source_exact.f90",
    _FORTRAN / "src/compose_eos_cowling_sharp_root.f90",
    _FORTRAN / "src/compose_eos_cowling_sharp_source_root.f90",
    _FORTRAN / "src/compose_eos_cowling_tov_dvode.f90",
    _FORTRAN / "vendor/odepack/vode.f",
    _FORTRAN / "vendor/odepack/dgefa.f",
    _FORTRAN / "vendor/odepack/dgesl.f",
    _FORTRAN / "vendor/odepack/dgbfa.f",
    _FORTRAN / "vendor/odepack/dgbsl.f",
]
_LOCK = threading.RLock()
_DLL = None
_ACTIVE = None
_DOUBLE_ARRAY = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
_INT_ARRAY = np.ctypeslib.ndpointer(dtype=np.int32,flags="C_CONTIGUOUS")
_FORTRAN_ARRAY = np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")


def _native_library():
    """Compile the independent native source once, keyed by byte identities."""
    global _DLL
    if _DLL is not None:
        return _DLL
    missing = [str(p) for p in _SOURCES if not p.is_file()]
    if missing:
        raise RuntimeError("Fortran Cowling source unavailable: " + str(missing))
    key = sha256(b"".join(p.read_bytes() for p in _SOURCES)
                 + b"gfortran -O2 -fPIC -lblas").hexdigest()[:18]
    directory = Path(tempfile.gettempdir()) / "ns_cowling_native_dvode"
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / ("libcowling_" + key + ".so")
    if not target.exists():
        temporary = target.with_name(target.name + ".tmp")
        cmd = ["gfortran", "-O2", "-fPIC", "-shared", "-ffree-line-length-none",
               *map(str, _SOURCES), "-lblas", "-o", str(temporary)]
        run = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if run.returncode:
            raise RuntimeError("Fortran DVODE compiler failure: " + run.stderr[-3000:])
        temporary.replace(target)
    dll = ct.CDLL(str(target))
    dll.cowling_native_pchip_slopes.argtypes = [
        ct.c_int,_DOUBLE_ARRAY,_FORTRAN_ARRAY,_FORTRAN_ARRAY,
        ct.POINTER(ct.c_int)]
    dll.cowling_native_pchip_slopes.restype = None
    dll.cowling_native_pchip_values.argtypes = [
        ct.c_int,_DOUBLE_ARRAY,_FORTRAN_ARRAY,_FORTRAN_ARRAY,
        ct.c_int,_DOUBLE_ARRAY,_FORTRAN_ARRAY,ct.POINTER(ct.c_int)]
    dll.cowling_native_pchip_values.restype = None
    dll.cowling_native_setup.argtypes = [
        ct.c_int, _DOUBLE_ARRAY, _FORTRAN_ARRAY, _FORTRAN_ARRAY,
        ct.POINTER(ct.c_int)]
    dll.cowling_native_setup.restype = None
    dll.cowling_native_shot.argtypes = [
        ct.c_double, ct.c_int, ct.c_double, ct.c_double, ct.c_double,
        ct.c_int, ct.c_int, _DOUBLE_ARRAY, ct.c_int, _DOUBLE_ARRAY,
        _FORTRAN_ARRAY, _FORTRAN_ARRAY,
        ct.POINTER(ct.c_double), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
    dll.cowling_native_shot.restype = None
    dll.cowling_native_scan.argtypes = [
        ct.c_int, _DOUBLE_ARRAY, ct.c_int, ct.c_double, ct.c_double,
        ct.c_double, ct.c_int, ct.c_double, ct.c_double, _DOUBLE_ARRAY,
        _INT_ARRAY, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
    dll.cowling_native_scan.restype = None
    dll.cowling_native_root.argtypes = [
        ct.c_double, ct.c_double, ct.c_int, ct.c_double, ct.c_double,
        ct.c_double, ct.c_int, ct.c_double, ct.c_double, ct.c_double,
        ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
    dll.cowling_native_root.restype = None
    dll.cowling_native_sharp_onejump_brent.argtypes = [
        ct.c_double,ct.c_double,ct.c_int,ct.c_double,ct.c_double,
        ct.c_double,ct.c_int,
        ct.c_int,_DOUBLE_ARRAY,_FORTRAN_ARRAY,_FORTRAN_ARRAY,
        ct.c_int,_DOUBLE_ARRAY,_FORTRAN_ARRAY,_FORTRAN_ARRAY,
        ct.c_double,_DOUBLE_ARRAY,ct.c_int,ct.c_double,
        ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_int),ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    dll.cowling_native_sharp_onejump_brent.restype = None
    dll.cowling_native_sharp_source_brent.argtypes = [
        ct.c_double,ct.c_double,ct.c_int,ct.c_double,ct.c_double,
        ct.c_double,ct.c_int,
        ct.c_int,_DOUBLE_ARRAY,_FORTRAN_ARRAY,_FORTRAN_ARRAY,
        ct.c_int,_DOUBLE_ARRAY,_FORTRAN_ARRAY,_FORTRAN_ARRAY,
        ct.c_double,ct.c_double,_DOUBLE_ARRAY,ct.c_int,
        ct.c_double,ct.c_double,
        ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_int),ct.POINTER(ct.c_int),
        ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    dll.cowling_native_sharp_source_brent.restype = None
    dll.cowling_native_source_eos_setup.argtypes = [
        ct.c_int,_FORTRAN_ARRAY,ct.c_int,_FORTRAN_ARRAY,
        _DOUBLE_ARRAY,_DOUBLE_ARRAY,ct.c_int,ct.POINTER(ct.c_int)]
    dll.cowling_native_source_eos_setup.restype = None
    dll.cowling_native_source_cs2_eval.argtypes = [
        ct.c_double,ct.c_int,ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    dll.cowling_native_source_cs2_eval.restype = None
    dll.cowling_native_source_select_restarts.argtypes = [
        ct.c_int,ct.c_double,ct.c_double,ct.c_double,ct.c_double,
        ct.c_int,_DOUBLE_ARRAY,ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    dll.cowling_native_source_select_restarts.restype = None
    dll.cowling_native_source_phase_setup.argtypes = [
        ct.c_int,_DOUBLE_ARRAY,_FORTRAN_ARRAY,_FORTRAN_ARRAY,
        ct.c_double,ct.c_int,ct.POINTER(ct.c_int)]
    dll.cowling_native_source_phase_setup.restype = None
    dll.cowling_native_source_register.argtypes = [ct.POINTER(ct.c_int)]
    dll.cowling_native_source_register.restype = None
    dll.cowling_native_source_unregister.argtypes = [ct.POINTER(ct.c_int)]
    dll.cowling_native_source_unregister.restype = None
    dll.cowling_native_segment.argtypes = [
        ct.c_double,ct.c_int,ct.c_double,ct.c_double,ct.c_double,
        ct.c_int,ct.c_int,_DOUBLE_ARRAY,_DOUBLE_ARRAY,_FORTRAN_ARRAY,
        ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    dll.cowling_native_segment.restype = None
    dll.cowling_native_segment_restarts.argtypes = [
        ct.c_double,ct.c_int,ct.c_double,ct.c_double,ct.c_double,
        ct.c_int,ct.c_int,_DOUBLE_ARRAY,ct.c_int,_INT_ARRAY,
        _DOUBLE_ARRAY,_FORTRAN_ARRAY,ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    dll.cowling_native_segment_restarts.restype = None
    dll.cowling_native_set_source_callback.argtypes = [
        ct.c_void_p,ct.POINTER(ct.c_int)]
    dll.cowling_native_set_source_callback.restype = None
    dll.cowling_native_interface_transfer.argtypes = [
        ct.c_double,ct.c_double,ct.c_double,ct.c_double,ct.c_double,
        ct.c_double,ct.c_double,ct.c_double,ct.c_int,_DOUBLE_ARRAY,
        _DOUBLE_ARRAY,ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),
        ct.POINTER(ct.c_int)]
    dll.cowling_native_interface_transfer.restype = None
    dll.tov_native_profile.argtypes = [
        ct.c_int, _DOUBLE_ARRAY, _DOUBLE_ARRAY, _DOUBLE_ARRAY,
        _DOUBLE_ARRAY, _DOUBLE_ARRAY,
        ct.c_double, ct.c_double, ct.c_double, ct.c_double,
        ct.c_double, ct.c_double, ct.c_int, _FORTRAN_ARRAY,
        ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
    dll.tov_native_profile.restype = None
    _DLL = dll
    return dll


def _nodes(a):
    a = np.asarray(a)
    threshold = 1e-9*np.max(abs(a))
    a = a[abs(a)>threshold]
    return int(np.count_nonzero(a[:-1]*a[1:]<0))


class NativeCowlingDVODE:
    """A fixed smooth background, with native eigenfunction shots.

    Coefficient rows are (r_km, Phi, lambda_metric, dPhi_dr_km,
    cs2_over_c2) and are PCHIP-interpolated in s=sqrt(h_start-h).
    The tabulation is an approximation to the Python exact-contract
    background; test sample resolution before promoting results.
    """

    def __init__(self, background, settings=CowlingSettings(), *,
                 samples=8193, method="adams"):
        settings.validate(background)
        jump = bool(getattr(background, "sharp_interfaces", False))
        barotrope = getattr(getattr(background, "eos", None), "barotrope", None)
        jump = jump or bool(getattr(barotrope, "sharp_interfaces", False))
        if jump:
            raise ValueError("sharp Maxwell phase interfaces require separate jump matching")
        if samples < 33:
            raise ValueError("samples must be at least 33")
        if method not in ("adams", "bdf"):
            raise ValueError("native DVODE method must be adams or bdf")
        self.background = background
        self.settings = settings
        self.method = method
        self.mf = {"adams": 10, "bdf": 22}[method]
        self.surface_s = float(np.sqrt(background.h_start-settings.surface_h))
        self.match_s = float(np.sqrt(background.h_start
                                     -settings.match_h_fraction*background.h_c))
        base_s = np.linspace(0, self.surface_s, samples)
        # Log-enthalpy clustering resolves the near-surface sound-speed
        # behavior without needing O(h_start/surface_h) uniform s knots.
        if settings.surface_h > 0.:
            hg=np.geomspace(settings.surface_h, background.h_start,
                             max(1025, samples//2))
            base_s=np.unique(np.r_[base_s,
                                    np.sqrt(np.maximum(background.h_start-hg,0.))])
        self.break_s = np.array([], dtype=float)
        if hasattr(background, "break_h"):
            break_h = np.asarray(background.break_h, dtype=float)
            break_h = break_h[(break_h > settings.surface_h)
                              & (break_h < background.h_start)]
            self.break_s = np.sort(np.sqrt(background.h_start-break_h))
        # EOS-polytrope boundaries are only C0 in sound speed. One-sided
        # samples avoid smearing O(1) jumps across a global background cell.
        ds = self.surface_s/(samples-1)
        if len(self.break_s):
            widths = np.minimum(0.01*ds,
                                np.maximum(1e-10, 2e-12/np.maximum(self.break_s,1e-4)))
            neighbors = np.r_[self.break_s-widths,self.break_s,
                              self.break_s+widths]
            neighbors = neighbors[(neighbors>0)&(neighbors<self.surface_s)]
            self.s = np.unique(np.r_[base_s,neighbors])
        else:
            self.s = base_s
        hvalues = background.h_start-self.s*self.s
        hvalues[0] = background.h_start
        hvalues[-1] = settings.surface_h
        records = []
        if hasattr(background, "coefficients_on_s"):
            vals = np.asarray(background.coefficients_on_s(self.s), dtype=np.float64)
            if vals.shape != (len(self.s), 5):
                raise ValueError("native background coefficient shape mismatch")
        else:
            for h in hvalues:
                st = background.state_at_enthalpy(float(h))
                cs = float(st["cs2_over_c2"])
                if np.isposinf(cs):
                    cs = 1e100
                records.append([st["r_km"], st["Phi"],
                                st["lambda_metric"], st["dPhi_dr_km"], cs])
            vals = np.asarray(records, dtype=np.float64)
        if not np.all(np.isfinite(vals)):
            raise ValueError("nonfinite background coefficients")
        self._pchip = PchipInterpolator(self.s, vals, axis=0, extrapolate=False)
        self._s = np.ascontiguousarray(self.s)
        self._values = np.asfortranarray(vals.T)
        self._slopes = np.asfortranarray(self._pchip(self.s, nu=1).T)
        self._lib = None

    def _setup(self):
        global _ACTIVE
        dll = _native_library()
        if _ACTIVE is not self:
            status = ct.c_int(-1)
            dll.cowling_native_setup(len(self._s), self._s, self._values,
                                     self._slopes, ct.byref(status))
            if status.value != 0:
                raise RuntimeError(f"Fortran Cowling background rejected: {status.value}")
            _ACTIVE = self
        self._lib = dll

    def _shoot(self, frequency_hz, ni=2, no=2):
        if not np.isfinite(frequency_hz) or frequency_hz <= 0:
            raise ValueError("frequency_hz must be finite and positive")
        if ni < 2 or no < 2:
            raise ValueError("each native trajectory requires two or more points")
        si = np.ascontiguousarray(np.linspace(0., self.match_s, ni))
        so = np.ascontiguousarray(np.linspace(self.surface_s, self.match_s, no))
        ui = np.empty((2,ni), dtype=float, order="F")
        uo = np.empty((2,no), dtype=float, order="F")
        residual, nfev, status = ct.c_double(), ct.c_int(), ct.c_int()
        omega_geom = 2*np.pi*float(frequency_hz)/C_KM_S
        # At the outer end |dh/ds|=2s. Cap native s-steps so that
        # the maximum h-step cannot exceed the reference max_step_h.
        max_ds = min(self.surface_s/200,
                     self.settings.max_step_h/(2*self.surface_s))
        with _LOCK:
            self._setup()
            self._lib.cowling_native_shot(
                omega_geom, self.settings.ell,
                self.settings.rtol, self.settings.atol, max_ds,
                self.mf, ni, si, no, so, ui, uo,
                ct.byref(residual), ct.byref(nfev), ct.byref(status))
        if status.value:
            raise RuntimeError("Fortran DVODE Cowling shot rejected: "
                               +str(status.value))
        return si, so, ui, uo, float(residual.value), int(nfev.value)

    def matching_residual(self, frequency_hz):
        return self._shoot(frequency_hz)[4]

    def scan_brackets(self, frequencies_hz):
        """Frontend marshaling only; Fortran evaluates all frequency shots."""
        fs = np.ascontiguousarray(np.asarray(frequencies_hz, dtype=np.float64))
        if (fs.ndim != 1 or len(fs) < 2 or not np.all(np.isfinite(fs))
            or fs[0] <= 0 or np.any(np.diff(fs) <= 0)):
            raise ValueError("scan frequencies must be strictly increasing and positive")
        residuals = np.empty(len(fs), dtype=np.float64)
        crossings = np.empty(len(fs)-1, dtype=np.int32)
        calls, status = ct.c_int(0), ct.c_int(-1)
        max_ds = min(self.surface_s/200,
                     self.settings.max_step_h/(2*self.surface_s))
        with _LOCK:
            self._setup()
            self._lib.cowling_native_scan(
                len(fs), fs, self.settings.ell, self.settings.rtol,
                self.settings.atol, max_ds, self.mf,
                self.match_s, self.surface_s, residuals, crossings,
                ct.byref(calls), ct.byref(status))
        if status.value:
            raise RuntimeError("Fortran DVODE bulk scan rejected: "
                               +str(status.value))
        return ([(float(fs[i]),float(fs[i+1])) for i in range(len(fs)-1)
                 if crossings[i]],
                [dict(frequency_hz=float(f),residual=float(v))
                 for f,v in zip(fs,residuals)])

    def solve_mode(self, bracket_hz, profile_points=1200, *,
                   root_backend="fortran"):
        lo,hi = map(float,bracket_hz)
        if not (0 < lo < hi and np.isfinite(lo) and np.isfinite(hi)):
            raise ValueError("positive ordered frequency bracket required")
        if profile_points < 20:
            raise ValueError("profile_points must be at least 20")
        if root_backend not in ("fortran", "python_brentq"):
            raise ValueError("unknown Cowling root_backend")
        if root_backend == "python_brentq":
            root, info = brentq(self.matching_residual,lo,hi,
                               xtol=self.settings.root_xtol_hz,
                               rtol=4*np.finfo(float).eps, full_output=True)
            root_iterations = int(info.iterations)
            root_converged = bool(info.converged)
        else:
            root_out, residual_out = ct.c_double(0), ct.c_double(0)
            iterations = ct.c_int(0)
            calls, status = ct.c_int(0), ct.c_int(-1)
            max_ds = min(self.surface_s/200,
                         self.settings.max_step_h/(2*self.surface_s))
            with _LOCK:
                self._setup()
                self._lib.cowling_native_root(
                    lo,hi,self.settings.ell,self.settings.rtol,
                    self.settings.atol,max_ds,self.mf,self.match_s,
                    self.surface_s,self.settings.root_xtol_hz,
                    ct.byref(root_out),ct.byref(residual_out),
                    ct.byref(iterations),ct.byref(calls),ct.byref(status))
            if status.value == 3:
                raise ValueError("native Cowling root bracket lacks sign change")
            if status.value:
                raise RuntimeError("Fortran DVODE root polishing failed: "
                                   +str(status.value))
            root = float(root_out.value)
            root_iterations = int(iterations.value)
            root_converged = True
        si,so,inner,outer,residual,nfev = self._shoot(
            root, profile_points//2, profile_points//2)
        at_match = inner[:,-1]
        outer_match = outer[:,-1]
        scale = np.dot(at_match,outer_match)/np.dot(outer_match,outer_match)
        s_path = np.r_[si,so[::-1][1:]]
        states = np.column_stack((inner,(outer*scale)[:,::-1][:,1:]))
        radius = self._pchip(s_path)[:,0]
        ell = self.settings.ell
        W = radius**(ell+1)*states[0]
        V = radius**ell*states[1]
        norm = np.max(abs(W))
        if not np.isfinite(norm) or norm <= 0:
            raise RuntimeError("zero/nonfinite native eigenfunction")
        h_path = self.background.h_start-s_path*s_path
        profile = [dict(h=float(h),r_km=float(r),W=float(w/norm),
                        V_km_inverse=float(v/norm))
                   for h,r,w,v in zip(h_path,radius,W,V)]
        # The surface BC is imposed by the native initial condition; retain
        # an independent endpoint residual check of the returned eigenvector.
        r,phi,lam,gp,_ = self._values[:,-1]
        omega = 2*np.pi*root/C_KM_S
        terms = np.array([omega**2*np.exp(lam-2*phi)*r*states[1,-1],
                          gp*states[0,-1]])
        return dict(
            frequency_hz=float(root), ell=ell, bracket_hz=[lo,hi],
            root_converged=root_converged,root_iterations=root_iterations,
            root_backend=("fortran_bracketed_brent" if root_backend == "fortran"
                          else "python_brentq"),
            matching_residual=float(residual),
            relative_join_mismatch=float(np.linalg.norm(at_match-scale*outer_match)
                                         /np.linalg.norm(at_match)),
            surface_BC_relative=float(abs(sum(terms))/sum(abs(terms))),
            W_nodes=_nodes(W),V_nodes=_nodes(V),
            r_surface_mode_km=float(radius[-1]),
            profile=profile,nfev=nfev,backend="Fortran DVODE "+self.method.capitalize(),
            background_interpolation="PCHIP(sqrt(h_start-h))",
            background_samples=len(self.s))
