#!/usr/bin/env python3
"""Thin process-local frontend for original-source Fortran TOV/Cowling scans.

All evolution and frequency residuals use native Fortran DVODE. The existing
EOS interpolation/PCHIP setup and presentation remain in the Python frontend.
No multi-thread calls enter DVODE's COMMON-based numerical core.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.stellar_background import ContractEnthalpyEOS
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_cowling_native import NativeCowlingDVODE
from compose_eos.stellar_cowling_tov_native import solve_background_native

SOURCE = ROOT / "data/core_support_20261117/primary_v1/REPRESENTATION.json"
SOURCE_SHA = "591ac27aa99448512e291684c98f8273efc729d17c63be88e331c20e83bd2653"


def run(density: float, destination: Path, count: int) -> None:
    if not np.isfinite(density) or density <= 0:
        raise ValueError("density must be finite positive fm^-3")
    if not 8 <= count <= 4096:
        raise ValueError("native scan point count outside [8,4096]")
    if hashlib.sha256(SOURCE.read_bytes()).hexdigest() != SOURCE_SHA:
        raise RuntimeError("original EOS source identity changed")
    eos = ContractEnthalpyEOS.from_path(SOURCE)
    bg = solve_background_native(
        eos, density, surface_h=1e-4, rtol=1e-10,
        atol=1e-12, max_step_h=.002, native_samples=4097)
    settings = CowlingSettings(
        surface_h=1e-4, rtol=1e-9, atol=1e-12, max_step_h=.002)
    solver = NativeCowlingDVODE(bg, settings, samples=2049)
    freqs = np.linspace(1800., 7200., count)
    brackets, scan = solver.scan_brackets(freqs)
    result = dict(
        schema="ns_cowling.c9.real_source_native_bulk_scan.v1",
        source_sha256=SOURCE_SHA,
        central_density_fm3=float(density),
        background="Fortran DVODE TOV",
        shooting="Fortran DVODE native batch",
        expected_cowling_closure="Gamma1=Gamma_eq; smooth equilibrium barotrope",
        mass_msun=float(bg.mass_msun),
        radius_km=float(bg.radius_km),
        nfrequencies=len(scan),
        frequency_range_hz=[1800.,7200.],
        brackets_hz=brackets,
        scan=scan,
        status="SCOPED_REAL_SMOOTH_STAR_NATIVE_SCAN")
    destination.parent.mkdir(parents=True, exist_ok=True)
    candidate = destination.with_suffix(destination.suffix + ".part")
    candidate.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    os.replace(candidate, destination)
    print("NATIVE_REAL_SCAN",density,"mass",bg.mass_msun,
          "brackets",len(brackets), "records",len(scan),
          "path",destination,flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--central-density",type=float,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--count",type=int,default=241)
    args = ap.parse_args()
    run(args.central_density,args.output,args.count)
