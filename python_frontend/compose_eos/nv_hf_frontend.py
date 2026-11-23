"""Process-isolated frontend for the native NV orbital Hartree-Fock solver.

This module transports a declared experiment and validates process/result
identity. It contains no EOS, nuclear-functional, or minimization arithmetic.
An accepted process result is a fixed-composition cell, not an admitted EOS.
"""
from __future__ import annotations

import hashlib
import math
import os
from pathlib import Path
import subprocess
import time


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_native_cell(path: Path) -> dict:
    """Read the native receipt and verify its declared profile shape."""
    fields: dict[str, str] = {}
    profiles = 0
    last_radius = 0.0
    profile_started = False
    for line in Path(path).read_text().splitlines():
        if line.startswith("# r_fm "):
            if profile_started:
                raise ValueError("Duplicate native profile header")
            profile_started = True
            continue
        if not line or line.startswith("#"):
            continue
        if profile_started:
            row = line.split()
            if len(row) != 7 or not all(math.isfinite(float(v)) for v in row):
                raise ValueError("Invalid or nonfinite native radial profile")
            if float(row[0]) <= last_radius or float(row[1]) <= 0:
                raise ValueError("Invalid native radial geometry")
            last_radius = float(row[0])
            profiles += 1
            continue
        key, separator, value = line.partition("=")
        if not separator or key in fields:
            raise ValueError("Malformed or duplicate native receipt field")
        fields[key] = value.strip()
    for key in ("status", "converged", "spectral_complete", "neutrons", "protons", "mesh", "seed_kind", "coefficient_mode", "nb_fm3"):
        if key not in fields:
            raise ValueError("Missing native receipt field: " + key)
    if fields["status"] == "ok":
        for key in ("iterations", "radius_fm", "tau_residual", "spin_residual", "selfconsistent_residual", "energy_per_baryon_MeV", "pressure_MeV_fm3",
                    "mu_n_kin_MeV", "mu_p_kin_MeV", "mu_e_total_MeV", "beta_residual_MeV",
                    "density_residual", "energy_residual_MeV", "particle_error", "dilation_pressure_error_MeV_fm3"):
            if key not in fields:
                raise ValueError("Successful native receipt missing required diagnostic: " + key)
    for key in ("converged", "spectral_complete", "outer_smoothing", "energy_valid"):
        if fields.get(key) not in ("T", "F"):
            raise ValueError("Invalid native logical field: " + key)
    for key, value in fields.items():
        if key in ("status", "converged", "spectral_complete", "outer_smoothing", "energy_valid"):
            continue
        if not math.isfinite(float(value)):
            raise ValueError("Nonfinite native receipt: " + key)
    if profile_started and profiles != int(fields["mesh"]):
        raise ValueError("Native profile row count differs from mesh")
    if fields["status"] == "ok" and (not profile_started or fields["converged"] != "T" or fields["spectral_complete"] != "T" or fields["energy_valid"] != "T"):
        raise ValueError("Native success lacks a complete converged orbital profile")
    return {"fields": fields, "profile_rows": profiles}


def run_native_cell(executable: Path, source: Path, output: Path, config: dict, *, timeout: float = 300) -> dict:
    """Run one explicit (N,Z,n_b) configuration with one native thread."""
    allowed = {"neutrons", "protons", "nb", "mesh", "max_iterations", "mixing", "seed_kind", "outer_smoothing", "coefficient_mode"}
    if set(config) != allowed:
        raise ValueError("HF configuration must declare exactly: " + ", ".join(sorted(allowed)))
    for key in ("neutrons", "protons", "mesh", "max_iterations", "seed_kind", "coefficient_mode"):
        if type(config[key]) is not int:
            raise ValueError("HF integer input required: " + key)
    for key in ("outer_smoothing",):
        if type(config[key]) is not bool:
            raise ValueError("HF logical input must be boolean: " + key)
    if any(not math.isfinite(float(config[key])) for key in ("nb", "mixing")):
        raise ValueError("HF numeric input must be finite")
    executable, source, output = map(lambda p: Path(p).resolve(), (executable, source, output))
    if output.exists() or output.with_suffix(".log").exists():
        raise FileExistsError("Refusing to overwrite frozen HF output: " + str(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    before = {str(p): sha256(p) for p in (executable, source)}
    command = [str(executable), str(source), str(output), *[
        str(config[key]) for key in ("neutrons", "protons", "nb", "mesh", "max_iterations", "mixing", "seed_kind")
    ], ".true." if config["outer_smoothing"] else ".false.", str(config["coefficient_mode"])]
    env = dict(os.environ)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
        env[key] = "1"
    start = time.perf_counter()
    try:
        process = subprocess.run(command, text=True, capture_output=True, timeout=timeout, env=env)
        log = process.stdout + process.stderr
        returncode = process.returncode
    except subprocess.TimeoutExpired as exc:
        log = str(exc.stdout or "") + str(exc.stderr or "")
        returncode = None
    output.with_suffix(".log").write_text(log)
    after = {str(p): sha256(p) for p in (executable, source)}
    result = dict(command=command, config=config, returncode=returncode,
                  wall_seconds=time.perf_counter() - start, inputs_before=before,
                  inputs_after=after, inputs_unchanged=before == after,
                  output_path=str(output), log_path=str(output.with_suffix(".log")),
                  status="RUNTIME_FAILURE", source_eos_admitted=False)
    if output.is_file():
        result["output_sha256"] = sha256(output)
        try:
            parsed = read_native_cell(output)
            result["native"] = parsed
            fields = parsed["fields"]
            for key in ("neutrons", "protons", "mesh", "seed_kind", "coefficient_mode"):
                if int(fields[key]) != config[key]:
                    raise ValueError("Native configuration identity mismatch: " + key)
            if float(fields["nb_fm3"]) != float(config["nb"]):
                raise ValueError("Native density identity mismatch")
            for key in ("outer_smoothing",):
                if (fields[key] == "T") != config[key]:
                    raise ValueError("Native logical configuration identity mismatch: " + key)
            if returncode == 0 and fields["status"] == "ok":
                result["status"] = "FIXED_COMPOSITION_CELL_CONVERGED"
            elif fields["status"] == "HF_ITERATION_NOT_CONVERGED":
                result["status"] = "NUMERICAL_SCF_NOT_CONVERGED"
            else:
                result["status"] = "NATIVE_CELL_REJECTED"
        except (ValueError, OverflowError) as exc:
            result["status"] = "PROTOCOL_FAILURE"
            result["error"] = str(exc)
    if returncode is None:
        result["status"] = "RUNTIME_TIMEOUT"
    if before != after:
        result["status"] = "INPUT_IDENTITY_CHANGED"
    return result
