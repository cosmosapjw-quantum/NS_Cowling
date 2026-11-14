"""Build and invoke the Fortran EOS core through a stable subprocess protocol."""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

class NativeExecutionError(RuntimeError):
    pass

class NativeFeatureUnavailable(NotImplementedError):
    pass

def project_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "fortran_eos" / "build.sh").is_file():
            return parent
    override = os.environ.get("COMPOSE_EOS_FORTRAN_ROOT")
    if override and (Path(override) / "fortran_eos" / "build.sh").is_file():
        return Path(override).resolve()
    installed = Path(sys.prefix) / "share" / "compose-eos"
    if (installed / "fortran_eos" / "build.sh").is_file():
        return installed
    raise NativeExecutionError(
        "Fortran sources were not found. Set COMPOSE_EOS_FORTRAN_ROOT to the project root."
    )

def executable() -> Path:
    root = project_root()
    override = os.environ.get("COMPOSE_EOS_FORTRAN_EXE")
    if override:
        exe = Path(override).expanduser().resolve()
        if exe.is_file():
            return exe
        raise NativeExecutionError(f"COMPOSE_EOS_FORTRAN_EXE does not exist: {exe}")
    exe = root / "build" / "compose_eos_fortran" / "compose_eos_fortran"
    if exe.is_file():
        return exe
    if shutil.which("gfortran") is None:
        raise NativeExecutionError(
            "gfortran is required to build the Fortran EOS backend; no Python numerical fallback is enabled."
        )
    build = subprocess.run(
        ["bash", str(root / "fortran_eos" / "build.sh"), str(exe.parent)],
        cwd=root, capture_output=True, text=True,
    )
    if build.returncode or not exe.is_file():
        raise NativeExecutionError(
            "Fortran EOS backend build failed.\n"
            f"stdout:\n{build.stdout}\nstderr:\n{build.stderr}"
        )
    return exe

def run_native(action: str, *args: str | os.PathLike[str]) -> dict[str, Any]:
    exe = executable()
    command = [str(exe), action, *(str(a) for a in args)]
    completed = subprocess.run(command, capture_output=True, text=True)
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    payload: dict[str, Any] = {}
    for line in reversed(lines):
        try:
            candidate = json.loads(line)
            if isinstance(candidate, dict):
                payload = candidate
                break
        except json.JSONDecodeError:
            continue
    if completed.returncode:
        raise NativeExecutionError(
            f"Fortran operation {action!r} failed with exit code {completed.returncode}: "
            f"{payload.get('message', completed.stderr.strip() or completed.stdout.strip())}"
        )
    if payload.get("status") == "error":
        raise NativeExecutionError(str(payload.get("message", "unknown Fortran error")))
    return payload


def evaluate_ddrmf(model_id: str, rho_b_fm3: float) -> dict[str, Any]:
    """Evaluate a literature-parameterized cold DD-RMF beta-equilibrium state natively.

    The Fortran backend solves the scalar field, rearrangement self-energy,
    charge-neutral npe/mu branch and two explicitly distinct symmetry-energy
    definitions. Python only validates arguments and decodes the native JSON.
    """
    if not isinstance(model_id, str) or not model_id.strip():
        raise ValueError("model_id must be a non-empty string")
    try:
        density = float(rho_b_fm3)
    except (TypeError, ValueError) as exc:
        raise ValueError("rho_b_fm3 must be a finite positive number") from exc
    if not (density > 0.0):
        raise ValueError("rho_b_fm3 must be finite and positive")
    if density != density or density in (float("inf"), float("-inf")):
        raise ValueError("rho_b_fm3 must be finite and positive")
    return run_native("rmf-evaluate", model_id.strip(), format(density, ".17g"))
