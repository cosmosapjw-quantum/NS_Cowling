"""Build and invoke the Fortran EOS core through a stable subprocess protocol."""
from __future__ import annotations
import json
import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any

class NativeExecutionError(RuntimeError):
    """A Fortran backend operation failed.

    Carries the structured `category` field emitted by the Fortran layer when
    available (e.g. the maxwell-transition command emits a stable §3/§4
    terminal-outcome identifier alongside the free-form message). This lets
    callers branch on a stable status without parsing exception text.
    """
    def __init__(self, message: str, *, category: str | None = None,
                 operation: str | None = None) -> None:
        super().__init__(message)
        self.category = category
        self.operation = operation

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

# The script owns content hashes, compiler options and cross-process locking.
# This process-local cache skips repeated shell startup only while every input's
# metadata (including ctime) and build artifact remain unchanged. Any edit,
# replacement, chmod or compiler/configuration change reruns --ensure.
_BUILD_LOCK = threading.Lock()
_VERIFIED_BUILDS: dict[Path, tuple] = {}
_BUILD_ENV = ("FC", "FFLAGS", "SLSQP_FFLAGS", "LDFLAGS", "NATIVE_LIBS")


def _build_state(root: Path, exe: Path, compiler: str) -> tuple:
    src = root / "fortran_eos"
    files = [src / "build.sh", src / "run-generator.sh", Path(compiler),
             exe, exe.parent / ".source-signature"]
    for folder in ("src", "app", "vendor"):
        files.extend(p for p in (src / folder).rglob("*") if p.is_file() and
                     (p.suffix in (".f90", ".f") or p.name == "PROVENANCE.json" or
                      "LICENSE" in p.name))
    records = []
    for path in sorted(set(files)):
        try:
            st = path.stat()
            state = (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns,
                     st.st_ctime_ns, st.st_mode)
        except FileNotFoundError:
            state = None
        records.append((str(path), state))
    return tuple((name, os.environ.get(name)) for name in _BUILD_ENV), tuple(records)


def executable() -> Path:
    override = os.environ.get("COMPOSE_EOS_FORTRAN_EXE")
    if override:
        exe = Path(override).expanduser().resolve()
        if exe.is_file() and os.access(exe, os.X_OK):
            return exe
        raise NativeExecutionError(f"COMPOSE_EOS_FORTRAN_EXE is not an executable file: {exe}")
    root = project_root()
    exe = root / "build" / "compose_eos_fortran" / "compose_eos_fortran"
    compiler = shutil.which(os.environ.get("FC", "gfortran"))
    if compiler is None:
        raise NativeExecutionError(
            "A Fortran compiler is required to verify the source build; "
            "use COMPOSE_EOS_FORTRAN_EXE only for an explicitly selected prebuilt executable."
        )
    with _BUILD_LOCK:
        state = _build_state(root, exe, compiler)
        if (_VERIFIED_BUILDS.get(exe) == state and exe.is_file() and
                os.access(exe, os.X_OK)):
            return exe
        build = subprocess.run(
            ["bash", str(root / "fortran_eos" / "build.sh"), "--ensure", str(exe.parent)],
            cwd=root, capture_output=True, text=True,
        )
        if build.returncode or not exe.is_file() or not os.access(exe, os.X_OK):
            raise NativeExecutionError(
                "Fortran EOS backend build/verification failed.\n"
                f"stdout:\n{build.stdout}\nstderr:\n{build.stderr}"
            )
        _VERIFIED_BUILDS[exe] = _build_state(root, exe, compiler)
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
    # the structured category the Fortran layer attaches to structured failures
    # (e.g. maxwell-transition). None when the Fortran layer did not emit one.
    category = payload.get("category")
    op = payload.get("operation", action)
    if completed.returncode:
        raise NativeExecutionError(
            f"Fortran operation {action!r} failed with exit code {completed.returncode}: "
            f"{payload.get('message', completed.stderr.strip() or completed.stdout.strip())}",
            category=category, operation=op,
        )
    if payload.get("status") == "error":
        raise NativeExecutionError(
            str(payload.get("message", "unknown Fortran error")),
            category=category, operation=op,
        )
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
