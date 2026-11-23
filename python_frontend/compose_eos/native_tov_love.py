"""A pure-input/output Python facade for the standalone Fortran TOV/Love solver.

Numerical integration, EOS interpolation, surface and sharp Maxwell matching
are all done by the executable. Python does not solve any differential equation
or change physical source arrays.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Mapping, Sequence


class NativeTovError(RuntimeError):
    """Native Fortran returned a nonzero exit or a malformed result."""


class NativeTovLove:
    """Run Fortran TOV-Love for the supplied already-materialized P,epsilon CSV.

    Input fields are P and epsilon in MeV/fm^3; output radius km, mass Msun,
    tidal k2/Lambda dimensionless. Only the native executable evaluates EOS
    derivatives or stellar ODEs. No Python fallback exists.
    """

    VALID = frozenset({"dvode_adams", "dvode_bdf", "rk45"})

    def __init__(
        self,
        *,
        binary: str | os.PathLike[str] | None = None,
        integrator: str = "dvode_adams",
        rtol: float = 2e-9,
        timeout_seconds: float = 120.0,
    ) -> None:
        if integrator not in self.VALID:
            raise ValueError("Unknown Fortran integrator; choose Adams, BDF or RK45")
        if not math.isfinite(rtol) or not (1e-12 <= rtol <= 1e-5):
            raise ValueError("Fortran DVODE rtol must be in [1e-12, 1e-5]")
        if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
            raise ValueError("positive finite subprocess deadline required")
        self.integrator = integrator
        self.rtol = rtol
        self.timeout_seconds = timeout_seconds
        candidate = binary or os.environ.get("COMPOSE_EOS_FORTRAN_EXE")
        if not candidate:
            from .native import executable
            candidate = executable()
        self.binary = Path(candidate).resolve()
        if not self.binary.is_file():
            raise FileNotFoundError(self.binary)

    @staticmethod
    def _path(path: str | os.PathLike[str]) -> Path:
        item = Path(path).resolve(strict=True)
        if not item.is_file():
            raise ValueError("native input must be a regular file")
        return item

    def _run(self, command: Sequence[str], *, expected: int) -> list[dict[str, Any]]:
        env = os.environ.copy()
        env["NS_COWLING_TOV_INTEGRATOR"] = self.integrator
        env["NS_COWLING_TOV_RTOL"] = f"{self.rtol:.17e}"
        env["OMP_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        try:
            result = subprocess.run(
                [str(self.binary), "tov-love", *command],
                check=False, capture_output=True, text=True,
                env=env, timeout=self.timeout_seconds,
            )
        except subprocess.TimeoutExpired as e:
            raise NativeTovError("native Fortran TOV exceeded execution deadline") from e
        if result.returncode != 0:
            raise NativeTovError(
                f"native Fortran TOV rejected ({result.returncode}): "
                f"{(result.stdout + result.stderr)[-900:]}"
            )
        try:
            rows = [json.loads(line) for line in result.stdout.splitlines() if line.strip()]
        except json.JSONDecodeError as e:
            raise NativeTovError("native Fortran output is not valid JSON lines") from e
        if len(rows) != expected or any(
            row.get("status") != "ok" or not row.get("converged") for row in rows
        ):
            raise NativeTovError("native Fortran output count or convergence mismatch")
        return rows

    def smooth(self, source: str | os.PathLike[str], central_pressure: float) -> dict[str, Any]:
        if not math.isfinite(central_pressure):
            raise ValueError("finite central pressure required")
        return self._run(
            ["smooth", str(self._path(source)), f"{central_pressure:.17e}"],
            expected=1,
        )[0]

    def layered(
        self, directory: str | os.PathLike[str], central_pressure: float, max_dlogp: float = 0.02,
    ) -> dict[str, Any]:
        """Integrate through labelled CompOSE interfaces to the first pressure row."""
        return self._run(
            ["layered", str(Path(directory).resolve(strict=True)), str(central_pressure), str(max_dlogp)],
            expected=1,
        )[0]

    def layered_target(
        self, directory: str | os.PathLike[str], mass_msun: float, min_pressure: float,
        max_pressure: float, max_dlogp: float = 0.02,
    ) -> dict[str, Any]:
        """Locate the first sampled rising mass crossing in the supplied pressure range."""
        return self._run(
            ["layered-target", str(Path(directory).resolve(strict=True)), str(mass_msun),
             str(min_pressure), str(max_pressure), str(max_dlogp)], expected=1,
        )[0]

    def maxwell(
        self,
        high: str | os.PathLike[str],
        low: str | os.PathLike[str],
        transition_pressure: float,
        delta_epsilon: float,
        central_pressure: float,
    ) -> dict[str, Any]:
        if not all(
            math.isfinite(x) for x in
            (transition_pressure, delta_epsilon, central_pressure)
        ):
            raise ValueError("finite Maxwell transition and central pressure inputs required")
        return self._run(
            ["maxwell", str(self._path(high)), str(self._path(low)),
             f"{transition_pressure:.17e}", f"{delta_epsilon:.17e}",
             f"{central_pressure:.17e}"],
            expected=1,
        )[0]

    def batch_smooth(
        self,
        source: str | os.PathLike[str],
        central_pressures: Sequence[float],
    ) -> list[dict[str, Any]]:
        values = tuple(float(x) for x in central_pressures)
        if not values or len(values) > 10000 or not all(map(math.isfinite, values)):
            raise ValueError("batch requires 1 to 10000 finite central pressures")
        with tempfile.TemporaryDirectory(prefix="ns_native_tov_batch_") as tempdir:
            center_file = Path(tempdir) / "centers.txt"
            center_file.write_text(
                "".join(f"{value:.17e}\n" for value in values),
                encoding="ascii",
            )
            return self._run(
                ["batch-smooth", str(self._path(source)), str(center_file)],
                expected=len(values),
            )
