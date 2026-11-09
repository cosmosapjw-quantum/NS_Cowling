"""EOS-preserving selectors for finite-pressure TOV truncation diagnostics.

This module never edits, sorts, interpolates, resamples, or writes EOS rows.
It selects only pressure values already present in a source CSV, in source order.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence


_REQUIRED_COLUMNS = (
    "pressure_pa",
    "energy_mev_fm3",
    "n_b_fm3",
    "mu_b_mev",
)


@dataclass(frozen=True)
class NativeCutoffKnot:
    """One unmodified EOS row selected as a finite-pressure integration limit."""

    row_index: int
    source_row: str
    target_offset_dex: float
    actual_offset_dex: float
    pressure_pa: float
    energy_mev_fm3: float
    n_b_fm3: float
    mu_b_mev: float


def _positive_finite(row: Mapping[str, object], key: str, row_index: int) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"row {row_index + 1} has invalid {key}") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"row {row_index + 1} {key} must be finite and positive")
    return value


def load_native_eos_rows(path: str | Path) -> list[dict[str, str]]:
    """Read rows without reordering or rewriting the source table."""
    source = Path(path)
    with source.open("r", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fields = set(reader.fieldnames or ())
        missing = [key for key in _REQUIRED_COLUMNS if key not in fields]
        if missing:
            raise ValueError(f"missing required EOS columns: {', '.join(missing)}")
        rows = list(reader)
    if len(rows) < 2:
        raise ValueError("native EOS table must contain at least two rows")
    for i, row in enumerate(rows):
        for key in _REQUIRED_COLUMNS:
            _positive_finite(row, key, i)
    return rows


def select_native_pressure_knots(
    rows: Sequence[Mapping[str, object]],
    target_offsets_dex: Sequence[float],
) -> list[NativeCutoffKnot]:
    """Select nearest existing pressure knots relative to the first source row.

    A target offset is only a request for a nearby native knot. The returned
    pressure and every other field are copied from that row exactly after
    numeric parsing; no value is generated at the target pressure.
    """
    if len(rows) < 2:
        raise ValueError("native EOS table must contain at least two rows")
    if not target_offsets_dex:
        raise ValueError("at least one target pressure offset is required")

    targets = [float(x) for x in target_offsets_dex]
    if any(not math.isfinite(x) or x < 0.0 for x in targets):
        raise ValueError("target pressure offsets must be finite and nonnegative")
    if any(b <= a for a, b in zip(targets, targets[1:])):
        raise ValueError("target pressure offsets must be strictly increasing")

    pressure = [_positive_finite(row, "pressure_pa", i) for i, row in enumerate(rows)]
    if any(b <= a for a, b in zip(pressure, pressure[1:])):
        raise ValueError("pressure must be strictly increasing in source order")
    logp = [math.log10(x) for x in pressure]
    edge_logp = logp[0]

    selected: list[NativeCutoffKnot] = []
    used_indices: set[int] = set()
    previous_index = -1
    for target in targets:
        target_logp = edge_logp + target
        if target_logp > logp[-1]:
            raise ValueError("target pressure offset lies outside the native EOS table")
        index = min(range(len(logp)), key=lambda i: abs(logp[i] - target_logp))
        if index in used_indices or index <= previous_index:
            raise ValueError(
                "target pressure offsets must resolve to distinct native knots"
            )
        row = rows[index]
        source_row = str(row.get("source_row", index + 1))
        selected.append(
            NativeCutoffKnot(
                row_index=index,
                source_row=source_row,
                target_offset_dex=target,
                actual_offset_dex=logp[index] - edge_logp,
                pressure_pa=pressure[index],
                energy_mev_fm3=_positive_finite(row, "energy_mev_fm3", index),
                n_b_fm3=_positive_finite(row, "n_b_fm3", index),
                mu_b_mev=_positive_finite(row, "mu_b_mev", index),
            )
        )
        used_indices.add(index)
        previous_index = index
    return selected


@dataclass(frozen=True)
class NativeDensityCutoffKnot:
    """One existing EOS row selected from a low-density native segment."""

    row_index: int
    source_row: str
    target_n_b_fm3: float
    target_log10_n_b_fm3: float
    actual_offset_dex: float
    pressure_pa: float
    energy_mev_fm3: float
    n_b_fm3: float
    mu_b_mev: float


def select_native_density_knots(
    rows: Sequence[Mapping[str, object]],
    *,
    target_count: int = 8,
    maximum_density_fm3: float = 1.0e-3,
) -> list[NativeDensityCutoffKnot]:
    """Select log-spaced density targets using only existing source rows.

    The target grid is a selection rule, not a resampled EOS.  Returned
    pressure, energy density, baryon density and chemical potential are read
    from the selected native row.  The EOS source is never modified.

    The default upper bound isolates the very-low-density segment before the
    first inner-crust-like knots in the current projected cohort.  It is an
    operational cutoff for this diagnostic, not a claim about the physical
    crust-core transition density.
    """
    if len(rows) < 2:
        raise ValueError("native EOS table must contain at least two rows")
    if isinstance(target_count, bool) or int(target_count) != target_count or target_count < 2:
        raise ValueError("target_count must be an integer >= 2")
    maximum_density = float(maximum_density_fm3)
    if not math.isfinite(maximum_density) or maximum_density <= 0.0:
        raise ValueError("maximum_density_fm3 must be finite and positive")

    pressure = [_positive_finite(row, "pressure_pa", i) for i, row in enumerate(rows)]
    density = [_positive_finite(row, "n_b_fm3", i) for i, row in enumerate(rows)]
    if any(b <= a for a, b in zip(pressure, pressure[1:])):
        raise ValueError("pressure must be strictly increasing in source order")
    if any(b <= a for a, b in zip(density, density[1:])):
        raise ValueError("baryon density must be strictly increasing in source order")

    eligible = [i for i, n_b in enumerate(density) if n_b <= maximum_density]
    if len(eligible) < target_count:
        raise ValueError(
            f"low-density segment has {len(eligible)} native rows; "
            f"{target_count} distinct cutoffs were requested"
        )
    if eligible != list(range(eligible[-1] + 1)):
        raise ValueError("low-density native segment is not a contiguous source prefix")

    n_min = density[eligible[0]]
    n_max = density[eligible[-1]]
    log_min, log_max = math.log10(n_min), math.log10(n_max)
    if not log_max > log_min:
        raise ValueError("low-density segment must span a nonzero density interval")

    selected: list[NativeDensityCutoffKnot] = []
    used_indices: set[int] = set()
    previous_index = -1
    for j in range(int(target_count)):
        fraction = j / (int(target_count) - 1)
        target_log = log_min + fraction * (log_max - log_min)
        target_n = 10.0 ** target_log
        index = min(
            eligible,
            key=lambda i: abs(math.log10(density[i]) - target_log),
        )
        if index in used_indices or index <= previous_index:
            raise ValueError(
                "log-density targets do not resolve to distinct increasing native knots"
            )
        row = rows[index]
        selected.append(
            NativeDensityCutoffKnot(
                row_index=index,
                source_row=str(row.get("source_row", index + 1)),
                target_n_b_fm3=target_n,
                target_log10_n_b_fm3=target_log,
                actual_offset_dex=math.log10(density[index] / n_min),
                pressure_pa=pressure[index],
                energy_mev_fm3=_positive_finite(row, "energy_mev_fm3", index),
                n_b_fm3=density[index],
                mu_b_mev=_positive_finite(row, "mu_b_mev", index),
            )
        )
        used_indices.add(index)
        previous_index = index
    return selected
