"""B3 regression: euler_tolerance=5e-4 must admit all source tables.
Fails if a source table's max Euler residual exceeds 5e-4 (a regression
indicating a genuine thermodynamic inconsistency was introduced)."""
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos.api import read_compose

EULER_TOLERANCE = 5e-4


def _max_euler_residual(d):
    n = np.asarray(d.nb, float); mn = float(d.mn)
    Q1 = np.asarray(d.Q1, float); Q3 = np.asarray(d.Q3, float); Q7 = np.asarray(d.Q7, float)
    P = Q1 * n; eps = n * (1.0 + Q7) * mn
    mu_euler = (eps + P) / n; mu_source = (1.0 + Q3) * mn
    return float(np.max(np.abs(mu_source - mu_euler) / np.maximum(np.abs(mu_euler), 1.0)))


def _needs_backend():
    import os
    return os.environ.get("COMPOSE_EOS_FORTRAN_EXE") is None


@pytest.mark.skipif(_needs_backend(), reason="requires the native Fortran backend")
def test_all_source_tables_pass_euler_tolerance():
    """Every source EOS table must pass the v1 euler_tolerance=5e-4 gate.
    The worst observed is RG_KDE0v1 at 1.15e-4 (4.3x margin). A regression
    that pushes any table above 5e-4 indicates a genuine inconsistency."""
    eos_root = ROOT / "eos_generator" / "eos_data"
    worst = 0.0; worst_name = ""
    for d in sorted(eos_root.iterdir()):
        if not d.is_dir() or not (d / "eos.thermo").exists():
            continue
        if d.name in ("GPPVA_DD2_reexpressed", "GPPVA_DD2_matched"):
            continue
        eos = read_compose(d)
        r = _max_euler_residual(eos)
        assert r <= EULER_TOLERANCE, (
            f"{d.name} max Euler residual {r:.3e} exceeds euler_tolerance {EULER_TOLERANCE}; "
            f"genuine thermodynamic inconsistency introduced")
        if r > worst:
            worst = r; worst_name = d.name
    # pin the known worst case as a sanity anchor
    assert worst_name == "RG_KDE0v1", f"worst-case model changed: {worst_name} (expected RG_KDE0v1)"
    assert 1.0e-4 <= worst <= 1.3e-4, f"worst-case residual {worst:.3e} outside expected band"


@pytest.mark.skipif(_needs_backend(), reason="requires the native Fortran backend")
def test_core_tables_pass_at_machine_precision():
    """Core tables (BBB/APR/GM) print Q3 to 17 sig figs and pass the Euler
    relation at machine precision (~1e-16), confirming their providers
    computed Q3 algebraically self-consistent with Q1/Q7."""
    for name in ("BBB_BHF2", "APR_APR", "GM_GM1"):
        eos = read_compose(ROOT / "eos_generator" / "eos_data" / name)
        r = _max_euler_residual(eos)
        assert r < 1e-9, f"core {name} Euler residual {r:.3e} not at machine precision"
