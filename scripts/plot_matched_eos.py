"""Generate 5 physics plots for each surviving matched EOS.

For each surviving css-matched (crust, core) pair, reads the matched CompOSE
EOS and produces 5 plots:
  1. eps vs P        (energy density vs pressure, MeV/fm^3)
  2. nb vs P         (baryon density vs pressure, fm^-3 vs MeV/fm^3)
  3. nb vs cs^2      (baryon density vs sound speed squared, dimensionless)
  4. nb vs mu_b      (baryon density vs baryon chemical potential, MeV)
  5. eps vs mu_b     (energy density vs baryon chemical potential)

Physics (CompOSE convention in this codebase):
  P     = Q1 * nb              (pressure, MeV/fm^3)
  mu_B  = (Q3 + 1) * mn        (baryon chemical potential, MeV)
  eps   = nb * (Q7 + 1) * mn   (energy density, MeV/fm^3)
  cs^2  = dP / d(eps)          (dimensionless, via finite differences)

All plots use a log x-axis where appropriate and mark the crust-core
transition density (the density at which the matching bridge begins).
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import read_compose  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
fm.fontManager.addfont("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

MATCHED_ROOT = ROOT / "download" / "matched"
PLOT_ROOT = ROOT / "download" / "matched_plots"


def compute_cs2(nb: np.ndarray, P: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """cs^2 = dP/d(eps) via central finite differences, clipped to [0, 1]."""
    cs2 = np.zeros_like(nb)
    for i in range(1, len(nb) - 1):
        deps = eps[i + 1] - eps[i - 1]
        if deps > 0:
            cs2[i] = (P[i + 1] - P[i - 1]) / deps
    # forward/backward at endpoints
    if len(nb) >= 2 and (eps[1] - eps[0]) > 0:
        cs2[0] = (P[1] - P[0]) / (eps[1] - eps[0])
    if len(nb) >= 2 and (eps[-1] - eps[-2]) > 0:
        cs2[-1] = (P[-1] - P[-2]) / (eps[-1] - eps[-2])
    return np.clip(cs2, 0.0, 1.0)


def plot_pair(eos, pair_id: str, out_dir: Path) -> list[Path]:
    nb = np.asarray(eos.nb, float)
    Q1 = np.asarray(eos.Q1, float)
    Q3 = np.asarray(eos.Q3, float)
    Q7 = np.asarray(eos.Q7, float)
    mn = float(eos.mn)
    P = Q1 * nb
    mu_b = (Q3 + 1.0) * mn
    eps = nb * (Q7 + 1.0) * mn
    cs2 = compute_cs2(nb, P, eps)

    # sort by nb for clean plots
    order = np.argsort(nb)
    nb, P, mu_b, eps, cs2 = nb[order], P[order], mu_b[order], eps[order], cs2[order]

    plots = []

    def _save(fig, name):
        p = out_dir / name
        fig.savefig(p, dpi=120, bbox_inches="tight")
        plt.close(fig)
        plots.append(p)

    # 1. eps vs P
    fig, ax = plt.subplots(figsize=(5.0, 3.8), constrained_layout=True)
    ax.plot(eps, P, "-", lw=1.6, color="#1f77b4")
    ax.set_xlabel(r"$\varepsilon$ [MeV fm$^{-3}$]")
    ax.set_ylabel(r"$P$ [MeV fm$^{-3}$]")
    ax.set_title(f"{pair_id}: $\\varepsilon$–$P$")
    _save(fig, "eps_P.png")

    # 2. nb vs P
    fig, ax = plt.subplots(figsize=(5.0, 3.8), constrained_layout=True)
    ax.plot(nb, P, "-", lw=1.6, color="#ff7f0e")
    ax.set_xlabel(r"$n_B$ [fm$^{-3}$]")
    ax.set_ylabel(r"$P$ [MeV fm$^{-3}$]")
    ax.set_title(f"{pair_id}: $n_B$–$P$")
    _save(fig, "nb_P.png")

    # 3. nb vs cs^2
    fig, ax = plt.subplots(figsize=(5.0, 3.8), constrained_layout=True)
    ax.plot(nb, cs2, "-", lw=1.6, color="#2ca02c")
    ax.axhline(1.0, ls="--", lw=0.8, color="grey", alpha=0.6)
    ax.axhline(1.0 / 3.0, ls=":", lw=0.8, color="grey", alpha=0.6)
    ax.set_xlabel(r"$n_B$ [fm$^{-3}$]")
    ax.set_ylabel(r"$c_s^2$ [dimensionless]")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title(f"{pair_id}: $n_B$–$c_s^2$")
    _save(fig, "nb_cs2.png")

    # 4. nb vs mu_b
    fig, ax = plt.subplots(figsize=(5.0, 3.8), constrained_layout=True)
    ax.plot(nb, mu_b, "-", lw=1.6, color="#d62728")
    ax.set_xlabel(r"$n_B$ [fm$^{-3}$]")
    ax.set_ylabel(r"$\mu_B$ [MeV]")
    ax.set_title(f"{pair_id}: $n_B$–$\\mu_B$")
    _save(fig, "nb_mub.png")

    # 5. eps vs mu_b
    fig, ax = plt.subplots(figsize=(5.0, 3.8), constrained_layout=True)
    ax.plot(eps, mu_b, "-", lw=1.6, color="#9467bd")
    ax.set_xlabel(r"$\varepsilon$ [MeV fm$^{-3}$]")
    ax.set_ylabel(r"$\mu_B$ [MeV]")
    ax.set_title(f"{pair_id}: $\\varepsilon$–$\\mu_B$")
    _save(fig, "eps_mub.png")

    return plots


def main() -> int:
    if PLOT_ROOT.exists():
        import shutil
        shutil.rmtree(PLOT_ROOT)
    PLOT_ROOT.mkdir(parents=True)

    manifest = json.loads((MATCHED_ROOT / "matching_manifest.json").read_text(encoding="utf-8"))
    survivors = manifest["survivors"]
    print(f"Generating plots for {len(survivors)} surviving EOS...")

    plot_index = []
    for s in survivors:
        pair_id = s["pair_id"]
        compose_dir = Path(s["output_dir"])
        out_dir = PLOT_ROOT / pair_id
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            eos = read_compose(compose_dir)
            paths = plot_pair(eos, pair_id, out_dir)
            plot_index.append({
                "pair_id": pair_id,
                "crust_model": s["crust_model"],
                "core_model": s["core_model"],
                "combined_n": s["combined_n"],
                "plot_dir": str(out_dir),
                "plots": [p.name for p in paths],
            })
            print(f"  {pair_id}: 5 plots saved")
        except Exception as exc:
            print(f"  {pair_id}: FAILED - {type(exc).__name__}: {exc}")
            plot_index.append({"pair_id": pair_id, "error": str(exc)})

    (PLOT_ROOT / "plot_index.json").write_text(
        json.dumps(plot_index, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nPlots saved to {PLOT_ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
