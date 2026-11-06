"""Diagnostic figures for a cold EOS table.

All figures are written as vector PDF (suitable for inclusion in a LaTeX
article) and, optionally, as raster PNG.  The routine :func:`make_figures`
produces the standard set used in the accompanying report.
"""
from __future__ import annotations

import os
from typing import Optional

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .eos_data import EOSData
from .consistency import check_consistency

# professional colour palette
_C = dict(orig="#4c72b0", recon="#dd8452", fill="#55a868",
          gap="#d8d8d8", accent="#c44e52", purple="#8172b3",
          orange="#ccb974", red="#937860")


def _style(ax):
    ax.tick_params(direction="in", which="both", top=True, right=True)
    ax.grid(True, ls=":", alpha=0.4)


def plot_pressure(data: EOSData, ax=None, gap_index: Optional[int] = None,
                  show_fill: bool = False, gap_indices=None):
    """Pressure versus baryon density; support one or many gap-edge indices.

    ``gap_index`` remains accepted for backwards compatibility.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    ax.loglog(data.nb, data.p, ".", ms=2.5, color=_C["orig"], label="table")
    indices = [] if gap_indices is None else list(gap_indices)
    if gap_index is not None:
        indices.append(gap_index)
    for idx in sorted(set(int(i) for i in indices)):
        if idx < 0 or idx >= data.n - 1:
            raise IndexError(f"gap edge index {idx} is outside [0, {data.n - 2}]")
        ax.axvspan(data.nb[idx], data.nb[idx + 1],
                   alpha=0.18, color=_C["gap"])
        ax.axvline(data.nb[idx], ls=":", color="grey", lw=0.7)
        ax.axvline(data.nb[idx + 1], ls=":", color="grey", lw=0.7)
    ax.set_xlabel(r"$n_b\;[\mathrm{fm}^{-3}]$")
    ax.set_ylabel(r"$p\;[\mathrm{MeV\,fm}^{-3}]$")
    _style(ax)
    return ax


def plot_sound_speed(data: EOSData, ax=None, report: Optional[dict] = None):
    """Smooth-branch sound speed squared, exposing both sub- and superluminal values."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    report = report or check_consistency(data, verbose=False)
    cs2 = np.asarray(report["cs2"], dtype=float).copy()
    mask = np.asarray(report["derivative_mask"], dtype=bool)
    cs2[~mask] = np.nan
    ax.semilogx(data.nb, cs2, "-", color=_C["orig"], lw=1.0)
    ax.axhline(1.0, ls="--", color=_C["accent"], lw=0.8, label="causal limit")
    ax.axhline(1 / 3, ls=":", color="grey", lw=0.7)
    ax.set_xlabel(r"$n_b\;[\mathrm{fm}^{-3}]$")
    ax.set_ylabel(r"$c_s^2$")
    finite = cs2[np.isfinite(cs2)]
    if finite.size:
        ymin = min(-0.05, float(np.min(finite)) * 1.08)
        ymax = max(1.05, float(np.max(finite)) * 1.08)
        ax.set_ylim(ymin, ymax)
    _style(ax)
    return ax


def plot_composition(data: EOSData, ax=None):
    """Particle fractions versus baryon density."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    ax.semilogx(data.nb, data.Ye, "-", color=_C["accent"], lw=1, label=r"$Y_e$")
    ax.semilogx(data.nb, data.Yp, "-", color=_C["orig"], lw=1, label=r"$Y_p$")
    ax.semilogx(data.nb, data.Ymu, "-", color=_C["purple"], lw=1, label=r"$Y_\mu$")
    ax.semilogx(data.nb, data.Yn_eff, "-", color=_C["fill"], lw=1,
                label=r"$Y_n$ (free)")
    ax.set_xlabel(r"$n_b\;[\mathrm{fm}^{-3}]$")
    ax.set_ylabel("fraction")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=7, ncol=2, loc="center left")
    _style(ax)
    return ax


def plot_gibbs_duhem(report: dict, ax=None):
    """Gibbs-Duhem residual versus baryon density."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    nb = report["nb"]
    ax.semilogx(nb, np.abs(report["gd_residual"]), "-", color=_C["orig"], lw=1)
    ax.axhline(1e-4, ls="--", color="grey", lw=0.7)
    ax.set_xlabel(r"$n_b\;[\mathrm{fm}^{-3}]$")
    ax.set_ylabel(r"$|\Delta_{\rm GD}|/|\mu_b|$")
    ax.set_ylim(1e-10, 1e-2)
    _style(ax)
    return ax


def plot_mass_radius(mr_path: str, ax=None):
    """Neutron-star mass-radius relation from a CompOSE ``eos.mr`` file."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 3.4), constrained_layout=True)
    mr = np.loadtxt(mr_path, ndmin=2)
    ax.plot(mr[:, 0], mr[:, 1], "-", color=_C["orig"], lw=1.5)
    imax = mr[:, 1].argmax()
    ax.scatter([mr[imax, 0]], [mr[imax, 1]], color=_C["accent"], s=25, zorder=5)
    ax.axhline(1.4, ls=":", color="grey", lw=0.7)
    ax.axhline(2.0, ls=":", color="grey", lw=0.7)
    ax.set_xlabel("radius [km]")
    ax.set_ylabel(r"gravitational mass $[M_\odot]$")
    _style(ax)
    return ax


def make_figures(data: EOSData, outdir: str,
                 mr_path: Optional[str] = None,
                 formats: tuple = ("pdf",)) -> list:
    """Produce the standard figure set and return the list of file paths.

    Parameters
    ----------
    data : EOSData
        The (reconstructed) table to plot.
    outdir : str
        Output directory (created if absent).
    mr_path : str, optional
        Path to ``eos.mr`` for the mass-radius curve; skipped if ``None``.
    formats : tuple
        Output formats, e.g. ``("pdf", "png")``.
    """
    os.makedirs(outdir, exist_ok=True)
    nb = data.nb
    report = check_consistency(data, verbose=False)
    report["nb"] = nb
    igaps = report["grid"].get("gap_indices", [])

    paths = []
    specs = [
        ("pressure", lambda ax: plot_pressure(data, ax, gap_indices=igaps)),
        ("sound_speed", lambda ax: plot_sound_speed(data, ax, report=report)),
        ("composition", lambda ax: plot_composition(data, ax)),
        ("gibbs_duhem", lambda ax: plot_gibbs_duhem(report, ax)),
    ]
    if mr_path is not None and os.path.exists(mr_path):
        specs.append(("mass_radius", lambda ax: plot_mass_radius(mr_path, ax)))

    for name, fn in specs:
        fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
        fn(ax)
        for fmt in formats:
            p = os.path.join(outdir, f"eos_{name}.{fmt}")
            fig.savefig(p, dpi=150)
            paths.append(p)
        plt.close(fig)

    # combined 2x2 overview panel
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.0), constrained_layout=True)
    plot_pressure(data, axes[0, 0], gap_indices=igaps)
    axes[0, 0].set_title("Pressure", fontsize=10)
    plot_sound_speed(data, axes[0, 1], report=report)
    axes[0, 1].set_title(r"Sound speed $c_s^2$", fontsize=10)
    plot_composition(data, axes[1, 0])
    axes[1, 0].set_title("Composition", fontsize=10)
    plot_gibbs_duhem(report, axes[1, 1])
    axes[1, 1].set_title("Gibbs-Duhem residual", fontsize=10)
    for fmt in formats:
        p = os.path.join(outdir, f"eos_overview.{fmt}")
        fig.savefig(p, dpi=150)
        paths.append(p)
    plt.close(fig)
    return paths


__all__ = ["make_figures", "plot_pressure", "plot_sound_speed",
           "plot_composition", "plot_gibbs_duhem", "plot_mass_radius"]
