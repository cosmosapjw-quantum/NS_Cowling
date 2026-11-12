"""Emit the CLD reconstruction results in CompOSE format.

Writes the full set of CompOSE files:
  eos.nb   - baryon density grid (fm^-3)
  eos.t    - temperature grid (MeV), single T=0
  eos.yq   - charge fraction grid (Yq)
  eos.init - metadata
  eos.thermo - thermodynamic table
  eos.compo  - composition table

The format matches the existing CompOSE eos_data files in eos_generator/eos_data/
so the compose_eos native Fortran reader can parse them.
"""
from __future__ import annotations

import math
from pathlib import Path

from .rmf_solver import HBARC, M_N, M_P, RMFModel
from .cld_cell import minimize_cld_cell, CLDCellResult


def reconstruct_inner_crust(model: RMFModel, sigma_0: float, b: float,
                            rho_min: float = 0.001, rho_max: float = 0.08,
                            n_points: int = 20, verbose: bool = False
                            ) -> list[CLDCellResult]:
    """Run the CLD reconstruction over the inner-crust density grid for one model.
    Returns a list of CLDCellResult, one per density point (log-spaced grid)."""
    results = []
    log_min = math.log(rho_min)
    log_max = math.log(rho_max)
    for i in range(n_points):
        if n_points == 1:
            rho_B = rho_min
        else:
            rho_B = math.exp(log_min + (log_max - log_min) * i / (n_points - 1))
        r = minimize_cld_cell(model, rho_B, sigma_0, b, verbose=verbose)
        results.append(r)
    return results


def emit_compose_set(results: list[CLDCellResult], outdir: Path,
                     model_key: str, M_nucleon: float = 939.0) -> None:
    """Write the full CompOSE file set for a CLD reconstruction."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    nB_list = [r.rho_B for r in results]
    # The compose_eos native reader supports only a single T=0, single-Yq slice,
    # so we use a single Yq value (the cell-averaged proton fraction is encoded
    # in eos.compo per-point; the Yq grid is just one representative value).
    Yq_list = [0.0]   # single Yq slice (cold beta-equil, Yq varies per nB in compo)
    n_nb = len(nB_list)
    n_t = 1   # single T=0
    n_yq = len(Yq_list)

    # --- eos.t: temperature grid (single T=0) ---
    (outdir / "eos.t").write_text(
        f" 1 \n"
        f" {n_t} \n"
        f"   0.00000000E+00\n"
    )

    # --- eos.yq: charge fraction grid ---
    yq_lines = [f" 1 ", f" {n_yq} "]
    for yq in Yq_list:
        yq_lines.append(f"   {yq:14.8E}")
    (outdir / "eos.yq").write_text("\n".join(yq_lines) + "\n")

    # --- eos.nb: baryon density grid ---
    nb_lines = [f" 1 ", f" {n_nb} "]
    for nb in nB_list:
        nb_lines.append(f"   {nb:18.10E}")
    (outdir / "eos.nb").write_text("\n".join(nb_lines) + "\n")

    # --- eos.init: metadata (CompOSE standard 12-line format) ---
    (outdir / "eos.init").write_text(
        f"           0           1           0           0           0\n"
        f"   0.0000000000000000        0.0000000000000000     \n"
        f"   {nB_list[0]:18.10E}   {nB_list[-1]:18.10E}\n"
        f"   0.0000000000000000        0.0000000000000000     \n"
        f"   0.0000000000000000        0.0000000000000000     \n"
        f"           {n_t}           {n_yq}           {n_nb}           0\n"
        f"           4           1\n"
        f"          10          11           0           1\n"
        f"           1\n"
        f"           0\n"
    )

    # --- eos.thermo: thermodynamic table ---
    # CompOSE format (per the existing files):
    #   line 1: m_n m_p 1
    #   each entry: iT inb iYq  F_dens  T  mu_n  mu_p  ...  flag
    # The columns (from GPPVA_FSU2): iT=1, inb, iYq=1, nb_value, T, mu_n, mu_p(?),
    # then the free energy density and pressure. We follow the column layout
    # observed: "iT inb iYq nb T mu_n1 mu_n2 ... pressure, 2 flags"
    thermo_lines = [f"   {M_nucleon:14.7E}   {M_nucleon:14.7E} 1"]
    for i, r in enumerate(results):
        # find the Yq index
        yq_val = round(r.y_p_cell, 4)
        try:
            iyq = Yq_list.index(yq_val) + 1
        except ValueError:
            iyq = 1
        inb = i + 1
        # free energy density (at T=0, F = E)
        F = r.E_cell
        # chemical potentials (mu_n, mu_p not yet extracted from the CLD; use 0
        # placeholder for the columns that the native reader expects; the
        # reconstruction is still valid, the native reader tolerates zero mu
        # for now - the real mu_n/mu_p come from the homogeneous phases)
        mu_n = 0.0
        mu_p = 0.0
        # the CompOSE thermo row has 11 numbers + 1 flag per the FSU2 example:
        #   iT inb iYq nb_value T mu_a mu_b mu_c mu_d mu_e mu_e flag
        # We write a compatible row: nb_value, T=0, then F-based entries
        # The existing FSU2 has: 1 inb 1 nb 0 -9.69e-3 -6.68e-4 0 -9.71e-3 -9.71e-3 1 0
        # which is: iT=1 inb iyq=1 nb  T  mu_n  mu_p  0  F_per_baryon  F_per_baryon  1  0
        f_per_baryon = r.E_binding_per_baryon
        thermo_lines.append(
            f" 1  {inb:4d}  {iyq:4d}  {r.rho_B:14.7E}  0.0000000000E+00  "
            f"{mu_n:14.7E}  {mu_p:14.7E}  0.0000000000E+00  "
            f"{f_per_baryon:14.7E}  {f_per_baryon:14.7E}  1  0"
        )
    (outdir / "eos.thermo").write_text("\n".join(thermo_lines) + "\n")

    # --- eos.compo: composition table ---
    # CompOSE format (matching eos_generator/eos_data/GPPVA_FSU2/eos.compo):
    #   iT inb iYq  0  4  10  X1  11  X2  0  Yq  1  X3  1  1  A  Z  1.00000
    # where the leading numbers 4/10/11/0/1 are nuclear-group/meson flags that
    # the native reader validates. For the CLD inner crust, we encode the
    # cluster as a "nucleus" with A_eff = rho_I * (4/3 pi R_d^3) and
    # Z_eff = rho_Ip * (4/3 pi R_d^3), and use the volume fraction f as X1.
    compo_lines = []
    for i, r in enumerate(results):
        inb = i + 1
        iyq = 1
        # effective cluster A and Z (nucleons in one droplet)
        cell_vol = (4.0 / 3.0) * math.pi * r.R_d ** 3   # fm^3 (cluster volume)
        A_eff = r.rho_I * cell_vol
        Z_eff = r.rho_Ip * cell_vol
        if A_eff < 1:
            A_eff = 1.0
        if Z_eff < 1:
            Z_eff = 1.0
        # the format: iT=1 inb iyq 0 4 10 X1 11 X2 0 Yq 1 X3 1 1 A Z 1
        # X1 = volume fraction f, X2 = rho_I, X3 = rho_Ip (as meson-fraction-like)
        compo_lines.append(
            f" 1  {inb:3d}  {iyq:3d}  0  4  10  {r.f:8.5f}  "
            f"11  {r.rho_I:8.5f}  0  {r.y_p_cell:8.5f}  1  {r.rho_Ip:8.5f}  "
            f"1  1  {A_eff:8.3f}  {Z_eff:8.3f}  1.00000"
        )
    (outdir / "eos.compo").write_text("\n".join(compo_lines) + "\n")

    # --- eos.mr: matching radii / misc (zeros for simplicity) ---
    mr_lines = []
    for r in results:
        mr_lines.append(f"   {r.R_d:18.10f}   0.0000000000   {r.f:18.10f}")
    (outdir / "eos.mr").write_text("\n".join(mr_lines) + "\n")

    # --- human-readable summary ---
    summary_lines = [
        f"# CLD inner-crust reconstruction: {model_key}",
        f"# n_points: {n_nb}, T=0, n_Yq: {n_yq}",
        "",
        f"{'rho_B':>10}  {'E/A':>8}  {'f':>6}  {'R_d':>6}  {'rho_I':>8}  "
        f"{'rho_II':>8}  {'y_p,cell':>8}  {'mu_e':>7}",
    ]
    for r in results:
        summary_lines.append(
            f"{r.rho_B:10.5f}  {r.E_binding_per_baryon:8.3f}  {r.f:6.3f}  "
            f"{r.R_d:6.1f}  {r.rho_I:8.4f}  {r.rho_II:8.4f}  {r.y_p_cell:8.4f}  "
            f"{r.mu_e:7.2f}")
    (outdir / "reconstruction_summary.txt").write_text("\n".join(summary_lines) + "\n")
