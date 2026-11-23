"""Plot source-level diagnostics for the reconstructed crust campaign."""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/eos_reconstructed_matching_20261124"
PROFILES = ROOT / "data/eos_remaining_generated_20261124/nv_profiles"
NV = ROOT / "data/eos_original_papers_20261123/generated_sources/NV1973_author_fit/generated"


def read(path):
    return np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("output", type=Path)
    args = cli.parse_args()
    anchors = read(PROFILES / "printed_anchors.csv")
    d03 = read(PROFILES / "d03x100_anchors.csv")
    orbital = read(SOURCE / "orbital_hf/d03_triplet_native/SOURCE_VALUES.csv")
    reference = anchors["source_ground_E_minus_mn_MeV"]
    n = anchors["nb_fm3"]
    assert np.allclose(n, d03["nb_fm3"]) and np.allclose(n, orbital["n_b_fm3"])
    bps = read(SOURCE / "generated_sources/BPS_MTW1964_MS1966/generated/myers_bps_states.csv")
    inner = read(NV / "nv1973_states.csv")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "savefig.facecolor": "white"})
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.9, 3.9), layout="constrained")
    left.set_title("a   Energies at the 11 published densities", loc="left", fontweight="bold", fontsize=11)
    left.axhline(0, color="0.32", lw=.9, ls=":")
    for data, name, c, mk in [
        (anchors["E_minus_mn_MeV"], "Printed coefficients · constrained profiles", "#b14825", "o"),
        (d03["E_minus_mn_MeV"], "Adjusted D₀₃ · constrained profiles", "#007f79", "s"),
        (orbital["energy_per_baryon_offset_mev"], "Adjusted D₀₃ · orbital charge triplets", "#4c55a4", "^"),
    ]:
        left.plot(n, data-reference, color=c, marker=mk, ms=4.3, lw=1.1, label=name)
    left.set_xscale("log")
    left.set_xlabel(r"Baryon density $n_B$ [fm$^{-3}$]")
    left.set_ylabel(r"$E/A-(E/A)_{\rm NV73,Table\ 3}$ [MeV]")
    left.legend(fontsize=8, frameon=False, loc="best")
    left.grid(which="major", alpha=.15)
    right.set_title("b   Chemical potentials near neutron drip", loc="left", fontweight="bold", fontsize=11)
    p0, p1 = inner["pressure_mev_fm3"][0], bps["pressure"].max()
    p = np.linspace(p0, p1, 160)
    bmu = np.interp(p, bps["pressure"], bps["mu_b"])
    nmu = np.interp(p, inner["pressure_mev_fm3"], inner["mu_b_mev"])
    mn = bps["mu_b"][-1]
    right.plot(p*1e4, bmu-mn, color="#3F5C8A", lw=1.8, label="Reconstructed BPS outer crust")
    right.plot(p*1e4, nmu-mn, color="#af491e", lw=1.8, label="NV73 inner crust · published fit")
    right.fill_between(p*1e4, bmu-mn, nmu-mn, color="#9BA9AE", alpha=.14)
    right.set_xlim(p0*1e4, p1*1e4)
    right.set_xlabel(r"Pressure [$10^{-4}$ MeV fm$^{-3}$]")
    right.set_ylabel(r"$\mu_B-m_nc^2$ [MeV]")
    right.legend(fontsize=8, frameon=False, loc="upper left")
    right.grid(alpha=.15)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
