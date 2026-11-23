"""Plot native fixed-mass TOV/Love results for the public full-crust cohort."""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RESULTS = Path(__file__).resolve().parents[1] / "data/eos_native_upgrade_20261124/layered_stars"
METHODS = {
    "first_law_constrained_hermite_pressure_bridge": ("First-law Hermite", "o"),
    "cubic_F": ("Cubic F", "s"),
    "quintic_F": ("Quintic F", "^"),
    "cubic_P": ("Cubic P", "D"),
    "quintic_P": ("Quintic P", "v"),
}
# Each status is supplied by the native solver and endpoint audit.
STATUSES = {
    "no_failed_interval_in_table": ("Whole table unexcluded", "#12685F"),
    "below_first_failed_interval": ("First failure above $P_c$", "#316AA0"),
    "above_first_failed_interval": ("Failed interval traversed", "#AD441F"),
    "no_numerical_solution": ("No sampled crossing", "#666A70"),
}
GRID = "#DDE0DF"
INK = "#263631"


def load_results(path: Path):
    grouped = defaultdict(dict)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["request"] != "M1.4" or float(row["max_dlogp"]) != 0.01:
                continue
            group = grouped[row["core_id"]]
            if row["method"] in group:
                raise ValueError(f"duplicate method for {row['core_id']}")
            group[row["method"]] = row
    if not grouped:
        raise ValueError("no M1.4 results at max_dlogp 0.01")
    for core, group in grouped.items():
        if group.keys() != METHODS.keys():
            raise ValueError(f"incomplete method group for {core}")
        labels = {row["causal_interval_status"] for row in group.values()}
        if len(labels) != 1 or not labels.issubset(STATUSES):
            raise ValueError(f"unsupported or mixed method status in {core}: {labels}")
        for row in group.values():
            expected = ("no sampled ascending mass crossing in supplied interval"
                        if row["causal_interval_status"] == "no_numerical_solution" else "ok")
            if row["status"] != expected:
                raise ValueError(f"unexpected solver outcome in {core}: {row['status']}")
    return grouped


def draw(path: Path, output: Path):
    grouped = load_results(path)
    status_order = list(STATUSES)
    representative = {core: next(iter(group.values())) for core, group in grouped.items()}
    def sort_key(core):
        row = representative[core]
        return status_order.index(row["causal_interval_status"]), float(row["R_km"] or "inf"), core
    cores = sorted(grouped, key=sort_key)
    counts = Counter(row["causal_interval_status"] for group in grouped.values() for row in group.values())
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
                         "axes.labelcolor": INK, "xtick.color": INK, "pdf.fonttype": 42,
                         "axes.unicode_minus": False, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(11.0, 10.15))
    layout = fig.add_gridspec(1, 4, width_ratios=(2.5, 3.45, 3.45, 1.34),
                              left=.038, right=.976, top=.807, bottom=.114, wspace=.095)
    labels, radius, tidal, spans = [fig.add_subplot(layout[0,i]) for i in range(4)]
    for ax in (labels, radius, tidal, spans):
        ax.set_ylim(len(cores) - .48, -.52)
        ax.set_yticks([])
        ax.spines[["left","right","top","bottom"]].set_visible(False)
    for ax in (labels, spans):
        ax.set_xlim(0, 1)
        ax.set_xticks([])
    radius.set_xlim(10.98, 14.78)
    radius.set_xticks([11,12,13,14])
    tidal.set_xlim(125, 1320)
    tidal.set_xticks([200,400,600,800,1000,1200])
    for ax in (radius, tidal):
        ax.grid(axis="x", color=GRID, linewidth=.55)
        ax.tick_params(axis="x", labelsize=8, length=3, color=GRID)
        ax.set_axisbelow(True)
    labels.set_title("Core table", loc="left", fontsize=10, weight="bold", pad=10)
    radius.set_title(r"Cutoff radius  $R_{\rm cut}$  (km)", loc="left", fontsize=10, weight="bold", pad=10)
    tidal.set_title(r"Tidal deformability  $\Lambda_{\rm cut}$", loc="left", fontsize=10, weight="bold", pad=10)
    spans.set_title("Method spans", loc="left", fontsize=10, weight="bold", pad=10)
    spans.text(.16, 1.002, r"$\Delta R$ (m)", transform=spans.transAxes, ha="center", fontsize=7)
    spans.text(.80, 1.002, r"$\Delta\Lambda$", transform=spans.transAxes, ha="center", fontsize=7)
    previous = None
    for index, core in enumerate(cores):
        row = representative[core]
        status = row["causal_interval_status"]
        color = STATUSES[status][1]
        if previous is not None and previous != status:
            for ax in (labels, radius, tidal, spans):
                ax.axhline(index-.5, lw=.9, color="#A1ACA9")
        previous = status
        name = row["core_model"].replace("_", " ")
        display_name = {
            "2020-MVCD-D1MS-full": "MVCD D1MS",
            "DDMEX ORIGINAL TABLE RHO0": "DDMEX original density",
            "DDMEX PROVIDER ROUNDED N0": "DDMEX rounded density",
            "NL3 PROVIDER PRECISION": "NL3 provider precision",
            "PK1 CANONICAL": "PK1 canonical",
            "TM1 CANONICAL": "TM1 canonical",
        }.get(name, name)
        labels.text(.003, index, core.removeprefix("core_"), va="center", color=color, fontsize=8.4, weight="bold")
        labels.text(.165, index, display_name, va="center", fontsize=8.5)
        if status == "no_numerical_solution":
            for ax in (radius, tidal):
                ax.text(.5, index, "no sampled crossing", ha="center", va="center",
                        transform=ax.get_yaxis_transform(), fontsize=8, color=color)
            spans.text(.48, index, "—", va="center", ha="center", color=color)
            continue
        for method_index, (method, (_, marker)) in enumerate(METHODS.items()):
            value = grouped[core][method]
            y = index + .123*(method_index-2)
            for ax, column in ((radius,"R_km"), (tidal,"Lambda")):
                ax.plot(float(value[column]), y, ls="none", marker=marker, ms=3.2,
                        markerfacecolor=color, markeredgecolor="white", markeredgewidth=.28)
        rvalues = [float(x["R_km"]) for x in grouped[core].values()]
        lvalues = [float(x["Lambda"]) for x in grouped[core].values()]
        for xpos, delta in ((.15, 1000*(max(rvalues)-min(rvalues))), (.79, max(lvalues)-min(lvalues))):
            spans.text(xpos, index, f"{delta:.3f}", ha="center", va="center", fontsize=8.2)
    fig.text(.038, .962, r"Radius and tides at $1.4\,M_\odot$", fontsize=23, weight="bold", va="top")
    fig.text(.039, .916, f"Reconstructed BPS–NV crust  ·  {len(cores)} public core tables  ·  five crust–core matching methods",
             fontsize=11)
    for xpos, (status, (name, color)) in zip((.039, .307, .545, .795), STATUSES.items()):
        fig.text(xpos, .878, "●", color=color, fontsize=10, va="center")
        fig.text(xpos+.015, .878, f"{name} ({counts[status]})", va="center", fontsize=9.1)
    handles = [Line2D([], [], color=INK, ls="none", marker=marker, markersize=4, label=label)
               for label, marker in METHODS.values()]
    fig.legend(handles=handles, ncol=5, loc="center left", bbox_to_anchor=(.035,.850),
               frameon=False, columnspacing=2.1, handlelength=1.1, borderaxespad=0, fontsize=8.4)
    fig.text(.039, .065, "Rows group the five matching methods. The last two columns give their maximum minus minimum. Core IDs distinguish source versions.",
             fontsize=8.6)
    fig.text(.039, .045, r"Cutoff at the lowest tabulated pressure ($5.63\times10^{-15}$ MeV fm$^{-3}$). "
             r"Ascending mass crossing searched over $P_c=5$–$300$ MeV fm$^{-3}$ within table support.", fontsize=8.6)
    fig.text(.039, .025, "Color tracks the nominal cold-EOS causal endpoint test. The native integration uses a maximum logarithmic pressure step of 0.01.",
             fontsize=8.6)
    metadata = {"Title": "Cutoff radius and tidal deformability at 1.4 solar masses", "CreationDate": None} if output.suffix.lower() == ".pdf" else {}
    fig.savefig(output, dpi=180, metadata=metadata)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=RESULTS / "public_stars.csv")
    parser.add_argument("--output", type=Path, default=RESULTS / "public_fixed_mass.pdf")
    args = parser.parse_args()
    draw(args.input, args.output)
