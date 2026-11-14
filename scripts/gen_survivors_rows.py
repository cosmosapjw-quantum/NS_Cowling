"""Generate the survivors_rows.tex table-body fragment for the RevTeX report.

Each row ends with `\\%` (commented newline) so that \\input does not inject
a stray newline that triggers a Misplaced \\noalign error inside tabular.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
stats = json.loads((ROOT / "download" / "matched" / "summary_stats.json").read_text())

lines = []
for i, s in enumerate(stats, 1):
    crust = s["crust"].replace("_", r"\_")
    core = s["core"].replace("_", r"\_")
    nb_max = f"{s['nb_max']:.3f}"
    p_max = f"{s['P_max']:.1f}"
    mu_max = f"{s['mu_b_max']:.0f}"
    cs2_mean = f"{s['cs2_mean']:.3f}"
    n = s["n"]
    readback = s["n"]
    lines.append(
        f"{i} & {crust} & {core} & {n} & {readback} & "
        f"{nb_max} & {p_max} & {mu_max} & {cs2_mean} \\%")

out = Path(ROOT / "download" / "report" / "survivors_rows.tex")
out.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"wrote {len(lines)} rows to {out}")
