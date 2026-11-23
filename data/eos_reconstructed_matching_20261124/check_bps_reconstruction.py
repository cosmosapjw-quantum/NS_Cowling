"""Check the generated historical BPS outer-crust rows and compare at reference densities.

Run from any directory. An optional previous bundle supplies an independent numerical
baseline. The default comparison only needs the generator output and local source files.
"""
import argparse
import bisect
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "source_data" / "bps"
MODEL = "BPS_MTW1964_MS1966"


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log_hermite(x, x0, x1, y0, y1, slope0, slope1):
    """Interpolate log(y) with supplied endpoint dlog(y)/dlog(x)."""
    interval = math.log(x1 / x0)
    t = math.log(x / x0) / interval
    h01 = t*t*(3-2*t)
    delta = (h01*math.log(y1/y0) +
             interval*(t*(1-t)**2*slope0 + t*t*(t-1)*slope1))
    return y0 * math.exp(delta)


def evaluate(states, grid, n):
    i = bisect.bisect_left(grid, n)
    if i < len(grid) and grid[i] == n:
        return "pure", states[i]["pressure"], states[i]["epsilon"]
    if i == 0 or i == len(grid):
        return "outside_support", None, None
    a, b = states[i-1], states[i]
    na, nb = grid[i-1], grid[i]
    if a["phase"] != b["phase"]:
        if a["pressure"] != b["pressure"]:
            raise ValueError("non-isobaric coexistence pair")
        return ("density_jump", a["pressure"],
                a["epsilon"] + (n-na)/(nb-na)*(b["epsilon"]-a["epsilon"]))
    p = log_hermite(n, na, nb, a["pressure"], b["pressure"], a["Gamma"], b["Gamma"])
    eps = log_hermite(n, na, nb, a["epsilon"], b["epsilon"],
                      a["mu_b"]*na/a["epsilon"], b["mu_b"]*nb/b["epsilon"])
    return "pure", p, eps


def compose_rows(max_pressure):
    source = SOURCE / "compose68"
    nraw = (source / "eos.nb").read_text().split()
    first, last = map(int, nraw[:2])
    thermo = [s.split() for s in (source / "eos.thermo").read_text().splitlines()]
    mn = float(thermo[0][0])
    if len(nraw)-2 != last-first+1 or len(thermo)-1 != last-first+1:
        raise ValueError("CompOSE axis length inconsistent")
    out = []
    for idx, (n, raw) in enumerate(zip(map(float, nraw[2:]), thermo[1:]), first):
        if int(raw[1]) != idx:
            raise ValueError("CompOSE index mismatch")
        p = n*float(raw[3])
        if p <= max_pressure*(1+1e-6):
            out.append((idx, n, p, n*mn*(1+float(raw[9]))))
    return out


def check_prior(directory, states, generated):
    expected = read_rows(directory / "eos.csv")
    mt = read_rows(directory / "masses.csv")
    fg = read_rows(generated / "myers_mass_catalog.csv")
    et = read_rows(directory / "transitions.csv")
    ft = read_rows(generated / "myers_bps_transitions.csv")
    if (len(mt), len(fg), len(ft), len(et), len(states), len(expected)) != (130,130,11,11,623,623):
        raise ValueError("mass, EOS or jump row count mismatch")
    for r, f in [*zip(expected, states), *zip(mt, fg)]:
        if (int(r["A"]),int(r["Z"])) != (int(f["A"]),int(f["Z"])):
            raise ValueError("nuclide mismatch")
    rel = {}
    for key, ref in (("n_b","nB_fm3"),("pressure","P_MeV_fm3"),("epsilon","eps_MeV_fm3"),("mu_b","g_MeV")):
        rel[key] = max(abs(f[key]/float(r[ref])-1) for f, r in zip(states, expected))
    rel["transition_pressure"] = max(abs(float(f["pressure"])/float(r["P_MeV_fm3"])-1) for f, r in zip(ft,et))
    mass_abs = max(abs(float(m["atomic_mass_excess_MeV"])-float(f["atomic_mass_excess_MeV"])) for m,f in zip(mt,fg))
    if max(rel.values()) > 1e-9 or mass_abs > 1e-8:
        raise ValueError("Fortran/Python numerical cross-check exceeded tolerance")
    return dict(implementation="independent previous Python/SciPy reconstruction",
                reference_sha256={p.name:sha(p) for p in sorted(directory.glob("*.csv"))},
                mass_rows=130,states=623,phase_count=12,transition_pairs=11,
                nuclides_identical_rowwise=True, maximum_relative_difference=rel,
                max_atomic_mass_excess_abs_difference_MeV=mass_abs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior-model", type=Path, help="previous historical_1964_1966 model directory")
    parser.add_argument("--native", type=Path, default=ROOT / "generated_sources" / MODEL / "generated")
    parser.add_argument("--out", type=Path, default=ROOT / "bps_validation")
    args = parser.parse_args()
    states = [{k:float(v) for k,v in r.items()} for r in read_rows(args.native / "myers_bps_states.csv")]
    grid = [s["n_b"] for s in states]
    if not all(b>a for a,b in zip(grid,grid[1:])):
        raise ValueError("non-monotone baryon density")
    species = [(int(s["phase"]),int(s["A"]),int(s["Z"])) for s in states]
    layers = [s for i,s in enumerate(species) if i==0 or s != species[i-1]]
    if [p for p,a,z in layers]!=list(range(201,213)):
        raise ValueError("phase chain mismatch")
    table5 = [(int(s["row"]), float(s["nb_fm3"]), float(s["P_MeV_fm3"]), float(s["eps_MeV_fm3"]))
              for s in read_rows(SOURCE / "Table5_BPS_outer.csv")]
    refs = {"BPS1971_Table5": table5, "CompOSE_68_APR": compose_rows(max(p for _,_,p,_ in table5))}
    if [len(r) for r in refs.values()] != [43,49]:
        raise ValueError("reference density selection inconsistent")
    rows, summary = [], {}
    for name, ref in refs.items():
        current = []
        for i,n,p,eps in ref:
            tag, pred_p, pred_eps = evaluate(states, grid, n)
            row = dict(reference=name, row=i, n_b_fm3=n, state=tag,
                       P_ref_MeV_fm3=p, epsilon_ref_MeV_fm3=eps,
                       P_rec_MeV_fm3=pred_p, epsilon_rec_MeV_fm3=pred_eps,
                       dP_over_P=None, deps_over_eps=None)
            if pred_p is not None:
                row["dP_over_P"] = pred_p/p-1
                row["deps_over_eps"] = pred_eps/eps-1
                current.append(row)
            rows.append(row)
        if len(current)!=43 or sum(r["state"]=="density_jump" for r in current)!=3:
            raise ValueError("reference/model support or coexistence count changed")
        residual = [r["dP_over_P"] for r in current]
        worst = max(current,key=lambda r:abs(r["dP_over_P"]))
        summary[name] = dict(reference_rows=len(ref), overlapping=len(current), density_jump_rows=3,
                             median_abs_delta_P_over_P=statistics.median(map(abs,residual)),
                             max_abs_delta_P_over_P=abs(worst["dP_over_P"]),
                             worst_row=worst["row"], worst_density_fm3=worst["n_b_fm3"],
                             worst_signed_delta_P_over_P=worst["dP_over_P"])
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out / "reference_residuals.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,tuple(rows[0]));writer.writeheader();writer.writerows(rows)
    report = dict(model=MODEL, table_status="reconstructed", units="MeV and fm", phases_A_Z=layers,
                  reference_sha256={str(p.relative_to(SOURCE)):sha(p) for p in sorted(SOURCE.rglob("*")) if p.is_file() and p.suffix!=".md"},
                  same_density_references=summary, low_n=grid[0], drip_n=grid[-1],
                  drip_pressure_MeV_fm3=states[-1]["pressure"])
    (args.out / "reference_check.json").write_text(json.dumps(report,indent=2)+"\n")
    if args.prior_model is not None:
        previous=check_prior(args.prior_model, states, args.native)
        (args.out / "previous_independent_implementation.json").write_text(json.dumps(previous,indent=2)+"\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
