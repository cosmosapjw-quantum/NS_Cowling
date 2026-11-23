#!/usr/bin/env python3
"""Collect native NV HF receipts, reconstruct the fixed-A triplet table and check its support."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

MN = 939.56542052
MP = 938.27208816
ENERGY = "energy_per_baryon_MeV"
PRESSURE = "pressure_MeV_fm3"
RAW = ["nb_fm3", ENERGY, PRESSURE, "mu_n_kin_MeV", "mu_p_kin_MeV", "mu_e_total_MeV"]
CONFIG = ["nb_fm3", "neutrons", "protons", "mesh", "seed_kind", "coefficient_mode",
          "outer_smoothing", "lmax", "radial_states", "radius_fm"]
FIELDS = [ENERGY, PRESSURE, *RAW[3:], "mu_n_addition_kin_MeV", "mu_p_addition_kin_MeV",
          "beta_residual_MeV", "beta_interval_lower_MeV", "beta_interval_upper_MeV",
          "selfconsistent_residual", "dilation_pressure_error_MeV_fm3", "particle_error",
          "minimum_effective_mass_coefficient"]
CAMPAIGNS = {"site96": "PLAN.json", "literal_alt": "LITERAL_ALT_PLAN.json",
             "charge_ext": "CHARGE_EXT_PLAN.json", "followup": "FOLLOWUP_PLAN.json"}


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_cells(root, output):
    rows, receipts = [], {}
    plans = output / "plans"
    plans.mkdir(exist_ok=True)
    for campaign, planname in CAMPAIGNS.items():
        folder = root / campaign
        result = json.loads((folder / "SUMMARY.json").read_text())
        plan = (root / planname).read_bytes()
        assert result["execution"] == "COMPLETED" and result["input_identities_unchanged"]
        ids = {j["id"] for j in json.loads(plan)["jobs"]}
        records = [json.loads(line) for line in (folder / "RESULTS.jsonl").read_text().splitlines()]
        assert ids == {r["id"] for r in records} and len(records) == len(ids)
        hashes = {Path(k).name: v for k, v in result["identities_after"].items()}
        assert hashlib.sha256(plan).hexdigest() == hashes[planname]
        (plans / planname).write_bytes(plan)
        receipts[campaign] = {"jobs": len(ids), "statuses": result["statuses"], "input_sha256": hashes}
        for rec in records:
            f = rec["native"]["fields"]
            valid = rec["status"] == "FIXED_COMPOSITION_CELL_CONVERGED"
            assert valid == (f["status"] == "ok" and f["converged"] == "T"
                             and f["energy_valid"] == "T" and f["spectral_complete"] == "T")
            row = {"site": int(rec["id"][1:3]), "campaign": campaign, "id": rec["id"],
                   "purpose": rec["purpose"], "native_status": f["status"],
                   **{k: f[k] for k in CONFIG},
                   "iterations": f["iterations"], "max_iterations": rec["config"]["max_iterations"],
                   "mixing": rec["config"]["mixing"],
                   "maximum_raw_neutron_density_fm3": f["maximum_raw_neutron_density_fm3"],
                   "maximum_raw_proton_density_fm3": f["maximum_raw_proton_density_fm3"],
                   **{k: f[k] if valid else "" for k in FIELDS},
                   "native_output_sha256": rec["output_sha256"]}
            rows.append(row)
    rows.sort(key=lambda r: (r["site"], int(r["coefficient_mode"]), r["campaign"], int(r["protons"]), r["id"]))
    write_rows(output / "all_native_cells.csv", rows)
    return rows, receipts


def collect_triplets(rows, output):
    triplets = [r for r in rows if r["campaign"] == "site96" or r["purpose"] == "literal_fixed_A_triplets"]
    assert len(triplets) == 66
    write_rows(output / "native_triplets_both_modes.csv", triplets)
    winners = []
    charge = []
    for site in range(1, 12):
        cells = [r for r in triplets if r["site"] == site and r["coefficient_mode"] == "1"]
        assert len(cells) == 3 and all(r["native_status"] == "ok" for r in cells)
        winner = min(cells, key=lambda r: float(r[ENERGY]))
        winners.append(winner)
        tested = {int(r["protons"]): r for r in rows if r["site"] == site and r["coefficient_mode"] == "1"
                  and r["seed_kind"] == "0" and r["mesh"] == "96" and r["native_status"] == "ok"}
        z = min(tested, key=lambda k: float(tested[k][ENERGY]))
        extended = tested[z]
        lo = float(extended["beta_interval_lower_MeV"])
        hi = float(extended["beta_interval_upper_MeV"])
        row = {"site": site, "A": int(winner["neutrons"])+int(winner["protons"]), "nb_fm3": winner["nb_fm3"],
               "triplet_Z": winner["protons"], "extended_Z": z, "tested_Z": " ".join(map(str, sorted(tested))),
               "triplet_energy_offset_MeV": winner[ENERGY], "extended_energy_offset_MeV": extended[ENERGY],
               "extended_beta_interval_lower_MeV": lo, "extended_beta_interval_upper_MeV": hi,
               "extended_beta_interval_contains_zero": lo <= 0 <= hi}
        for side, label in [(-1, "lower"), (1, "upper")]:
            neighbor = tested.get(z + side)
            row[f"adjacent_{label}_minus_extended_energy_MeV_per_cell"] = (
                (float(neighbor[ENERGY])-float(extended[ENERGY]))*row["A"] if neighbor else "")
        charge.append(row)
    write_rows(output / "fixed_A_charge_minima.csv", charge)
    return winners


def check_projection(winners, normalized, output):
    with normalized.open() as stream:
        eos = list(csv.DictReader(stream))
    assert len(eos) == len(winners) == 11
    points, segments = [], []
    differences = []
    for i, (r, f) in enumerate(zip(eos, winners), 1):
        n, p, e, mu = (float(r[k]) for k in ["n_b_fm3", "pressure_mev_fm3", "energy_mev_fm3", "mu_b_mev"])
        assert n == float(f["nb_fm3"]) and p > 0
        assert math.isclose(e, n*(MN+float(f[ENERGY])), rel_tol=2e-15)
        assert math.isclose(p, float(f[PRESSURE]), rel_tol=2e-15)
        assert abs(mu-MN-float(f["mu_n_kin_MeV"])) < 3e-12
        assert abs(MN*float(r["q5_mu_l_over_mn"])-float(f["beta_residual_MeV"])) < 3e-12
        g = (e+p)/n
        points.append((n,p,e,g))
        differences.append(mu-g)
    for i, ((n,p,e,g),(n2,p2,e2,g2)) in enumerate(zip(points, points[1:]), 1):
        secant = (e2-e)/(n2-n)
        slope = (p2-p)/(e2-e)
        segments.append(dict(left_site=i, right_site=i+1, left_nb_fm3=n, right_nb_fm3=n2,
                             energy_secant_MeV=secant, left_euler_g_MeV=g, right_euler_g_MeV=g2,
                             pressure_energy_secant=slope, euler_bracketed=g <= secant <= g2,
                             increasing_pressure=p2>p, increasing_energy=e2>e, subluminal_secant=0<slope<1))
    write_rows(output/"thermodynamic_secants.csv", segments)
    return {"density_fm3": [points[0][0],points[-1][0]], "pressure_energy_secant_range": [
                min(r["pressure_energy_secant"] for r in segments), max(r["pressure_energy_secant"] for r in segments)],
            "all_secants_euler_bracketed": all(r["euler_bracketed"] for r in segments),
            "all_secants_subluminal": all(r["subluminal_secant"] for r in segments),
            "all_energy_secants_increasing": all(a["energy_secant_MeV"] < b["energy_secant_MeV"]
                                                 for a,b in zip(segments,segments[1:])),
            "occupied_mu_n_minus_euler_g_range_MeV": [min(differences),max(differences)]}


def check_partial_literal(rows, output):
    branch = [r for r in rows if r["campaign"] == "literal_alt" and r["mesh"] == "96"
              and r["seed_kind"] == "1" and r["native_status"] == "ok"]
    assert [r["site"] for r in branch] == [1, 2, 3, 10]
    write_rows(output/"literal_partial_cells.csv", branch)
    intervals = []
    for a, b in zip(branch, branch[1:]):
        if b["site"] != a["site"] + 1:
            continue
        n1, n2 = float(a["nb_fm3"]), float(b["nb_fm3"])
        p1, p2 = float(a[PRESSURE]), float(b[PRESSURE])
        e1, e2 = n1*(MN+float(a[ENERGY])), n2*(MN+float(b[ENERGY]))
        secant = (e2-e1)/(n2-n1)
        g1, g2 = (e1+p1)/n1, (e2+p2)/n2
        intervals.append(dict(left_site=a["site"], right_site=b["site"], left_nb_fm3=n1,
                              right_nb_fm3=n2, energy_secant_MeV=secant, left_euler_g_MeV=g1,
                              right_euler_g_MeV=g2, pressure_energy_secant=(p2-p1)/(e2-e1),
                              euler_bracketed=g1 <= secant <= g2))
    write_rows(output/"literal_partial_secants.csv", intervals)
    return {"converged_uniform_seed_mesh96_sites": [r["site"] for r in branch],
            "bracketed_adjacent_intervals": sum(r["euler_bracketed"] for r in intervals),
            "adjacent_intervals_tested": len(intervals), "separated_site": 10,
            "native_cold_matching_support": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--importer", required=True, type=Path)
    args = parser.parse_args()
    out = args.output
    out.mkdir(exist_ok=True, parents=True)
    rows, provenance = read_cells(args.campaign_root, out)
    winners = collect_triplets(rows, out)
    raw = "# n_fm^-3 E_per_baryon_minus_mn_MeV P_MeV_fm^-3 occ_mu_n_minus_mn_MeV occ_mu_p_minus_mp_MeV mu_e_total_MeV\n"
    raw += "".join(" ".join(f[k] for k in RAW)+"\n" for f in winners)
    source = out / "d03_triplet_source.dat"
    source.write_text(raw)
    native = out / "d03_triplet_native"
    subprocess.run([str(args.importer.resolve()), "nv_fixed_cell_frontiers", str(source), str(native),
                    str(MN), str(MP)], check=True)
    if "beta_equilibrated=F" not in (native/"eos_data.nml").read_text():
        raise ValueError("Lost fixed-composition metadata")
    thermo = check_projection(winners, native/"normalized_eos.csv", out)
    grid = []
    for winner in winners:
        for x in rows:
            if x["purpose"] != "triplet_winner_seed_mesh_check" or x["site"] != winner["site"]:
                continue
            grid.append(dict(site=x["site"], Z=winner["protons"], mesh=x["mesh"], seed=x["seed_kind"],
                             status=x["native_status"], delta_E_per_A_MeV=float(x[ENERGY])-float(winner[ENERGY])
                             if x[ENERGY] else "", delta_pressure_fraction=float(x[PRESSURE])/float(winner[PRESSURE])-1
                             if x[PRESSURE] else ""))
    write_rows(out/"triplet_winner_refinement.csv", grid)
    counts = Counter(r["native_status"] for r in rows)
    summary = {"status": "CONDITIONAL_FIXED_A_ORBITAL_DIAGNOSTIC", "certified_beta_equilibrium": False,
               "literal_orbital_beta_EOS_available": False, "rows": len(rows), "native_statuses": counts,
               "selected_triplet_native_rows": [r["id"] for r in winners],
               "triplet_sites_with_zero_in_frontier_beta_interval": [r["site"] for r in winners if
                  float(r["beta_interval_lower_MeV"]) <= 0 <= float(r["beta_interval_upper_MeV"])],
               "thermodynamic_secants": thermo, "partial_literal_uniform_seed": check_partial_literal(rows, out),
               "native_importer_sha256":
                   hashlib.sha256(args.importer.read_bytes()).hexdigest(), "campaigns": provenance}
    (out / "orbital_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
