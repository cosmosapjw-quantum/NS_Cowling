#!/usr/bin/env python3
"""Read the native generator's CSV products and present the Fortran audit as JSON."""
import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PROFILE = (
    "nb_fm3", "epsilon_MeV_fm3", "pressure_MeV_fm3", "N", "Z", "E_minus_mn_MeV",
    "max_nucleon_density_fm3", "max_particle_error", "mu_L_MeV",
    "mu_L_finite_difference_error_MeV", "profile_sector", "dilation_P_rel_error",
)
UNIFORM = (
    "nb_fm3", "epsilon_MeV_fm3", "pressure_MeV_fm3", "Ye", "mu_n_MeV", "mu_p_MeV",
    "mu_e_MeV", "beta_residual_MeV", "beta_curvature_MeV_fm3", "first_law_P_rel_error",
    "cs2_over_c2_FD", "minimum_sampled_charge_curvature_MeV_fm3",
)


def send_columns(stream, path, names):
    with path.open(newline="") as file:
        reader = csv.DictReader(file)
        if any((reader.fieldnames or []).count(name) != 1 for name in names):
            raise ValueError(f"Missing or duplicate required column in {path}")
        records = []
        for record in reader:
            try:
                # Only transport numeric tokens. CSV order and quoting are immaterial.
                records.append(" ".join(repr(float(record[name])) for name in names))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid CSV field at {path}:{reader.line_num}") from exc
    print(len(records), file=stream)
    print(*records, sep="\n", file=stream)


def run(primary, refined, executable=None):
    primary, refined = Path(primary).resolve(), Path(refined).resolve()
    if primary.samefile(refined):
        raise ValueError("Use independent coarse and fine profile runs")
    if executable is None:
        executable = os.environ.get("NV_PROFILE_VERIFY_EXE")
    if executable is None:
        build = Path(os.environ.get("COMPOSE_EOS_NATIVE_BUILD_DIR", ROOT / "build/compose_eos_fortran")).resolve()
        subprocess.run(["bash", str(ROOT / "fortran_eos/build.sh"), "--ensure", str(build)],
                       stdout=subprocess.DEVNULL, check=True)
        executable = build / "nv_profile_verify"
    stream = io.StringIO()
    for mode in ("printed", "d03x100"):
        send_columns(stream, primary / f"{mode}.csv", PROFILE)
        send_columns(stream, refined / f"{mode}.csv", ("nb_fm3", "pressure_MeV_fm3", "E_minus_mn_MeV"))
        send_columns(stream, primary / f"{mode}_anchors.csv",
                     ("pressure_MeV_fm3", "E_minus_mn_MeV", "source_ground_E_minus_mn_MeV"))
    send_columns(stream, primary / "uniform_beta.csv", UNIFORM)
    process = subprocess.run([str(Path(executable).resolve()), str(primary), str(refined)], input=stream.getvalue(),
                             text=True, stdout=subprocess.PIPE, check=True)
    report = {}
    for line in process.stdout.splitlines():
        tokens = line.split()
        if tokens[0] == "mesh_metadata_present":
            report["radial_mesh_provenance"] = [
                "generator_namelist" if flag == "T" else "legacy_campaign_assumption" for flag in tokens[1:]
            ]
            continue
        group, key, *fields = tokens
        item = report if group == "root" else report.setdefault(group, {})
        if len(fields) > 1:
            item[key] = [float(v) for v in fields]
        elif group == "root" or key in ("rows", "clustered_rows", "negative_pressure_anchor_rows"):
            item[key] = int(fields[0])
        else:
            item[key] = float(fields[0])
    return report


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("primary", type=Path)
    cli.add_argument("refined", type=Path)
    cli.add_argument("--exe", type=Path, help="Use a prebuilt nv_profile_verify executable")
    args = cli.parse_args()
    print(json.dumps(run(args.primary, args.refined, args.exe), indent=2))
