#!/usr/bin/env python3
"""Run the Fortran direct-EOS provider admission audit over a normalized EOS directory."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile


FIELDS = (
    "status", "row_count", "interval_count", "admitted_intervals", "rejected_intervals",
    "mu_strict_increasing", "pressure_strict_increasing", "g_strict_increasing",
    "raw_first_law_closed", "fully_admitted", "raw_energy_closure_max_rel",
    "raw_trapezoid_first_law_max_rel", "min_dg_dlogp", "reconstructed_n_max_rel",
)


def parse_records(output: str) -> dict[str, str]:
    result = {}
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) == 3 and fields[0] not in {"name", "protocol_version", "workspace", "case_id", "request_token"}:
            result[fields[0]] = fields[2]
    return result


def run_case(driver: Path, case: str, table: Path, work: Path) -> dict[str, str]:
    request = work / "request.nml"
    request.write_text(
        "&request request_token='direct-eos-v3-audit' /\n"
        "&parameters\n"
        f" table_path='{table.resolve()}', query_kind='pressure', query_value=1.0e30, side=0\n"
        "/\n",
        encoding="utf-8",
    )
    try:
        completed = subprocess.run(
            [str(driver), case, str(request)],
            check=False, capture_output=True, text=True, timeout=60,
        )
    finally:
        request.unlink(missing_ok=True)
    if completed.returncode:
        raise RuntimeError(
            f"Fortran case {case} failed for {table}: rc={completed.returncode}\n"
            f"stdout={completed.stdout}\nstderr={completed.stderr}"
        )
    return parse_records(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, required=True, help="compiled Fortran case_driver")
    parser.add_argument("--source-dir", type=Path, required=True, help="directory containing normalized source EOS CSVs")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bad-fixture", type=Path, help="optional deliberately inconsistent EOS fixture")
    args = parser.parse_args()
    driver, source_dir, output_dir = args.driver.resolve(), args.source_dir.resolve(), args.output_dir.resolve()
    if not driver.is_file():
        raise FileNotFoundError(driver)
    tables = sorted(p for p in source_dir.glob("*.csv") if p.name != "spectrum_summary_234.csv")
    if not tables:
        raise FileNotFoundError(f"No EOS CSV tables in {source_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows, failures = [], []
    with tempfile.TemporaryDirectory(prefix="direct-eos-v3-") as td:
        work = Path(td)
        for index, table in enumerate(tables, 1):
            values = run_case(driver, "table_diagnostics", table, work)
            row = {"file": table.name, **{key: values.get(key, "") for key in FIELDS}}
            rows.append(row)
            valid = (
                values.get("status") == "ok"
                and values.get("raw_first_law_closed") == "true"
                and values.get("fully_admitted") == "true"
                and values.get("rejected_intervals") == "0"
            )
            if not valid:
                failures.append({"file": table.name, "values": row})
            if index % 50 == 0:
                print(f"audited {index}/{len(tables)}", flush=True)
        bad_fixture_result = None
        if args.bad_fixture:
            values = run_case(driver, "table_diagnostics", args.bad_fixture.resolve(), work)
            bad_fixture_result = {key: values.get(key, "") for key in FIELDS}
    csv_path = output_dir / "fortran_direct_eos_cohort.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["file", *FIELDS])
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema": "fortran-direct-eos-cohort-audit-v2",
        "provider": "monotone_PCHIP_epsilon_minus_pressure_of_logP_with_firstlaw_density_v1",
        "n_tables": len(rows),
        "n_admitted": sum(row["fully_admitted"] == "true" for row in rows),
        "n_raw_first_law_closed": sum(row["raw_first_law_closed"] == "true" for row in rows),
        "total_intervals": sum(int(row["interval_count"]) for row in rows),
        "total_rejected_intervals": sum(int(row["rejected_intervals"]) for row in rows),
        "max_raw_knot_closure": max(float(row["raw_energy_closure_max_rel"]) for row in rows),
        "max_reconstructed_n_relative_change": max(float(row["reconstructed_n_max_rel"]) for row in rows),
        "max_trapezoid_first_law_diagnostic": max(float(row["raw_trapezoid_first_law_max_rel"]) for row in rows),
        "failures": failures,
        "bad_fixture": bad_fixture_result,
    }
    report["all_admissible"] = len(failures) == 0 and report["n_admitted"] == len(rows)
    json_path = output_dir / "fortran_direct_eos_cohort.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_admissible": report["all_admissible"], "n_tables": len(rows),
                      "n_admitted": report["n_admitted"], "n_failed": len(failures),
                      "max_raw_knot_closure": report["max_raw_knot_closure"],
                      "max_reconstructed_n_relative_change": report["max_reconstructed_n_relative_change"],
                      "bad_fixture": bad_fixture_result}, indent=2))
    return 0 if report["all_admissible"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
