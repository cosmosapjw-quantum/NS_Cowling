"""Read-only C9 adversarial acceptance of Fortran-native TOV/Love changes.

Facts verified at the actual bounded source/physical scope, not a blanket
native Fortran Cowling claim. Physics, numerics, source and runtime are separate.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import math
import re
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/fortran_dvode_perf_20261118"
SCRATCH = Path("/home/cosmosapjw/ns_cowling_wt/scratch/fortran_dvode_perf_20261118")
BIN = SCRATCH / "dvode/build/compose_eos_fortran"
PK1 = ROOT / "data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose/eos.thermo"
SRC = ROOT / "fortran_eos/vendor/odepack"

EXPECTED_VENDOR = {
    "vode.f": "76b3858bb1e32f416a29363dcf56e0889ac25aceb5911289f3f25fe275634e95",
    "dgefa.f": "764bd0a87dba77464eeac5c232a264b9de0aee192c8db152799b939579d95317",
    "dgesl.f": "2ef2289c6cb929825fff0fc9a9fd57c556c3cde97b21351cfdff140be645fbda",
    "dgbfa.f": "5aa49203557431ba5dd15719985c9d31915f936f2e34d5e30afd467326071873",
    "dgbsl.f": "f46ceaa18dd52ded3257b9e1a94a85786835ec9af21c7ca3f83687f2d4b13ddf",
}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p: Path) -> dict:
    return json.loads(p.read_text())


def main() -> None:
    checks = []
    failed = []

    def check(label: str, condition) -> None:
        value = bool(condition)
        checks.append(dict(name=label, result="PASS" if value else "FAIL"))
        if not value:
            failed.append(label)

    check("base_commit_exact", subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True
    ).strip() == "990d07f4f0383f3777f65ed8d3454263b09e14f2")
    check("branch_is_isolated", subprocess.check_output(
        ["git", "-C", str(ROOT), "branch", "--show-current"], text=True
    ).strip() == "research/fortran-dvode-perf-accuracy-20261118")

    input_json = read(SCRATCH/"baseline/INPUT.json")
    check("native_input_PΚ1_sha_locked", sha(PK1) == input_json["source_table_hash"])
    check("input_rows_4097", input_json["densities"] == [4097, 1025, 257])
    for name, fingerprint in EXPECTED_VENDOR.items():
        check("netlib_unchanged_" + name, sha(SRC/name) == fingerprint)

    benchmark = read(DATA/"BENCHMARK.json")
    check("twelve_condition_real_benchmarks", len(benchmark["cases"]) == 12)
    check("all_benchmarks_fortran_only", all(
        c["native_mode"] in ("rk45","dvode_adams","dvode_bdf")
        and c["source_table_sha256"] == sha(
            SCRATCH/"baseline"/f"PK1_native_{c['source_table_rows']}.csv"
        )
        for c in benchmark["cases"]
    ))
    check("all_7_repeated_times", all(len(c["measurements"]) == 7 for c in benchmark["cases"]))
    for c in benchmark["cases"]:
        label = c["backend_label"] + "_" + str(c["source_table_rows"])
        times = [x["seconds"] for x in c["measurements"]]
        check("median_genuine_" + label,
              abs(statistics.median(times)-c["wall_median_seconds"])<1e-13)
        x = c["result"]
        check("physical_state_finite_" + label,
              all(math.isfinite(x[k]) and x[k]>0
                  for k in ("M_msun","R_km","compactness","k2","Lambda"))
              and 0<x["compactness"]<.5)
        if not label.startswith("legacy_binary"):
            expected = (2/3)*x["k2"]/x["compactness"]**5
            check("lambda_love_identity_" + label,
                  abs(expected-x["Lambda"])/abs(x["Lambda"])<1e-12)

    subset = {x["backend_label"]: x for x in benchmark["cases"]
              if x["source_table_rows"] == 4097}
    a=subset["corrected_dvode_adams"]
    b=subset["corrected_dvode_bdf"]
    rk=subset["corrected_rk45"]
    old=subset["legacy_binary"]
    check("Adams_faster_same_model_than_corrected_RK",
          a["wall_median_seconds"] < rk["wall_median_seconds"])
    check("Adams_faster_than_BDF_on_nonstiff_case",
          a["wall_median_seconds"] < b["wall_median_seconds"])
    for field,limit in (("M_msun",3e-7),("R_km",2e-6),
                        ("k2",1e-5),("Lambda",1e-5)):
        check("matched_Adams_vs_BDF_"+field,
              abs(a["result"][field]-b["result"][field])/
              abs(b["result"][field]) < limit)
    check("old_science_is_different_and_not_accuracy_reference",
          old["result"]["k2"] < a["result"]["k2"]/15)
    check("minimum_real_perf_gain",
          rk["wall_median_seconds"] / a["wall_median_seconds"] > 1.5)
    check("original_legacy_speed_gain_is_not_accuracy_matched",
          old["wall_median_seconds"] / a["wall_median_seconds"] > 5)

    batch=read(DATA/"BATCH_BENCHMARK.json")
    check("batch_exact_32_stars",
          batch["serial_vs_batch_all_32_states_exact"]
          and batch["central_pressures"] == 32)
    check("batch_faster_than_32_cli_invocations",
          batch["measured_batch_speedup"] > 1.0)

    parallel=read(DATA/"PARALLEL_STRESS.json")
    check("DVODE_separate_process_only",
          parallel["status"] == "PASS"
          and parallel["max_parallel_native_processes"] == 6
          and parallel["num_threads_per_process"] == 1)
    check("parallel_has_six_valid_stars",
          len(parallel["results"]) == parallel["number_distinct_stars"] == 6)
    check("parallel_process_speedup_positive",
          parallel["speedup_serial_over_parallel"] > 1.0)

    for name,needle in [
        ("TIDAL_RHS_PHYSICS.log","PASS tidal RHS independent Hinderer Q identity"),
        ("test_tov_love_literature_20261118_dvode_adams.log","PASS analytic Love k2 normalization"),
        ("test_tov_love_literature_20261118_dvode_bdf.log","PASS analytic Love k2 normalization"),
        ("test_tov_love_literature_20261118_rk45.log","PASS analytic Love k2 normalization"),
        ("test_tov_love_dvode_adams.log","TOV-Love smooth + Maxwell branches: PASS"),
        ("test_tov_love_dvode_bdf.log","TOV-Love smooth + Maxwell branches: PASS"),
        ("test_tov_love_rk45.log","TOV-Love smooth + Maxwell branches: PASS"),
    ]:
        file=DATA/name
        check("native_test_" + name, file.exists() and needle in file.read_text())

    all_tests = DATA/"FULL_PYTEST.log"
    suite_text = all_tests.read_text() if all_tests.is_file() else ""
    match = re.search(r"(\d+) passed, (\d+) skipped, (\d+) subtests passed", suite_text)
    check("all_frontend_pipeline_regressions_pass",
          match is not None and int(match.group(1)) >= 1500 and
          int(match.group(2)) == 9 and int(match.group(3)) == 4)
    exit_path = DATA/"FULL_PYTEST.exit"
    check("full_pytest_actual_exit_zero",
          exit_path.is_file() and exit_path.read_text().strip() == "0")
    check("compiled_fortran_production_binary_available",
          BIN.is_file() and BIN.stat().st_size > 100_000)

    # The new Python module is strictly a frontend: import/parse/launch only.
    facade=(ROOT/"python_frontend/compose_eos/native_tov_love.py").read_text()
    check("frontend_without_scipy_integrator",
          "solve_ivp" not in facade and "odeint" not in facade and
          "from scipy" not in facade)
    check("native_Cowling_NOT_migrated_no_false_claim",
          (ROOT/"python_frontend/compose_eos/stellar_cowling.py").is_file()
          and not (ROOT/"fortran_eos/src/compose_eos_cowling_native.f90").exists())

    output=dict(
        status="PASS_SCOPED" if not failed else "FAILED",
        count=len(checks),
        failed=failed,
        checked=checks,
        source_base_commit="990d07f4f0383f3777f65ed8d3454263b09e14f2",
        physics_scope="native Fortran homogeneous smooth EOS and two-phase Maxwell TOV+Love",
        native_calculation_only=True,
        legacy_Python_Cowling_not_migrated=True,
        correct_Love_gr_identity_test="independent Fortran literature RHS + k2 formula",
        benchmark_source="PK1 native homogeneous uniform EOS; truncated at n=0.015 fm^-3, not full surface crust",
        maxwell_fixture_scope="analytic causal two-branch EOS, not proof of full-pasta physical coexistence",
        compiler="gfortran 13.3, -O2, -fcheck=all, Netlib DVODE and LINPACK",
        full_pytest_count=dict(
            passed=int(match.group(1)) if match else None,
            skipped=int(match.group(2)) if match else None,
            subtests=int(match.group(3)) if match else None,
        ),
        measured_4097_speedup_vs_corrected_RK=
            rk["wall_median_seconds"] / a["wall_median_seconds"],
        measured_4097_speedup_vs_wrong_old_RK=
            old["wall_median_seconds"] / a["wall_median_seconds"],
        measured_32_star_batch_speedup=batch["measured_batch_speedup"],
        measured_6_native_process_parallel_speedup=
            parallel["speedup_serial_over_parallel"],
    )
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA/"INDEPENDENT_AUDIT.json").write_text(json.dumps(output,indent=2)+"\n")
    print("SCOPED_NATIVE_DVODE_AUDIT",output["status"],"checks",
          output["count"],"failed",len(failed),failed[:12])
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
