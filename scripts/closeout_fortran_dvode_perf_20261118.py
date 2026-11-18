"""C10 delivery manifest: commit-scoped, byte-addressed executable evidence.

Must be run after final full regression + C9 audit. Does not claim that the
Python Cowling eigensolver has been migrated to Fortran.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data/fortran_dvode_perf_20261118"
DOC = ROOT / "docs/fortran_dvode_perf_20261118"
BUILD = Path("/home/cosmosapjw/ns_cowling_wt/scratch/fortran_dvode_perf_20261118/dvode/build")
ORIGINAL = Path("/home/cosmosapjw/ns_cowling_wt/scratch/fortran_dvode_perf_20261118/baseline/build")
BASE = "990d07f4f0383f3777f65ed8d3454263b09e14f2"
BRANCH = "research/fortran-dvode-perf-accuracy-20261118"


def sha(f):
    return hashlib.sha256(f.read_bytes()).hexdigest()


def data(f):
    return json.loads(f.read_text())


def main():
    head = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                   text=True).strip()
    branch = subprocess.check_output(["git", "-C", str(ROOT), "branch",
                                     "--show-current"], text=True).strip()
    if head != BASE or branch != BRANCH:
        raise ValueError(f"wrong immutable base or branch {head} {branch}")

    full = (D / "FULL_PYTEST.log").read_text()
    check = re.search(r"(\d+) passed, (\d+) skipped, (\d+) subtests passed", full)
    if not check or (D/"FULL_PYTEST.exit").read_text().strip() != "0":
        raise ValueError("full pytest verification not admitted")
    counts = list(map(int, check.groups()))
    if counts[0] < 1531 or counts[1] != 9 or counts[2] != 4:
        raise ValueError("regression count disagrees with frozen evidence")

    audited = data(D/"INDEPENDENT_AUDIT.json")
    if audited["status"] != "PASS_SCOPED" or audited["failed"]:
        raise ValueError("independent C9 audit has unresolved failures")
    if audited["count"] < 71:
        raise ValueError("C9 coverage below recorded contracts")

    batch = data(D/"BATCH_BENCHMARK.json")
    parallel = data(D/"PARALLEL_STRESS.json")
    performance = data(D/"BENCHMARK.json")
    if not batch["serial_vs_batch_all_32_states_exact"]:
        raise ValueError("native batch does not preserve per-star results")
    if parallel["status"] != "PASS" or parallel["number_distinct_stars"] != 6:
        raise ValueError("native independent-process stress not validated")
    if len(performance["cases"]) != 12:
        raise ValueError("performance evidence incomplete")
    if not (DOC/"C1_C10_RESULTS_KO.md").exists():
        raise ValueError("C1-C10 report missing")

    signed_sources=[
        ROOT/"fortran_eos/build.sh",
        ROOT/"fortran_eos/src/compose_eos_tov.f90",
        ROOT/"fortran_eos/app/compose_eos_fortran.f90",
        ROOT/"fortran_eos/vendor/odepack/PROVENANCE.md",
        ROOT/"python_frontend/compose_eos/native_tov_love.py",
    ]
    files=[
        *signed_sources,
        *list((ROOT/"fortran_eos/vendor/odepack").glob("*.f")),
        *list((ROOT/"fortran_eos/tests").glob("*20261118.f90")),
        *list((ROOT/"python_frontend/tests").glob("*20261118.py")),
        *list((ROOT/"scripts").glob("*dvode*20261118.py")),
        *list((ROOT/"scripts").glob("*native_fortran*20261118.py")),
        *list(D.glob("*.json")),
        *list(D.glob("*.log")),
        *list(DOC.glob("*.md")),
    ]
    files=sorted(set(f for f in files if f.is_file() and f not in (
        D/"DELIVERY_RECEIPT.json", D/"CLOSEOUT_RUN.log"
    )))
    if not files:
        raise ValueError("no deliverables")
    refs={str(f.relative_to(ROOT)): sha(f) for f in files}
    outcome=dict(
        schema="ns_cowling.fortran_dvode_performance_accuracy.v1",
        status="READY_TO_PUBLISH_SCOPED",
        base_commit=BASE,
        branch=BRANCH,
        actual_source_changed_files=[
            "fortran_eos/app/compose_eos_fortran.f90",
            "fortran_eos/build.sh",
            "fortran_eos/src/compose_eos_tov.f90",
        ],
        source_data_modified=False,
        implementation_verified={
            "Fortran_TOV_Love":True,
            "Fortran_Maxwell_TOV_Love":True,
            "Fortran_DVODE_Adams":True,
            "Fortran_DVODE_BDF":True,
            "Fortran_corrected_RK45":True,
            "Fortran_native_32star_batch":True,
            "Python_pure_IO_frontend":True,
            "Fortran_Cowling_f_p1":False,
            "previous_Python_Cowling_code_left_as_reference":True,
        },
        source_provenance={
            "PK1_native_eos_sha256":sha(ROOT/"data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose/eos.thermo"),
            "Netlib_VODE_original_sha256":sha(ROOT/"fortran_eos/vendor/odepack/vode.f"),
        },
        actual_runtime={
            "Fortran_compiler":"gfortran 13.3.0",
            "Fortran_build_flags":"-O2 -fcheck=all, with original Netlib and LINPACK",
            "baseline_binary_sha256":sha(ORIGINAL/"compose_eos_fortran"),
            "final_binary_sha256":sha(BUILD/"compose_eos_fortran"),
            "full_test_counts":dict(zip(("passed","skipped","subtests_passed"), counts)),
            "independent_audit_count":audited["count"],
        },
        measured_performance={
            "PK1_4097_speedup_vs_corrected_RK45":audited["measured_4097_speedup_vs_corrected_RK"],
            "PK1_4097_speedup_vs_invalid_old_physics_reference":audited["measured_4097_speedup_vs_wrong_old_RK"],
            "PK1_32star_batch_speedup":audited["measured_32_star_batch_speedup"],
            "PK1_6process_parallel_single_measurement_speedup":audited["measured_6_native_process_parallel_speedup"],
        },
        accepted_scope="Fortran TOV+Love on same native tabulated EOS only; benchmark low-density core truncation not a whole physical neutron star.",
        unresolved_next_action="COWLING-NATIVE-01 (Fortran f/p1 mode ODE + sharp slow/rapid conversion with physics-root and node gates)",
        artifacts_sha256=refs,
    )
    D.mkdir(parents=True,exist_ok=True)
    (D/"DELIVERY_RECEIPT.json").write_text(json.dumps(outcome,indent=2,ensure_ascii=False)+"\n")
    print("READY_TO_PUBLISH_SCOPED",len(files),"hashed files",audited["count"],
          "C9 checks","full pytest",counts)


if __name__ == "__main__":
    main()
