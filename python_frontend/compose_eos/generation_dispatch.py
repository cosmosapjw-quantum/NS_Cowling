"""Source-checkout dispatch for independently evaluated, explicitly scoped EOSs.

Success means the requested rows and native output receipts exist. It does not
promote a candidate to an original author calibration or a complete crust.
"""
from __future__ import annotations
import hashlib
import gzip
import json
import math
import os
from pathlib import Path
import subprocess
import sys

CORE_IDS = ("APR_APR", "BBB_BHF2", "GM_GM1", "RG_KDE0v1", "RG_SLY4")
SPG_IDS = tuple(f"SPG_M{i}" for i in range(1, 6))
TF_IDS = ("GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW", "GPPVA_FSU2",
          "GPPVA_FSU2R", "GPPVA_FSU2H", "GPPVA_NL3wrL55", "GPPVA_TM1e")
from .primary_model_contracts import PRIMARY_IDS, primary_crust_catalog
MODEL_IDS = CORE_IDS + SPG_IDS + TF_IDS + PRIMARY_IDS
from .native_model_contracts import SPG_INFERRED_VARIANT, SPG_INFERRED_SURFACE

SPG_VARIANTS = (SPG_INFERRED_VARIANT, "AME2408_candidate", "geometric_area__nine",
                "geometric_area__eight_measured", "fixed_symmetric_area__nine",
                "fixed_symmetric_area__eight_measured")
SCRIPTS = dict(core="generate_core_eos.py", spg="reoptimize_spg_surfaces.py",
               tf="build_gppva_tf.py")
NATIVE_FILES = ("eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo")


def generator_catalog():
    from .native_model_contracts import MODEL_CONTRACTS
    result = {}
    for model, contract in MODEL_CONTRACTS.items():
        result[model] = dict(family="core", region="core_only",
                            composition="nucleon and lepton fractions",
                            **contract)
    for model in SPG_IDS:
        result[model] = {
            "family": "spg", "region": "inner_crust_only",
            "representation": "spherical GDFM CLD with explicit candidate surface fit",
            "calibration": "candidate; original M1-M5 surface coefficients unresolved",
            "composition": "free nucleons, leptons and sharp-interface cluster",
            "variants": list(SPG_VARIANTS),
            "conditional_inferred_input": {"variant": SPG_INFERRED_VARIANT, **SPG_INFERRED_SURFACE},
        }
    for model in TF_IDS:
        result[model] = {
            "family": "tf", "region": "sampled_inner_crust_cells_only",
            "representation": "spatial RMF Thomas-Fermi geometry/mesh campaign",
            "composition": "total nucleon profiles and lepton fractions; native free/cluster masks unavailable for nonuniform cells",
            "geometry_scope": "finite geometry, seed, mesh and radius search; no global minimum or continuous interpolation claim",
        }
    result.update(primary_crust_catalog())
    return result


def repository_root(path=None):
    root = Path(path).expanduser().resolve() if path else Path(__file__).resolve().parents[2]
    required = [root / "scripts" / script for script in SCRIPTS.values()]
    required += [root / "python_frontend/compose_eos/api.py"]
    if not (root / ".git").exists() or not all(p.is_file() for p in required):
        raise ValueError("generation requires the NS_Cowling source checkout; use --repository-root")
    return root


def _finite_positive(values, label, integer=False):
    if not values or any(not math.isfinite(v) or v <= 0 or
                         (integer and (not isinstance(v, int) or isinstance(v, bool)))
                         for v in values):
        raise ValueError(f"{label} must contain finite positive values")
    if any(right <= left for left, right in zip(values, values[1:])):
        raise ValueError(f"{label} must be strictly increasing")


def generation_request(model, *, points=None, n_min=None, n_max=None,
                       variant="AME2408_candidate", densities=None,
                       resolutions=None, workers=1, include_lower_crust=False):
    catalog = generator_catalog()
    if model not in catalog:
        raise ValueError(f"unknown model: {model}")
    if not isinstance(workers, int) or isinstance(workers, bool) or workers < 1:
        raise ValueError("workers must be a positive integer")
    family = catalog[model]["family"]
    if include_lower_crust and family == "core":
        raise ValueError("lower crust segments are only available for crust models")
    request = {"model": model, "family": family, "contract": catalog[model]}
    if family == "core":
        if densities is not None or resolutions is not None or variant != "AME2408_candidate":
            raise ValueError("core models use --points/--n-min/--n-max")
        count = 161 if points is None else points
        if not isinstance(count, int) or isinstance(count, bool) or count < 3:
            raise ValueError("core points must be an integer >=3")
        lo, hi = catalog[model]["default_interval_fm3"]
        lo = lo if n_min is None else n_min
        hi = hi if n_max is None else n_max
        _finite_positive([lo, hi], "core density range")
        request.update(points=count, n_min=lo, n_max=hi)
    elif family == "primary_crust":
        if include_lower_crust or any(x is not None for x in (n_min, n_max, densities, resolutions)) or variant != "AME2408_candidate":
            raise ValueError("primary crust models use their original source support; only --points is configurable")
        count = catalog[model]["default_points"] if points is None else points
        if not isinstance(count, int) or isinstance(count, bool) or not catalog[model].get("minimum_points",3) <= count <= 100000:
            raise ValueError("primary crust points must be an integer in [3,100000]")
        if model == "BBP1971_table5" and count != 26:
            raise ValueError("the raw primary BBP table contains exactly26 knots; use pressure_integral for resampling")
        request.update(points=count, primary_source_path=catalog[model]["source_path"])
    elif family == "spg":
        if any(x is not None for x in (n_min, n_max, densities, resolutions)):
            raise ValueError("SPG uses its independently solved drip/endpoint grid; use --points/--variant")
        count = 129 if points is None else points
        if not isinstance(count, int) or isinstance(count, bool) or count < 17 or (count - 1) % 8:
            raise ValueError("SPG points must be 8k+1 and at least17")
        if variant not in SPG_VARIANTS:
            raise ValueError("unknown SPG surface variant")
        request.update(points=count, variant=variant, workers=workers)
    else:
        if any(x is not None for x in (points, n_min, n_max)) or variant != "AME2408_candidate":
            raise ValueError("TF models require --densities and --resolutions, without core/SPG grid options")
        _finite_positive(densities, "TF densities")
        _finite_positive(resolutions, "TF resolutions", integer=True)
        if len(densities) < 3 or min(resolutions) < 8:
            raise ValueError("TF native output requires >=3 densities and resolutions >=8")
        request.update(densities=list(densities), resolutions=list(resolutions), workers=workers)
    if include_lower_crust:
        request["include_lower_crust"] = True
    return request


def _command(root, output, request):
    from .native_generation import write_request
    request_file=write_request(root,output,request)
    return [str(root / "fortran_eos/run-generator.sh"), str(request_file), str(output)]


def _read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _native_present(path):
    missing = [name for name in NATIVE_FILES
               if not (path / name).is_file() or (path / name).stat().st_size == 0]
    if missing:
        raise ValueError("missing native output files: " + ", ".join(missing))


def _verify_hashes(directory, hashes, required):
    if not isinstance(hashes, dict) or not set(required).issubset(hashes):
        raise ValueError("missing native output checksums")
    for relative, expected in hashes.items():
        path = (directory / relative).resolve()
        if not path.is_relative_to(directory.resolve()) or _sha(path) != expected:
            raise ValueError("native output checksum mismatch: " + relative)


def verify_generation(output, request):
    """Check child evidence; never interpret process exit zero as completeness."""
    model = request["model"]
    directory = output / model
    if request["family"] == "core":
        summary = _read_json(output / "SUMMARY.json")
        record = _read_json(directory / "RECEIPT.json")
        rows = record["rows"]
        if summary.get("status") != "COMPLETE" or set(summary["models"]) != {model}:
            raise ValueError("core campaign is incomplete or model identity differs")
        if record.get("native_readback") != "PASS" or summary["models"][model].get("native_readback") != "PASS":
            raise ValueError("core native readback was not verified")
        if summary["models"][model]["points"] != len(rows):
            raise ValueError("core row count differs from receipt")
        expected = [request["n_min"] * (request["n_max"]/request["n_min"]) ** (i/(request["points"]-1))
                    for i in range(request["points"])]
        actual = [row["n_b"] for row in rows]
        _finite_positive(actual, "generated core density grid")
        if len(actual) < len(expected) or any(
                not any(math.isclose(a, n, rel_tol=5e-13) for a in actual) for n in expected):
            raise ValueError("requested core density rows are missing")
        # APR deliberately inserts the two coexistence endpoints.
        allowed = expected[:]
        if model == "APR_APR":
            allowed += [record["phase_context"][key]["n_b"] for key in ("left", "right")]
        if len(actual) > len(expected) + (2 if model == "APR_APR" else 0) or any(
                not any(math.isclose(a, n, rel_tol=5e-13) for n in allowed) for a in actual):
            raise ValueError("unexpected generated core density rows")
        _native_present(directory)
        _verify_hashes(directory, record["output_sha256"], NATIVE_FILES)
        return {"rows": len(rows), "native_directory": str(directory)}
    if request["family"] == "spg":
        summary = _read_json(output / "SUMMARY.json")
        directory /= request["variant"]
        record = _read_json(directory / "RECEIPT.json")
        case = summary.get("cases", [])
        if (summary.get("status") != "COMPLETED_SCOPED" or
            summary.get("source_unchanged_during_execution") is not True or
            len(case) != 1 or case[0].get("model") != model or
            case[0].get("label") != request["variant"]):
            raise ValueError("SPG campaign identity, completeness or source stability failed")
        if (record.get("status") != "REOPTIMIZED_NATIVE_COMPOSITION_VERIFIED" or
            record.get("model") != model or record.get("label") != request["variant"] or
            record.get("rows") != request["points"] or
            record.get("stationary_rows") != request["points"] or record.get("failures")):
            raise ValueError("SPG requested stationary rows are incomplete")
        exports = record.get("native_exports", [])
        if (len(exports) != 1 or exports[0].get("path") != "compose" or
            exports[0].get("status") != "NATIVE_READBACK_VERIFIED" or
            exports[0].get("rows") != request["points"]):
            raise ValueError("SPG complete native export is missing")
        _native_present(directory / "compose")
        _verify_hashes(directory, record["output_sha256"], ["compose/" + name for name in NATIVE_FILES])
        return {"rows": record["rows"], "native_directory": str(directory / "compose")}
    campaign = _read_json(output / "CAMPAIGN.json")
    record = _read_json(directory / "RECEIPT.json")
    cases = campaign.get("models", [])
    if (campaign.get("source_unchanged_during_execution") is not True or len(cases) != 1 or
        cases[0].get("model") != model or cases[0].get("rows") != len(request["densities"]) or
        campaign.get("densities_fm3") != request["densities"] or
        campaign.get("resolutions") != request["resolutions"] or record.get("model_id") != model):
        raise ValueError("TF campaign identity or requested rows are incomplete")
    stages = record["stages"]
    final = stages[-1]
    if ([stage.get("points") for stage in stages] != request["resolutions"] or
        final.get("points") != request["resolutions"][-1] or
        final.get("requested_points") != len(request["densities"]) or
        final.get("selected_stationary_points") != len(request["densities"]) or
        record["export"].get("status") != "EXPORTED_NATIVE_READBACK_VERIFIED" or
        record["export"].get("row_count") != len(request["densities"])):
        raise ValueError("TF all-row native export is not verified")
    _native_present(directory / "compose")
    if not (directory / "composition.json").is_file() or not (directory / "selected_states.json.gz").is_file():
        raise ValueError("TF companion composition/profile output is missing")
    required_hashes = ["compose/" + name for name in NATIVE_FILES]
    required_hashes += ["composition.json", "selected_states.json.gz"]
    if (directory / "INTERVAL_SUPPORT.json").exists():
        required_hashes.append("INTERVAL_SUPPORT.json")
    _verify_hashes(directory, record["export"].get("output_sha256"), required_hashes)
    with gzip.open(directory / "selected_states.json.gz", "rt") as stream:
        selected = json.load(stream)
    composition = _read_json(directory / "composition.json")
    if ([row["n_b"] for row in selected] != request["densities"] or
        composition.get("row_count") != len(request["densities"]) or
        [row["n_b"] for row in composition["rows"]] != request["densities"]):
        raise ValueError("TF profile/composition density rows differ from request")
    from .tf_selection import STATIONARY, feasibility
    for row, density in zip(selected, request["densities"]):
        try:
            finite_thermodynamics = all(math.isfinite(float(row[key])) for key in
                ("epsilon", "energy_per_baryon", "pressure", "mu_B", "mu_Q_lepton"))
        except (KeyError, TypeError, ValueError, OverflowError):
            finite_thermodynamics = False
        if (row.get("model_id") != model or row.get("points") != request["resolutions"][-1] or
            row.get("numerical_status") not in STATIONARY or
            row.get("muon_threshold_crossed") is not False or
            feasibility(row, density) != "FEASIBLE" or not finite_thermodynamics):
            raise ValueError(f"selected TF row lacks requested model/mesh or stationary feasible thermodynamic evidence at n={density}")
    return {"rows": len(request["densities"]), "native_directory": str(directory / "compose"),
            "composition": "sparse native masks; total species in companion profiles"}


def run_generation(model, output, *, root=None, **options):
    root = repository_root(root)
    request = generation_request(model, **options)
    output = Path(output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=False)
    child_output = output / "generated"
    command = _command(root, child_output, request)
    from .native_generation import source_files, surface_input, verify_native_generation
    tracked = [Path(__file__).resolve(), Path(__file__).with_name("native_generation.py"), *source_files(root)]
    if request["family"] == "primary_crust":
        tracked += [root / request["primary_source_path"], Path(__file__).with_name("primary_model_contracts.py")]
    if request["family"] == "spg":
        tracked.append(root / surface_input(root, model, request["variant"])["path"])
    initial = {str(p): _sha(p) for p in tracked}
    receipt = {"status": "RUNNING", "request": request, "command": command,
               "repository_root": str(root), "source_sha256": initial,
               "backend": "fortran",
               "completion_scope": "requested finite native rows/readback; source equivalence and full-crust completion remain separate"}
    if request["family"] == "primary_crust":
        source = root / request["primary_source_path"]
        receipt["primary_source_input"] = {"path": str(source), "sha256": _sha(source),
            "source_equivalence": request["contract"]["source_equivalence"]}
    if request["family"] == "spg":
        receipt["surface_input"] = surface_input(root, model, request["variant"])
    if request.get("include_lower_crust"):
        mass = root / "data/outer_crust/hfb22_ame2016_cu.csv"
        receipt["mass_input"] = {"path": str(mass.relative_to(root)), "sha256": _sha(mass)}
    receipt["request_sha256"] = _sha(output / "request.nml")
    path = output / "DISPATCH.json"
    save = lambda: path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    save()
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(root / "python_frontend") + os.pathsep + environment.get("PYTHONPATH", "")
    try:
        with (output / "execution.log").open("w") as log:
            child = subprocess.run(command, cwd=root, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT, check=False)
        receipt["child_returncode"] = child.returncode
        if child.returncode:
            receipt["status"] = "CHILD_FAILED"
            code = child.returncode if child.returncode > 0 else 128 - child.returncode
        else:
            if initial != {str(p): _sha(p) for p in tracked}:
                raise ValueError("dispatcher or runner source changed during execution")
            receipt["evidence"] = verify_native_generation(child_output, request)
            receipt["status"] = "REQUESTED_ROWS_NATIVE_VERIFIED"
            code = 0
    except Exception as exc:
        receipt.update(status="EVIDENCE_OR_EXECUTION_FAILED",
                       error_type=type(exc).__name__, error=str(exc))
        code = 1
    receipt["returncode"] = code
    save()
    print(json.dumps({"status": receipt["status"], "receipt": str(path), "returncode": code}))
    return code
