#!/usr/bin/env python3
"""Package completed A barotropes, portable consumers and actual provenance.

Only model RESULT records with completed readback and continuous-consumer checks
enter the accepted catalog.  A primary bundle contains one scientific baseline;
the full bundle additionally contains its declared construction sensitivities.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope

ACCEPTED = "CONSTRUCTED_CONTINUOUS_AND_READBACK_VERIFIED"
COMPLETE = "COMPLETE_STEP_A_BASELINE_AND_DECLARED_SENSITIVITIES"

EVALUATOR = '''#!/usr/bin/env python3
"""Evaluate the included piecewise EOS; NumPy and SciPy are required."""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"code/python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--list",action="store_true")
p.add_argument("--model",help="model ID from CATALOG.json; default is the primary baseline")
p.add_argument("--pressure",type=float,nargs="+",help="pressure(s) in MeV fm^-3, including zero")
p.add_argument("--transition-side",choices=("low","high"),default="high")
p.add_argument("--output",type=Path,help="new JSON file; existing files are never overwritten")
a=p.parse_args()
catalog=json.loads((ROOT/"CATALOG.json").read_text())
if a.list:
 print(json.dumps(catalog,indent=2));sys.exit(0)
if a.pressure is None:p.error("provide --pressure or --list")
selected=a.model or catalog["primary_model_id"]
byid={m["id"]:m for m in catalog["models"]}
if selected not in byid:p.error("model is not included in this bundle")
record=byid[selected]
contract=json.loads((ROOT/record["representation_path"]).read_text())
eos=WholeCrustBarotrope.from_contract(contract)
try:values=eos.state_at_pressure(a.pressure,transition_side=a.transition_side)
except ValueError as exc:p.error(str(exc))
result={"model":selected,"family":record["family"],"transition_side":a.transition_side,
 "units":{"n":"fm^-3","p":"MeV fm^-3","epsilon":"MeV fm^-3","mu":"MeV",
          "cs2_over_c2":"dimensionless","enthalpy":"dimensionless"},
 "states":[{k:float(values[k][i]) for k in values} for i in range(len(a.pressure))],
 "physical_support_n_fm3":record["physical_support_n_fm3"],
 "consumer_source":str(Path(sys.modules["compose_eos.whole_crust"].__file__).resolve())}
text=json.dumps(result,indent=2,allow_nan=False)+"\\n"
if a.output:
 with a.output.open("x") as f:f.write(text)
else:print(text,end="")
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, record):
    Path(path).write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")


def safe_relative(path):
    p = Path(path)
    if p.is_absolute() or ".." in p.parts:
        raise ValueError("repository-relative artifact path required")
    return p


def read_production(path):
    campaign = json.loads(path.read_text())
    if campaign.get("status") != COMPLETE:
        raise ValueError("packaging requires a completed A production campaign")
    models = []
    seen = set()
    for entry in campaign["models"]:
        name = entry["id"]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or name in seen:
            raise ValueError("unsafe or duplicate model ID")
        seen.add(name)
        directory = path.parent/"models"/name
        result = json.loads((directory/"RESULT.json").read_text())
        if result.get("status") != ACCEPTED:
            continue
        if (result.get("id") != name or result.get("readback", {}).get("status") != "PASS"
                or result.get("consumer_readback_validation", {}).get("status") != "PASS"):
            raise ValueError("accepted model lacks completed native readback or consumer validation")
        for field in ("method", "variant", "primary", "row_count"):
            if result[field] != entry[field]:
                raise ValueError("campaign/model metadata mismatch: "+field)
        for name_, digest in result["output_sha256"].items():
            if sha(directory/"compose"/safe_relative(name_)) != digest:
                raise ValueError("accepted CoMPOSE output hash mismatch")
        contract = json.loads((directory/"REPRESENTATION.json").read_text())
        consumer = WholeCrustBarotrope.from_contract(contract)
        zero = consumer.state_at_pressure(0.)
        if any(zero[k] != 0 for k in ("n", "p", "epsilon")):
            raise ValueError("representation lacks the advertised exact zero-pressure surface")
        provenance = json.loads((directory/"ROW_PROVENANCE.json").read_text())
        if len(provenance["provenance"]) != result["row_count"]:
            raise ValueError("row provenance is incomplete")
        models.append((result, consumer))
    if not models or sum(bool(r["primary"]) for r, _ in models) != 1:
        raise ValueError("exactly one accepted primary baseline is required")
    if len(models) != campaign.get("model_count"):
        raise ValueError("completed campaign's table count differs from accepted model records")
    # The executable-independent consumer is precisely the implementation used
    # by this production; an updated source must not be silently substituted.
    for name, digest in campaign.get("implementation_sha256", {}).items():
        if sha(ROOT/safe_relative(name)) != digest:
            raise ValueError("production implementation identity changed: "+name)
    return campaign, models


def catalog_for(campaign_path, campaign, models, *, package_kind):
    result = dict(schema="ns-cowling.whole-crust-deliverable.v1", package_kind=package_kind,
                  campaign_sha256=sha(campaign_path), scientific_baseline_count=1,
                  production_table_count=campaign["model_count"], packaged_table_count=len(models),
                  primary_model_id=next(r["id"] for r, _ in models if r["primary"]),
                  family=campaign["family"], next_step=campaign["next_step"], models=[])
    for record, _ in models:
        item = {k: record[k] for k in (
            "id", "family", "method", "variant", "primary", "row_count", "status",
            "physical_support_n_fm3", "compose_positive_density_range_fm3",
            "pressure_range_MeV_fm3", "upper_native_source_row", "lower_bridge_cs2_extrema",
            "surface_anchor_enthalpy", "composition_scope", "limitations")}
        prefix = "models/"+record["id"]
        item.update(compose_directory=prefix+"/compose", representation_path=prefix+"/REPRESENTATION.json",
                    result_path=prefix+"/RESULT.json", provenance_path=prefix+"/ROW_PROVENANCE.json",
                    surface_gamma=record["surface"]["gamma"])
        result["models"].append(item)
    return result


def catalog_csv(catalog):
    stream = io.StringIO()
    fields = ("id", "method", "variant", "primary", "row_count", "surface_gamma",
              "upper_native_source_row", "compose_directory", "representation_path")
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows({k: r[k] for k in fields} for r in catalog["models"])
    return stream.getvalue()


def continuous_curve(consumer, p_min, p_max, count=900):
    pressure = np.unique(np.r_[np.geomspace(p_min, p_max, count),
                               consumer.surface.p0, consumer.outer.p[-1], consumer.bridge.p_right])
    pressure = pressure[(pressure >= p_min) & (pressure <= p_max)]
    values = consumer.state_at_pressure(pressure)
    # A plot break preserves a declared density jump: no fictitious curve is
    # drawn through its forbidden density interval.
    for index in sorted(consumer.sharp_interfaces, reverse=True):
        pt = float(consumer.upper_p[index])
        if not p_min < pt < p_max:
            continue
        keep = abs(pressure-pt) > 16*np.finfo(float).eps*abs(pt)
        pressure = pressure[keep]
        values = {k: np.asarray(v)[keep] for k, v in values.items()}
        position = int(np.searchsorted(pressure, pt))
        low = consumer.state_at_pressure(pt, transition_side="low")
        high = consumer.state_at_pressure(pt, transition_side="high")
        pressure = np.insert(pressure, position, [pt, np.nan, pt])
        values = {k: np.insert(v, position, [low[k], np.nan, high[k]]) for k, v in values.items()}
    return values


def make_figure(models, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    selected = {(r["method"], r["variant"]): (r, c) for r, c in models}
    curves = []
    colors = {"baseline": "#146c94", "lower_endpoint_036": "#d07124", "lower_endpoint_037": "#598c35"}
    for method, style in (("quintic_P", "-"), ("css", "--")):
        record, eos = selected[(method, "baseline")]
        p_min = eos.surface.p0*1e-6**eos.surface.gamma
        data = continuous_curve(eos, p_min, eos.pressure_domain[1])
        axes[0, 0].loglog(data["n"], data["p"], style, lw=1.7,
                         label="Smooth baseline" if method == "quintic_P" else "CSS sensitivity")
        curves.append(dict(model_id=record["id"], panel="whole_pressure_density", values={
            k: [None if not np.isfinite(x) else float(x) for x in v] for k, v in data.items()}))
    axes[0, 0].set(title="Whole-crust and finite core support", xlabel=r"$n_B$ [fm$^{-3}$]",
                   ylabel=r"$P$ [MeV fm$^{-3}$]")
    axes[0, 0].text(.03, .06, "Exact P = 0 is supplied by the analytic surface", transform=axes[0, 0].transAxes, fontsize=9)
    for variant in ("baseline", "lower_endpoint_036", "lower_endpoint_037"):
        record, eos = selected[("quintic_P", variant)]
        label = "inner row "+str(record["upper_native_source_row"])
        data = continuous_curve(eos, eos.outer.p[-1]*.65, eos.bridge.p_right*1.6, count=400)
        axes[0, 1].plot(data["n"], data["p"], color=colors[variant], lw=1.7, label=label)
        data_cs = continuous_curve(eos, eos.outer.p[-1]*.5, eos.pressure_domain[1], count=1100)
        axes[1, 0].loglog(data_cs["n"], data_cs["cs2_over_c2"], color=colors[variant], lw=1.5, label=label)
        curves.append(dict(model_id=record["id"], panel="lower_bridge", values={
            k: [None if not np.isfinite(x) else float(x) for x in v] for k, v in data.items()}))
        curves.append(dict(model_id=record["id"], panel="sound_speed", values={
            k: [None if not np.isfinite(x) else float(x) for x in v] for k, v in data_cs.items()}))
    css_record, css_eos = selected[("css", "baseline")]
    css_data = continuous_curve(css_eos, css_eos.outer.p[-1]*.5, css_eos.pressure_domain[1], count=1100)
    axes[1, 0].loglog(css_data["n"], css_data["cs2_over_c2"], "--", color="#7b4895", lw=1.4,
                      label="CSS sensitivity")
    curves.append(dict(model_id=css_record["id"], panel="sound_speed", values={
        k: [None if not np.isfinite(x) else float(x) for x in v] for k, v in css_data.items()}))
    axes[0, 1].ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
    axes[0, 1].set(title="Declared lower-bridge sensitivity", xlabel=r"$n_B$ [fm$^{-3}$]",
                   ylabel=r"$P$ [MeV fm$^{-3}$]")
    axes[1, 0].axhline(1., color="#555555", linestyle=":", lw=1, label="causal bound")
    axes[1, 0].set(title="Continuous one-sided sound speed", xlabel=r"$n_B$ [fm$^{-3}$]",
                   ylabel=r"$c_s^2/c^2$")
    cutoffs = []
    ratio = np.logspace(-12, -2, 81)
    for variant, color in (("baseline", "#146c94"), ("surface_gamma_4over3", "#9b3c75")):
        record, eos = selected[("quintic_P", variant)]
        omitted_h = eos.surface.enthalpy(ratio*eos.surface.n0)
        axes[1, 1].loglog(ratio, omitted_h, color=color, lw=1.8,
                         label=r"$\gamma=5/3$" if variant == "baseline" else r"$\gamma=4/3$")
        cutoffs.append(dict(model_id=record["id"], ratio=ratio.tolist(), omitted_dimensionless_enthalpy=omitted_h.tolist()))
    axes[1, 1].set(title="Positive-density export cutoff convergence", xlabel=r"$n_{cut}/n_0$",
                   ylabel=r"Omitted $h=\int_0^{P_{cut}}dP/(\epsilon+P)$")
    for ax in axes.flat:
        ax.grid(True, which="major", alpha=.2)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle("A. Whole-crust connection — one baseline and construction sensitivities", fontsize=14)
    fig.savefig(output/"WHOLE_CRUST_OVERVIEW.png", dpi=180)
    fig.savefig(output/"WHOLE_CRUST_OVERVIEW.pdf")
    plt.close(fig)
    save(output/"PLOT_DATA.json", dict(schema="ns-cowling.whole-crust-plot-data.v1",
         evaluator="WholeCrustBarotrope.from_contract; no interpolation of exported plots",
         null_semantics="break at a declared density discontinuity", curves=curves, surface_cutoffs=cutoffs))


def copy_file(source, staging, target):
    source = Path(source)
    if source.is_symlink():
        raise ValueError("symbolic-link input is not an archived source file")
    destination = staging/safe_relative(target)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def bundle(campaign_path, campaign, models, output, *, primary_only):
    chosen = [(r, c) for r, c in models if r["primary"] or not primary_only]
    ids = {r["id"] for r, _ in chosen}
    catalog = catalog_for(campaign_path, campaign, chosen,
                          package_kind="primary_only" if primary_only else "full_campaign")
    archive_name = "whole_crust_primary.zip" if primary_only else "whole_crust_all_variants.zip"
    with tempfile.TemporaryDirectory(prefix="ns_whole_crust_bundle_") as temporary:
        staging = Path(temporary)
        for source in sorted(campaign_path.parent.rglob("*")):
            if not source.is_file():
                continue
            relative = source.relative_to(campaign_path.parent)
            parts = relative.parts
            if parts[0] == "models":
                if primary_only and parts[1] not in ids:
                    continue
                # Non-admitted results, if present, retain diagnostics only.
                if parts[1] not in ids and len(parts) > 2 and parts[2] == "compose":
                    continue
            if primary_only and parts[0] == "sources" and "__css" in str(relative):
                continue
            copy_file(source, staging, relative)
        for source in sorted((ROOT/"python_frontend/compose_eos").glob("*.py")):
            copy_file(source, staging, "code/python_frontend/compose_eos/"+source.name)
        for relative in ("fortran_eos/app/whole_crust_kernel.f90", "scripts/build_whole_crust_20261117.py",
                         "scripts/package_whole_crust_20261117.py", "scripts/diagnose_whole_crust_sources_20261117.py",
                         "python_frontend/tests/test_whole_crust.py",
                         "python_frontend/tests/test_whole_crust_native.py"):
            copy_file(ROOT/relative, staging, "code/"+relative)
        for name in ("SUMMARY.json", "NEXT_STEP_B.json", "RUN_RECORD.json"):
            source = ROOT/"data/whole_crust_20261117"/name
            if source.is_file():
                copy_file(source, staging, "provenance/"+name)
        for source in sorted((ROOT/"docs/whole_crust_20261117").glob("*.md")):
            if not primary_only or source.name == "USAGE_KO.md":
                copy_file(source, staging, "docs/"+source.name)
        if not primary_only:
            for directory in (ROOT/"data/whole_crust_20261117/source_audit", ROOT/"data/whole_crust_20261117/validation"):
                if directory.is_dir():
                    for source in sorted(directory.rglob("*")):
                        if source.is_file() and source.read_bytes()[:4] != b"\x7fELF":
                            copy_file(source, staging, "provenance/"+str(source.relative_to(directory.parent)))
            symbolic = ROOT/"data/whole_crust_20261117/SYMBOLIC_SURFACE.json"
            if symbolic.is_file():
                copy_file(symbolic, staging, "provenance/"+symbolic.name)
            for old in sorted(campaign_path.parent.parent.glob("production*")):
                previous = old/"CAMPAIGN.json"
                if old == campaign_path.parent or not previous.is_file():
                    continue
                old_record = json.loads(previous.read_text())
                if old_record.get("status") != "FAILED":
                    continue
                for source in sorted(old.rglob("*")):
                    relative = source.relative_to(old)
                    if source.is_file() and "compose" not in relative.parts and source.read_bytes()[:4] != b"\x7fELF":
                        copy_file(source, staging, "provenance/failed_attempts/"+old.name+"/"+str(relative))
        for name in ("WHOLE_CRUST_OVERVIEW.png", "WHOLE_CRUST_OVERVIEW.pdf", "PLOT_DATA.json"):
            copy_file(output/name, staging, "figures/"+name)
        save(staging/"CATALOG.json", catalog)
        (staging/"CATALOG.csv").write_text(catalog_csv(catalog))
        (staging/"evaluate.py").write_text(EVALUATOR)
        (staging/"requirements-consumer.txt").write_text("numpy\nscipy\n")
        (staging/"README.md").write_text(
            "# Whole-crust A: portable EOS bundle\n\n"
            f"Included tables: {len(chosen)}. Scientific baselines: 1. "
            f"Original production tables: {campaign['model_count']}.\n\n"
            "Read docs/USAGE_KO.md. The primary is a reference-assisted pressure-primary outer + "
            "SPG_M1–SLy4 smooth barotrope. Other included tables are declared construction sensitivities.\n\n"
            "Install NumPy and SciPy, then run `python evaluate.py --list` or "
            "`python evaluate.py --pressure 0 0.0001 1`. No native executable is needed for evaluation.\n\n"
            "REPRESENTATION.json and its consumer provide exact P=0. The CoMPOSE export contains only "
            "positive densities. Use the recorded piecewise laws and one-sided CSS interface choice. "
            "High-density extrapolation is rejected; upper density is 0.5 fm^-3.\n\n"
            "The 17 public outer isotope tabulation gaps are explicitly effective regularizations. "
            "The package does not contain TOV or oscillation results.\n\n"
            "SHA256SUMS.txt and MANIFEST.json cover every member except those two self-describing files. "
            "Native execution commands retain their historical paths as provenance; portable evaluation "
            "uses only relative paths. Full regeneration still requires the original repository.\n")
        inventory = {str(p.relative_to(staging)): dict(bytes=p.stat().st_size, sha256=sha(p))
                     for p in sorted(staging.rglob("*")) if p.is_file()}
        save(staging/"MANIFEST.json", dict(schema="ns-cowling.bundle-members.v1", files=inventory))
        (staging/"SHA256SUMS.txt").write_text("".join(v["sha256"]+"  "+k+"\n" for k, v in inventory.items()))
        archive = output/archive_name
        with zipfile.ZipFile(archive, "x", zipfile.ZIP_DEFLATED, compresslevel=6) as zipped:
            for path in sorted(staging.rglob("*")):
                if path.is_file():
                    zipped.write(path, str(path.relative_to(staging)))
        verified = verify_archive(archive, catalog)
    return dict(file=archive_name, bytes=archive.stat().st_size, sha256=sha(archive), **verified)


def verify_archive(archive, catalog):
    with zipfile.ZipFile(archive) as zipped:
        names = zipped.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate archive member")
        manifest = json.loads(zipped.read("MANIFEST.json"))["files"]
        if set(names) != set(manifest) | {"MANIFEST.json", "SHA256SUMS.txt"}:
            raise ValueError("archive member inventory mismatch")
        for name, record in manifest.items():
            data = zipped.read(name)  # zipfile checks CRC while reading.
            if len(data) != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
                raise ValueError("archive member byte/hash mismatch: "+name)
        accepted_tables = [n for n in names if n.startswith("models/") and n.endswith("/compose/eos.thermo")]
        if len(accepted_tables) != catalog["packaged_table_count"]:
            raise ValueError("archive accepted table count differs from catalog")
        with tempfile.TemporaryDirectory(prefix="ns_whole_crust_consumer_") as temporary:
            target = Path(temporary)
            zipped.extractall(target)
            env = dict(os.environ)
            env.pop("PYTHONPATH", None)
            check = subprocess.run([sys.executable, str(target/"evaluate.py"), "--pressure", "0", "0.0001", "1"],
                                   cwd=target, env=env, capture_output=True, text=True, timeout=45)
            if check.returncode:
                raise ValueError("standalone archive consumer failed: "+check.stderr)
            result = json.loads(check.stdout)
            if not Path(result["consumer_source"]).is_relative_to(target):
                raise ValueError("standalone verification imported the source checkout")
            if any(result["states"][0][k] != 0 for k in ("n", "p", "epsilon")):
                raise ValueError("packaged consumer did not reach the exact surface")
            outside = subprocess.run([sys.executable, str(target/"evaluate.py"), "--pressure", "1e20"],
                                     cwd=target, env=env, capture_output=True, text=True, timeout=45)
            if outside.returncode == 0 or "outside whole-crust contract support" not in outside.stderr:
                raise ValueError("packaged consumer did not reject high-density extrapolation")
    return dict(members=len(names), integrity="CRC_BYTE_COUNT_SHA256_VERIFIED",
                accepted_tables=len(accepted_tables), standalone_consumer="PASS",
                exact_zero_pressure="PASS", out_of_domain_rejection="PASS")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("new output directory required")
    campaign_path = args.campaign.resolve()
    if args.output.resolve().is_relative_to(campaign_path.parent):
        raise ValueError("deliverable directory must lie outside the production input directory")
    campaign, models = read_production(campaign_path)
    args.output.mkdir(parents=True)
    catalog = catalog_for(campaign_path, campaign, models, package_kind="full_campaign")
    save(args.output/"CATALOG.json", catalog)
    (args.output/"CATALOG.csv").write_text(catalog_csv(catalog))
    make_figure(models, args.output)
    archives = [bundle(campaign_path, campaign, models, args.output, primary_only=choice)
                for choice in (True, False)]
    receipt = dict(schema="ns-cowling.whole-crust-packaging.v1", status="VERIFIED",
                   created_at_utc=datetime.now(timezone.utc).isoformat(),
                   campaign_sha256=sha(campaign_path), package_script_sha256=sha(Path(__file__)),
                   table_count=len(models), scientific_baseline_count=1, archives=archives,
                   plot_evaluator="actual continuous piece consumer", plot_sha256={
                       name:sha(args.output/name) for name in
                       ("WHOLE_CRUST_OVERVIEW.png", "WHOLE_CRUST_OVERVIEW.pdf", "PLOT_DATA.json")})
    save(args.output/"PACKAGING_RECEIPT.json", receipt)
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
