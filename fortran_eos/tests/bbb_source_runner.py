#!/usr/bin/env python3
"""Compile/run the BBB native source diagnostic; Python only handles files/results."""
from pathlib import Path
import argparse, csv, hashlib, json, os, re, subprocess, time, tarfile
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/research_blockers_20261116/bbb"
SOURCE = ROOT / "data/research_generator_20261116/bbb2_published_tables.json"

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def cached_source(cache, name, manifest):
    entry = next(item for item in manifest["source_files"] if item["path"] == name)
    path = cache / name
    cache.mkdir(parents=True, exist_ok=True)
    if not path.is_file():
        if name == "ns_rev.tex":
            archive = cached_source(cache, "astro-ph_9707277.tex.tar", manifest)
            with tarfile.open(archive) as package:
                raw = package.extractfile("ns_rev.tex").read()
        else:
            with urlopen(entry["url"], timeout=60) as response:
                raw = response.read()
        if hashlib.sha256(raw).hexdigest() != entry["sha256"]:
            raise ValueError("downloaded source hash mismatch: " + name)
        path.write_bytes(raw)
    if digest(path) != entry["sha256"]:
        raise ValueError("cached source hash mismatch: " + name)
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=DEFAULT)
    p.add_argument("--source-cache", type=Path, default=ROOT / "build/blocker_bbb_sources")
    p.add_argument("--build", type=Path, default=Path("/tmp/ns_bbb_source_native_20261116"))
    args = p.parse_args()
    out, build = args.output.resolve(), args.build.resolve()
    out.mkdir(parents=True, exist_ok=True)
    build.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((DEFAULT / "SOURCE_MANIFEST.json").read_text())
    cache = args.source_cache.resolve()
    original = cached_source(cache, "ns_rev.tex", manifest)
    paper = original.read_text()
    block = paper.split(r"\begin{tabular}{|ccccc|}")[1].split(r"\end{tabular}")[0]
    rows = []
    for line in block.splitlines():
        if re.match(r"\s*\d", line):
            fields = [s.strip() for s in line.split(r"\\")[0].split("&")]
            if len(fields) == 5:
                rows.append(fields)
    source = json.loads(SOURCE.read_text())
    table3 = paper.split("% Table 3")[1].split(r"\end{tabular}")[0]
    parsed3 = []
    for line in table3.splitlines():
        if re.match(r"\s*\d", line):
            values = [t.strip() for t in line.split(r"\\")[0].split("&")]
            if len(values) == 4:
                parsed3.append(list(map(float, values)))
    assert parsed3 == list(map(list, zip(source["n"], source["table3_proton_fraction"],
        source["table3_mass_density_1e14_g_cm3"], source["table3_pressure_1e34_dyn_cm2"])))
    assert len(rows) == 16
    numeric = []
    for i, row in enumerate(rows):
        n, es, en = float(row[0]), float(row[2]), float(row[4])
        assert [n, es, en] == [source[k][i] for k in ("n", "symmetric_EA", "neutron_EA")]
        def halfwidth(token):
            digits = len(token.split(".")[1]) if "." in token else 0
            return 0.5 * 10**(-digits)
        numeric.append([n, es, en, halfwidth(row[2]), halfwidth(row[4]),
            source["table3_proton_fraction"][i], source["table3_pressure_1e34_dyn_cm2"][i],
            source["table3_mass_density_1e14_g_cm3"][i]])
    table = out / "published_source_with_precision.dat"
    table.write_text("".join(" ".join(format(x, ".17g") for x in row)+"\n" for row in numeric))
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    code = ROOT / "fortran_eos/tests/bbb_source_diagnostics.f90"
    deps = [ROOT/"fortran_eos/src/compose_eos_types.f90",
            ROOT/"fortran_eos/src/compose_eos_core_native.f90",code]
    exe = build/"bbb_source_diagnostics"
    command = ["gfortran","-std=f2018","-O2","-ffree-line-length-none","-fimplicit-none",
               "-fcheck=all","-J",str(build),"-I",str(build),*(str(x) for x in deps),"-o",str(exe)]
    start=time.time()
    result=subprocess.run(command,cwd=build,env=env,text=True,capture_output=True)
    (out/"compile.log").write_text(result.stdout+result.stderr)
    result.check_returncode()
    run=[str(exe),str(table),str(out/"native_source_diagnostics.csv")]
    result=subprocess.run(run,cwd=build,env=env,text=True,capture_output=True)
    (out/"native_run.log").write_text(result.stdout+result.stderr)
    result.check_returncode()
    records=[{k:float(v) for k,v in row.items()} for row in csv.DictReader((out/"native_source_diagnostics.csv").open())]
    nodes=records[:16]; microscopic=[r for r in records if r["n"]<=.8]
    with (out/"microscopic_source_only_native.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(["n_fm3","epsilon_MeV_fm3","pressure_MeV_fm3","Yp"])
        writer.writerows([r["n"],r["epsilon_truncated"],r["p_truncated"],r["yp_truncated"]] for r in sorted(microscopic,key=lambda r:r["n"]))
    contamination=[dict(n=r["n"],p_full=r["p"],p_source_only=r["p_truncated"],
       relative_pressure_difference=r["p_truncated"]/r["p"]-1,
       yp_full=r["yp"],yp_source_only=r["yp_truncated"],stencil_max=r["stencil_max"])
       for r in microscopic if r["n"]>=.7]
    # RNS is an independent delivery of final beta-equilibrium EOS, not source nuclear energies.
    rns_path=cached_source(cache, "rns_eosBBB2", manifest)
    rns_rows=[]
    for line in rns_path.read_text().splitlines()[1:]:
        v=list(map(float,line.split()))
        if len(v)==4 and v[3]*1e-39>=.079999:
            n=v[3]*1e-39
            nearest=min(nodes,key=lambda r:abs(r["n"]-n))
            rns_rows.append(dict(n=n,nearest_table3_n=nearest["n"],rns_pressure_dyn_cm2=v[1],
                table3_pressure_dyn_cm2=source["table3_pressure_1e34_dyn_cm2"][nodes.index(nearest)]*1e34,
                matches_table3_node=abs(nearest["n"]-n)<1e-5))
    summary=dict(
        claim="Published nuclear-energy reconstruction; no microscopic BHF kernel rerun.",
        source_json_sha256=digest(SOURCE),original_tex_sha256=digest(original),
        source_cache=str(cache),driver_sha256=digest(Path(__file__)),
        tex_table1_matches_locked_json=True,tex_table3_matches_locked_json=True,table1_nodes=len(nodes),
        production_max_abs_pressure_difference=max(abs(r["native_p_difference"]) for r in records),
        production_max_abs_energy_density_difference=max(abs(r["native_epsilon_difference"]) for r in records),
        pchip_vs_local_cubic_max_relative_pressure_difference=max(abs(r["p_local_cubic"]/r["p"]-1) for r in records),
        table3_pressure_max_relative_difference_all_nodes=max(abs(r["p"]/r["p_table3"]-1) for r in nodes),
        table3_pressure_max_relative_difference_excluding_008=max(abs(r["p"]/r["p_table3"]-1) for r in nodes[1:]),
        table3_yp_max_absolute_difference=max(abs(r["yp"]-r["yp_table3"]) for r in nodes),
        source_only_truncation=contamination,
        selected_rounding_sensitivity=[dict(n=r["n"],p=r["p"],p_corner_min=r["p_corner_min"],
            p_corner_max=r["p_corner_max"],yp_corner_min=r["yp_corner_min"],yp_corner_max=r["yp_corner_max"],
            corner_count=int(r["corner_count"])) for r in nodes if r["n"] in [.08,.16,.8,1.5]],
        rounding_interpretation="All vertices of each active local PCHIP input precision box; finite corner sample envelope, not rigorous nonlinear extrema, probability interval or physical EOS uncertainty.",
        interpolation_interpretation="Spread between two declared representations, not an upper error bound.",
        source_only_policy="First nine Table 1 rows, endpoint slope rebuilt at 0.8; no EOS evaluation above 0.8. Production representation unchanged.",
        full_source_policy="0.08..1.5 published source-energy reconstruction; high source rows extrapolated; PCHIP at 0.7<n<=0.8 uses the 0.9 source row.",
        native_output_rows=len(records),rns_core_rows=rns_rows,checks_stdout=result.stdout,
        execution=dict(command_compile=command,command_run=run,exit_code=result.returncode,
            elapsed_seconds=time.time()-start,thread_limits={k:env[k] for k in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"]}),
        dependency_sha256={str(x.relative_to(ROOT)):digest(x) for x in [*deps,SOURCE]},
        executable_sha256=digest(exe))
    (out/"result.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k not in ("rns_core_rows","execution","dependency_sha256")},indent=2))

if __name__=="__main__":
    main()
