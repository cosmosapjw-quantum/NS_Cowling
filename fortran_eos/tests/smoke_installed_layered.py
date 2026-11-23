"""Exercise the installed wheel in a venv using public cell_000001 and its stored controls."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

from compose_eos import native
from compose_eos.native_tov_love import NativeTovLove


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    parser.add_argument("public_archive", type=Path)
    parser.add_argument("public_stars_csv", type=Path)
    args = parser.parse_args()
    if any(os.environ.get(key) for key in (
        "COMPOSE_EOS_FORTRAN_ROOT", "COMPOSE_EOS_FORTRAN_EXE", "COMPOSE_EOS_NATIVE_BUILD_DIR",
    )):
        raise ValueError("Run with native path overrides unset.")
    prefix = Path(sys.prefix)
    assert sys.prefix != sys.base_prefix
    root = native.project_root()
    assert root == prefix / "share" / "compose-eos"
    assert Path(native.__file__).is_relative_to(prefix)
    model = NativeTovLove()
    cache = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
    assert model.binary.is_relative_to(cache.resolve() / "compose_eos_fortran" / "installed")
    with args.public_stars_csv.open(newline="") as handle:
        reference = {row["request"]: row for row in csv.DictReader(handle)
                     if row["census_cell_id"] == "cell_000001" and float(row["max_dlogp"]) == .04}
    fields = ("R_km", "M_msun", "compactness", "k2", "Lambda", "y_surface",
              "central_pressure_MeV_fm3", "central_eps_MeV_fm3", "n_steps", "transition_count",
              "surface_pressure_MeV_fm3", "surface_epsilon_MeV_fm3", "surface_jump_scale")
    with tempfile.TemporaryDirectory() as staging, zipfile.ZipFile(args.public_archive) as archive:
        table = Path(staging) / "cell_000001"
        table.mkdir()
        input_hashes = {}
        for name in ("eos.t", "eos.nb", "eos.yq", "eos.thermo", "eos.compo"):
            file = table / name
            file.write_bytes(archive.read(f"cell_000001/{name}"))
            input_hashes[name] = sha256(file)
        result = {"Pc25": model.layered(table, 25., .04),
                  "M1.4": model.layered_target(table, 1.4, 5., 300., .04)}
        for request, row in result.items():
            assert row["converged"] and row["status"] == reference[request]["status"] == "ok"
            for name in fields:
                assert row[name] == float(reference[request][name]), (request, name, row[name])
        cli = prefix / "bin" / "compose-eos"
        imported = json.loads(subprocess.check_output([cli, "read", table], text=True))
        assert imported["status"] == "ok" and imported["n"] == 1028
        catalog = json.loads(subprocess.check_output([cli, "generators"], text=True))
        assert "APR_APR" in catalog
    stored_vendor = root / "fortran_eos" / "vendor" / "odepack"
    required = [f"fortran_eos/vendor/odepack/{name}.f" for name in ("vode", "dgefa", "dgesl", "dgbfa", "dgbsl")]
    required += ["fortran_eos/run-nv-profile.sh", "fortran_eos/run-causal-audit.sh",
                 "fortran_eos/tools/cold_causal_endpoint_audit.f90", "fortran_eos/app/nv_profile_verify.f90",
                 "data/eos_original_papers_20261123/research_nv/nv1973_original_source.dat"]
    with zipfile.ZipFile(args.wheel) as archive:
        for rel in required:
            path, = [name for name in archive.namelist() if name.endswith("/share/compose-eos/" + rel)]
            assert (root / rel).read_bytes() == archive.read(path)
    compiler = shutil.which(os.environ.get("FC", "gfortran"))
    report = dict(
        wheel=args.wheel.name, wheel_size=args.wheel.stat().st_size, wheel_sha256=sha256(args.wheel),
        python=sys.version.split()[0], compiler=subprocess.check_output([compiler, "--version"], text=True).splitlines()[0],
        linker_libraries=os.environ.get("NATIVE_LIBS", "-llapack -lblas"), native_path_overrides=False,
        native_sources="sys.prefix/share/compose-eos", native_build="user cache", cli_rows=imported["n"],
        generator_catalog_models=len(catalog), public_matching_archive_sha256=sha256(args.public_archive),
        reference_csv_sha256=sha256(args.public_stars_csv), public_source="2020-MVCD-D1MS-full",
        cell="cell_000001", input_sha256=input_hashes, wheel_members_checked=required,
        wheel_odepack_sha256={p.name: sha256(p) for p in sorted(stored_vendor.glob("*.f"))},
        installed_oscillation_f90_count=len(list((root / "oscillation_source").rglob("*.f90"))),
        compared_fields=list(fields), equality="exact parsed numeric equality against stored CSV",
        max_dlogp=.04, result=result,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
