#!/usr/bin/env python3
"""Normalize legacy CompOSE nuclear-group abundance fields into official Y_I.

This utility never overwrites its input tree. It creates a new output tree and
writes a per-model provenance manifest. The reconstruction is accepted only
when baryon closure and charge neutrality both pass.

Supported scope: one cold density slice and at most one averaged nuclear-group
quadruple per grid point. Multiple nuclear groups are rejected explicitly.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path


TOL = 2.0e-3
NUCLEON_NEUTRON = 10
NUCLEON_PROTON = 11
LEPTON_ELECTRON = 0
LEPTON_MUON = 1


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _parse_compo_row(line: str, source: Path, row_number: int):
    fields = line.split()
    if not fields:
        return line, None

    if len(fields) < 6:
        raise ValueError(f"{source}:{row_number}: CompOSE row has fewer than 6 fields")

    try:
        indices = [int(fields[i]) for i in range(3)]
        phase = int(fields[3])
        npairs = int(fields[4])
    except Exception as exc:
        raise ValueError(f"{source}:{row_number}: invalid row header") from exc
    if npairs < 0:
        raise ValueError(f"{source}:{row_number}: negative Npairs")

    cursor = 5
    pairs: list[tuple[int, float]] = []
    seen = set()
    for _ in range(npairs):
        if cursor + 1 >= len(fields):
            raise ValueError(f"{source}:{row_number}: truncated particle-pair list")
        pid = int(fields[cursor])
        value = float(fields[cursor + 1])
        if pid in seen:
            raise ValueError(f"{source}:{row_number}: duplicate particle index {pid}")
        if not (value == value and abs(value) != float("inf")):
            raise ValueError(f"{source}:{row_number}: non-finite fraction for particle {pid}")
        seen.add(pid)
        pairs.append((pid, value))
        cursor += 2

    if cursor >= len(fields):
        raise ValueError(f"{source}:{row_number}: missing Nquad")
    nquad = int(fields[cursor])
    cursor += 1
    if nquad < 0:
        raise ValueError(f"{source}:{row_number}: negative Nquad")
    if nquad > 1:
        raise ValueError(
            f"{source}:{row_number}: Nquad={nquad}; current EOSData normalizer "
            "supports at most one averaged nuclear group per point"
        )
    if len(fields) != cursor + 4 * nquad:
        raise ValueError(f"{source}:{row_number}: inconsistent Nquad/trailing fields")

    # Empty composition rows are legal CompOSE syntax and must not be
    # interpreted as pure-neutron matter.
    if npairs == 0 and nquad == 0:
        return line, None

    groups = []
    for _ in range(nquad):
        gid = int(fields[cursor])
        Abar = float(fields[cursor + 1])
        Zbar = float(fields[cursor + 2])
        raw_y = float(fields[cursor + 3])
        groups.append((gid, Abar, Zbar, raw_y))
        cursor += 4

    # The all-zero quadruple is a legacy placeholder, not a nuclear group.
    if groups and groups[0][1:] == (0.0, 0.0, 0.0):
        groups = []

    pair_map = dict(pairs)
    Yn = pair_map.get(NUCLEON_NEUTRON, 0.0)
    Yp = pair_map.get(NUCLEON_PROTON, 0.0)
    Ye = pair_map.get(LEPTON_ELECTRON, 0.0)
    Ymu = pair_map.get(LEPTON_MUON, 0.0)
    adjustments = {"Yn": 0.0, "Yp": 0.0}
    normalization = None

    if groups:
        gid, Abar, Zbar, raw_y = groups[0]
        if not (Abar > 0 and 0 < Zbar <= Abar):
            raise ValueError(
                f"{source}:{row_number}: invalid nuclear group Abar={Abar}, Zbar={Zbar}"
            )

        # Use provider n/p abundances only when they close both the baryon sum
        # and charge-neutrality relation. Otherwise, reconstruct the missing
        # free-neutron fraction only when the supplied free-proton fraction is
        # zero within tolerance. Do not guess in the remaining ambiguous cases.
        Yh_candidate = 1.0 - Yn - Yp
        charge_candidate = Yp + (Zbar / Abar) * Yh_candidate - Ye - Ymu
        if (
            -TOL <= Yh_candidate <= 1.0 + TOL
            and abs(charge_candidate) <= TOL
        ):
            Yh = Yh_candidate
        elif Yp <= TOL:
            Yh = (Ye + Ymu - Yp) * (Abar / Zbar)
            if not (-TOL <= Yh <= 1.0 + TOL):
                raise ValueError(
                    f"{source}:{row_number}: charge-derived Yh={Yh:.8g} outside [0,1]"
                )
            Yn = 1.0 - Yp - Yh
            adjustments["Yn"] = Yn - pair_map.get(NUCLEON_NEUTRON, 0.0)
        else:
            raise ValueError(
                f"{source}:{row_number}: n/p fractions cannot be reconciled "
                f"(baryon-derived Yh={Yh_candidate:.8g}, charge residual={charge_candidate:.3e})"
            )

        YI = Yh / Abar
        baryon_residual = abs(Yn + Yp + Abar * YI - 1.0)
        charge_residual = abs(Yp + Zbar * YI - Ye - Ymu)
        normalization = {
            "group_id": gid,
            "Abar": Abar,
            "Zbar": Zbar,
            "raw_group_field": raw_y,
            "reconstructed_group_baryon_fraction": Yh,
            "normalized_compose_number_fraction": YI,
            "baryon_residual": baryon_residual,
            "charge_residual": charge_residual,
        }
        groups = [(gid, Abar, Zbar, YI)]
    else:
        baryon_residual = abs(Yn + Yp - 1.0)
        charge_residual = abs(Yp - Ye - Ymu)
        if baryon_residual > TOL or charge_residual > TOL:
            if Yp <= TOL:
                Yp = Ye + Ymu
                Yn = 1.0 - Yp
                adjustments["Yn"] = Yn - pair_map.get(NUCLEON_NEUTRON, 0.0)
                adjustments["Yp"] = Yp - pair_map.get(NUCLEON_PROTON, 0.0)
            else:
                raise ValueError(
                    f"{source}:{row_number}: homogeneous fractions fail closure "
                    f"(baryon={baryon_residual:.3e}, charge={charge_residual:.3e})"
                )
        baryon_residual = abs(Yn + Yp - 1.0)
        charge_residual = abs(Yp - Ye - Ymu)
        normalization = {
            "group_id": None,
            "normalized_compose_number_fraction": None,
            "baryon_residual": baryon_residual,
            "charge_residual": charge_residual,
        }

    if not (-TOL <= Yn <= 1.0 + TOL and -TOL <= Yp <= 1.0 + TOL):
        raise ValueError(f"{source}:{row_number}: reconstructed n/p fraction outside [0,1]")
    if baryon_residual > TOL or charge_residual > TOL:
        raise ValueError(
            f"{source}:{row_number}: post-normalization closure failed "
            f"(baryon={baryon_residual:.3e}, charge={charge_residual:.3e})"
        )

    # Update or add n/p species only when their reconstructed abundance differs
    # from the implicit zero. This preserves the provider ordering otherwise.
    rewritten_pairs = []
    pair_ids = {pid for pid, _ in pairs}
    for pid, value in pairs:
        if pid == NUCLEON_NEUTRON:
            value = Yn
        elif pid == NUCLEON_PROTON:
            value = Yp
        rewritten_pairs.append((pid, value))
    if NUCLEON_NEUTRON not in pair_ids and Yn > 1e-14:
        rewritten_pairs.append((NUCLEON_NEUTRON, Yn))
    if NUCLEON_PROTON not in pair_ids and Yp > 1e-14:
        rewritten_pairs.append((NUCLEON_PROTON, Yp))

    output = [str(indices[0]), str(indices[1]), str(indices[2]), str(phase), str(len(rewritten_pairs))]
    for pid, value in rewritten_pairs:
        output.extend((str(pid), f"{value:.16E}"))
    output.append(str(len(groups)))
    for gid, Abar, Zbar, YI in groups:
        output.extend((str(gid), f"{Abar:.16E}", f"{Zbar:.16E}", f"{YI:.16E}"))

    record = {
        "indices": indices,
        "phase": phase,
        "adjustments": adjustments,
        "normalization": normalization,
    }
    return " ".join(output), record


def normalize_compo_file(source: Path, target: Path) -> dict:
    lines = source.read_text(encoding="utf-8").splitlines()
    output_lines = []
    records = []
    max_baryon = 0.0
    max_charge = 0.0
    reconstructed_rows = 0
    nuclear_groups = 0
    homogeneous_rows = 0

    for row_number, line in enumerate(lines, 1):
        normalized, record = _parse_compo_row(line, source, row_number)
        output_lines.append(normalized)
        if record is None:
            continue
        norm = record["normalization"]
        if norm is None:
            continue
        max_baryon = max(max_baryon, float(norm["baryon_residual"]))
        max_charge = max(max_charge, float(norm["charge_residual"]))
        if norm["group_id"] is None:
            homogeneous_rows += 1
        else:
            nuclear_groups += 1
        if any(abs(v) > 0 for v in record["adjustments"].values()):
            reconstructed_rows += 1
        if norm["baryon_residual"] > TOL or norm["charge_residual"] > TOL:
            raise ValueError(f"{source}:{row_number}: closure residual exceeds tolerance")
        records.append(record)

    target.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return {
        "source_sha256": sha256(source),
        "normalized_sha256": sha256(target),
        "rows": len(lines),
        "nuclear_group_rows": nuclear_groups,
        "homogeneous_rows": homogeneous_rows,
        "rows_reconstructed": reconstructed_rows,
        "max_baryon_residual": max_baryon,
        "max_charge_residual": max_charge,
    }


def repair_single_point_grid(model_dir: Path, filename: str, target_index: int) -> dict | None:
    path = model_dir / filename
    if not path.is_file():
        return None
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 3:
        return None
    values = lines[2].split()
    if len(values) != 1:
        return None
    old = [lines[0].strip(), lines[1].strip()]
    if old == [str(target_index), str(target_index)]:
        return None
    lines[0] = str(target_index)
    lines[1] = str(target_index)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"file": filename, "old_indices": old, "new_index": target_index}


def normalize_tree(input_root: Path, output_root: Path) -> dict:
    input_root = input_root.resolve()
    output_root = output_root.resolve()
    if not input_root.is_dir():
        raise FileNotFoundError(f"input root does not exist: {input_root}")
    if input_root == output_root or input_root in output_root.parents:
        raise ValueError("output root must not be inside the input tree")
    if output_root == output_root.parent or output_root in input_root.parents:
        raise ValueError("output root must not be an ancestor of the input tree")
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"output root must not exist or must be empty: {output_root}")
    output_root.parent.mkdir(parents=True, exist_ok=True)
    if output_root.exists():
        output_root.rmdir()

    staging = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.stage-", dir=output_root.parent))
    stage_tree = staging / output_root.name
    try:
        shutil.copytree(input_root, stage_tree)
        models = []
        grid_repairs = []
        for source_model in sorted(p for p in input_root.iterdir() if p.is_dir()):
            target_model = stage_tree / source_model.name
            source_compo = source_model / "eos.compo"
            if source_compo.is_file():
                try:
                    result = normalize_compo_file(source_compo, target_model / "eos.compo")
                    result["model"] = source_model.name
                    result["status"] = "normalized"
                except Exception as exc:
                    # A partially rewritten composition file would falsely
                    # advertise complete microscopic coverage. Omit it from
                    # the normalized copy and retain the raw source untouched.
                    target_compo = target_model / "eos.compo"
                    if target_compo.exists():
                        target_compo.unlink()
                    result = {
                        "model": source_model.name,
                        "status": "composition_omitted_unresolved_rows",
                        "source_sha256": sha256(source_compo),
                        "rows": len(source_compo.read_text(encoding="utf-8").splitlines()),
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                models.append(result)

            # Match eos.t/eos.yq singleton grid indices to the actual thermo
            # row indices. Only one-point grids are changed; multi-point tables
            # are left untouched.
            thermo = source_model / "eos.thermo"
            if thermo.is_file():
                data_rows = [
                    line.split() for line in thermo.read_text(encoding="utf-8").splitlines()
                    if line.strip() and not line.lstrip().startswith("#")
                ]
                if len(data_rows) >= 2 and len(data_rows[1]) >= 3:
                    iT, iYq = int(float(data_rows[1][0])), int(float(data_rows[1][2]))
                    for fn, idx in (("eos.t", iT), ("eos.yq", iYq)):
                        repair = repair_single_point_grid(target_model, fn, idx)
                        if repair:
                            repair["model"] = source_model.name
                            grid_repairs.append(repair)

        manifest = {
            "schema": "compose-eos-composition-normalization-v1",
            "authority": "CompOSE Reference Manual v3.00, Sec. 4.2.3, Eqs. (4.20)-(4.22)",
            "manual_url": "https://compose.obspm.fr/download/pdf/manual_v3.00.pdf",
            "input_root_label": input_root.name,
            "output_root_label": output_root.name,
            "method": (
                "Reconstruct Y_h from baryon closure and verify charge neutrality; "
                "write official CompOSE nuclear-group number abundance Y_I=Y_h/Abar. "
                "Only infer missing free nucleon fractions under the explicitly coded "
                "closure rule. If a model contains unresolved rows, its normalized "
                "eos.compo is omitted as a whole; raw input remains untouched."
            ),
            "tolerance": TOL,
            "models": models,
            "grid_index_repairs": grid_repairs,
        }
        (stage_tree / "COMPOSITION_NORMALIZATION_MANIFEST.json").write_text(
            json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        output_root.parent.mkdir(parents=True, exist_ok=True)
        stage_tree.rename(output_root)
        return manifest
    finally:
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = normalize_tree(args.input_root, args.output_root)
    print(json.dumps({
        "output_root": str(args.output_root),
        "models_normalized": len(manifest["models"]),
        "grid_index_repairs": len(manifest["grid_index_repairs"]),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
