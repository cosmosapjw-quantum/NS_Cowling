"""Bounded APR single-branch matching demonstration; never erase source phase tags."""
import hashlib
import json
from dataclasses import fields
from pathlib import Path

import numpy as np

from compose_eos import read_compose, write_compose, match_crust_core
from compose_eos.native import executable

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).with_name("apr_low_branch")
OUT.mkdir(exist_ok=False)
CRUST = ROOT / "data/research_matching_20261117/spg/native_129_repaired/SPG_M1/generated"
CORE = ROOT / "data/research_matching_20261117/integration/real_pairs_v2/core/APR_APR/generated"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

crust, full = read_compose(CRUST), read_compose(CORE)
source_hashes = {side: {p.name: sha(p) for p in folder.glob("*") if p.is_file()}
                 for side, folder in (("crust", CRUST), ("core", CORE))}
end = int(np.flatnonzero(full.Iphase != full.Iphase[0])[0])
core = full.copy()
for field in fields(core):
    value = getattr(core, field.name)
    if isinstance(value, np.ndarray) and value.ndim == 1 and value.size == full.n:
        setattr(core, field.name, value[:end].copy())
assert core.source_dir == str(CORE.resolve())
assert np.array_equal(core.Iphase, full.Iphase[:end])
assert len(set(core.Iphase.tolist())) == 1
record = {
    "scope": "Contiguous APR low-phase subset only; not an assembled full multiphase APR stellar EOS",
    "binary": str(executable()), "binary_sha256": sha(executable()),
    "source_hashes": source_hashes,
    "core_rows_full": full.n, "selected_rows_zero_based": [0, end - 1],
    "selected_density_fm3": [float(core.nb[0]), float(core.nb[-1])],
    "selected_phase_id": int(core.Iphase[0]), "selected_rows": core.n,
    "phase_tags_modified": False, "source_energy_shift_MeV": 0.0,
    "full_multiphase_assembly": "UNSUPPORTED_PHASE_BOUNDARY_RETENTION",
    "methods": [],
}
for method in ("css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"):
    item = {"method": method}
    try:
        matched, info = match_crust_core(crust, core, method=method, n_bridge=40, drop_composition=True)
        i2 = info["i2"]
        count = core.n - i2
        original_f = (core.Q6[i2:] + 1.) * core.mn
        retained_f = (matched.Q6[-count:] + 1.) * matched.mn
        f_error = float(np.max(np.abs(retained_f-original_f)))
        assert f_error <= 64*np.finfo(float).eps*max(1., float(np.max(np.abs(original_f))))
        assert np.array_equal(matched.nb[-count:], core.nb[i2:])
        assert np.array_equal(matched.Q1[-count:], core.Q1[i2:])
        assert info["dm_crust"] == 0.0
        dest = OUT / method
        write_compose(dest, matched)
        back = read_compose(dest)
        for name in ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Iphase"):
            assert np.allclose(getattr(matched, name), getattr(back, name), rtol=5e-13, atol=1e-14)
        item.update(status="CONSTRUCTED_DIAGNOSTIC_ONLY" if info["diagnostic_only"] else "CONSTRUCTED_VERIFIED",
                    rows=matched.n, core_endpoint_index_in_selected_branch=i2,
                    core_endpoint_density_fm3=float(core.nb[i2]),
                    retained_core_F_max_abs_error_MeV=f_error,
                    retained_core_density_and_pressure_identity=True,
                    native_compose_readback="PASS", report=info)
    except Exception as exc:
        item.update(status="FAILED", error_type=type(exc).__name__, error=str(exc))
    record["methods"].append(item)
record["source_bytes_unchanged"] = all(sha(folder/name) == value
    for side, folder in (("crust", CRUST), ("core", CORE))
    for name, value in source_hashes[side].items())
record["status"] = "COMPLETE" if all(x["status"]=="CONSTRUCTED_VERIFIED" for x in record["methods"]) else "PARTIAL"
(OUT / "RECEIPT.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")
print(json.dumps({k:record[k] for k in ("status", "selected_rows_zero_based", "selected_density_fm3", "selected_phase_id", "source_bytes_unchanged")}))
print(json.dumps([{k:v for k,v in x.items() if k!="report"} for x in record["methods"]], indent=2))
