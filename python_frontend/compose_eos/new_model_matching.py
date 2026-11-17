"""Explicit new-cohort provenance policy and lossless XMLSLZ thermo projection.

Scientific kernels are unchanged. Provider Q1--Q7 tokens are never repaired.
The Nadd payload is retained as exact token strings and original source lines.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from compose_eos.api import read_compose, _thermo_only_view
from compose_eos.maxwell_v2 import _segment_eosdata

NEW_NATIVE = ("NL3_PROVIDER_PRECISION", "PK1_CANONICAL", "TM1_CANONICAL",
              "DDMEX_PROVIDER_ROUNDED_N0", "DDMEX_ORIGINAL_TABLE_RHO0")
NEW_PUBLIC = ("XMLSLZ_NL3", "XMLSLZ_PK1", "XMLSLZ_TM1", "XMLSLZ_DDME-X", "XMLSLZ_PKDD")
NEW_PUBLIC_CORES = NEW_PUBLIC[:4]
OLD_CRUSTS = tuple("SPG_M"+str(i) for i in range(1, 6)) + tuple(
    "GPPVA_"+s for s in ("DD2","DDME2","TW","FSU2","FSU2R","FSU2H","NL3wrL55","TM1e"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_origin_policy(crust, core):
    """Allow exactly the five explicitly authorized new campaign origin routes."""
    if crust["role"] != "crust" or core["role"] != "core":
        raise ValueError("crust/core role mismatch")
    c, k = crust["origin"], core["origin"]
    cm, km = crust["model"], core["model"]
    if c == "native_generated" and cm in OLD_CRUSTS and k == "new_native" and km in NEW_NATIVE:
        return "old_native_crust__new_native_core"
    if c == "new_public" and cm in NEW_PUBLIC:
        if k == "public_original" and km == "GM_GM1":
            return "new_public_crust__old_public_GM1"
        if k == "new_public" and km in NEW_PUBLIC_CORES:
            return "new_public_crust__new_public_core"
        if (k == "new_native" and km in NEW_NATIVE) or (k == "native_generated" and km == "GM_GM1"):
            return "reference_assisted_hybrid"
    raise ValueError("pair outside explicit new-model origin policy")


def parse_provider(directory):
    """Validate frozen cold slice rows, Nadd layout and phase/geometry agreement."""
    directory = Path(directory)
    thermo_lines = (directory/"eos.thermo").read_text().splitlines()
    header = thermo_lines[0].split()
    if len(header) != 3:
        raise ValueError("unexpected thermo header")
    rows = [line.split() for line in thermo_lines[1:] if line.strip()]
    compo = [line.split() for line in (directory/"eos.compo").read_text().splitlines() if line.strip()]
    nb_tokens = (directory/"eos.nb").read_text().split()
    first, last = map(int, nb_tokens[:2])
    n = np.asarray([float(s.replace("D","E")) for s in nb_tokens[2:]])
    if len(rows) != len(n) or len(compo) != len(n) or last-first+1 != len(n):
        raise ValueError("source row/axis count mismatch")
    expected = 5 if directory.name == "XMLSLZ_TM1" else 7
    phases = []
    for i, (row, cp) in enumerate(zip(rows, compo)):
        if len(row) < 11 or int(row[10]) != expected or len(row) != 11+expected:
            raise ValueError("Nadd/payload length differs from declared provider layout")
        if tuple(map(int,row[:3])) != (1,first+i,1) or row[:3] != cp[:3]:
            raise ValueError("source row index/order mismatch")
        geom = float(row[11+(0 if expected == 5 else 2)])
        phase = int(cp[3])
        if phase not in (3,2,1,0,-2,-3) or geom != phase:
            raise ValueError("provider phase/geometry mismatch")
        phases.append(phase)
    q = np.asarray([[float(s.replace("D","E")) for s in row[3:10]] for row in rows])
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(n)) or np.any(np.diff(n)<=0):
        raise ValueError("nonfinite or unordered provider table")
    phases = np.asarray(phases, dtype=int)
    uniform = np.flatnonzero(phases == 0)
    if not len(uniform) or not np.array_equal(uniform, np.arange(uniform[0],len(n))):
        raise ValueError("homogeneous phase0 must be the terminal contiguous provider segment")
    return {"header":header, "rows":rows, "raw_lines":thermo_lines[1:], "n":n, "q":q,
            "phases":phases, "nadd":expected, "core_start":int(uniform[0]), "first":first}


def project_provider(directory, destination, role):
    """Lossless thermodynamic projection; whole raw payload remains auditable."""
    directory, destination = Path(directory), Path(destination)
    parsed = parse_provider(directory)
    if role not in ("crust","core"):
        raise ValueError("role must be crust or core")
    destination.mkdir(parents=True, exist_ok=False)
    projection = destination/"projection"
    projection.mkdir()
    for name in ("eos.nb","eos.t","eos.yq"):
        shutil.copyfile(directory/name,projection/name)
    rows = parsed["rows"]
    (projection/"eos.thermo").write_text(" ".join(parsed["header"])+"\n"+
        "".join(" ".join(row[:10])+" 0\n" for row in rows))
    # Native reader carries phase-only eos.compo; abundances remain in source.
    phase_map = {p:i+1 for i,p in enumerate(sorted(set(parsed["phases"].tolist())))}
    (projection/"eos.compo").write_text("".join(
        " ".join(row[:3])+f" {phase_map[int(p)]} 0 0\n"
        for row,p in zip(rows,parsed["phases"])))
    sidecar = {
        "schema":"ns_cowling.xmlslz_exact_auxiliary.v1",
        "source":str(directory), "source_thermo_sha256":sha(directory/"eos.thermo"),
        "nadd":parsed["nadd"],
        "columns":(["d","Z","A","Rd","RW"] if parsed["nadd"]==5 else ["ce2","cs2","d","Z","A","Rd","RW"]),
        "working_phase_dictionary":{str(v):k for k,v in phase_map.items()},
        "thermodynamic_reconstruction":"NONE",
        "rows":[{"source_row_zero_based":i,"indices":row[:3],"q_tokens":row[3:10],
                 "auxiliary_tokens":row[11:],"raw_thermo_line":line,"source_iphase":int(p)}
                for i,(row,line,p) in enumerate(zip(rows,parsed["raw_lines"],parsed["phases"]))]}
    sidecar_path=destination/"SOURCE_AUXILIARY.json"
    sidecar_path.write_text(json.dumps(sidecar,indent=2)+"\n")
    data=read_compose(projection)
    if data.n != len(rows) or not np.array_equal(data.nb,parsed["n"]):
        raise ValueError("native projected readback density mismatch")
    for j in range(7):
        if not np.array_equal(getattr(data,"Q"+str(j+1)),parsed["q"][:,j]):
            raise ValueError("native projection changed provider Q columns")
    a,b=(0,parsed["core_start"]) if role=="crust" else (parsed["core_start"],data.n)
    out=_thermo_only_view(_segment_eosdata(data,(0,a,b)))
    out.beta_equilibrated=bool(np.max(np.abs(out.Q5))<=1e-8)
    # Provenance points to immutable original, not the transport projection.
    out.source_dir=str(directory)
    return out, {
        "source_start_zero_based":a,"source_stop_exclusive":b,
        "basis":"raw nonzero geometrical phases" if role=="crust" else "raw homogeneous Iphase0",
        "q1_q7_native_projection_bitwise_equal":True,
        "auxiliary_sidecar":str(sidecar_path),"auxiliary_sha256":sha(sidecar_path),
        "nadd_original":parsed["nadd"],"thermodynamic_reconstruction":"NONE",
        "muQ_vs_muE_max_abs_MeV":float(np.max(np.abs(out.mub-(out.eps+out.p)/out.nb))),
        "q5_max_abs":float(np.max(np.abs(out.Q5))),
        "public_numerical_precision":"source printed precision; no invented grid refinement"}
