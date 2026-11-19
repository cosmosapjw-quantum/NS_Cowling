#!/usr/bin/env python3
"""Inspect, fingerprint and isolate 23 UNMODIFIED nucleardatapy crust raw tables.

The source license is CC BY-NC-ND 4.0: no adapted original nuda data are
distributed. The model physics is not regenerated from its original
microscopic EDF. Report only numerical intake counts and status.
"""
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/"data/eos_missing_data_recalc_20261119/external_original"
WHEEL=Path("/home/cosmosapjw/ns_cowling_wt/scratch/all_origin_maxwell_css_20261118/nuda/nucleardatapy-1.0.2-py3-none-any.whl")
EXPECTED_WHEEL_SHA="6d8f57cf056696f39ee45fe4d516bc4eae7a2461c9aba5d83c2012b499116d61"

def sha(raw):return hashlib.sha256(raw).hexdigest()
def save(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
def parse_numeric_rows(raw):
    rows=[]
    for ln in raw.decode(errors="replace").splitlines():
        line=ln.strip()
        if not line or line.startswith("#"):continue
        try:values=[float(x.replace("D","E").replace("d","e")) for x in line.split()]
        except ValueError:continue
        if len(values)>=2 and all(np.isfinite(values)):
            rows.append(values)
    return rows

def main():
    binary=WHEEL.read_bytes()
    if sha(binary)!=EXPECTED_WHEEL_SHA:
        raise ValueError("original nucleardatapy-1.0.2 wheel SHA failed")
    DEST.mkdir(exist_ok=True,parents=True)
    rows=dict(schema="NS_COWLING_NUDA_1.0.2_RAW_UNCHANGED_SOURCE_AUDIT",
        wheel_sha256=sha(binary),license="CC BY-NC-ND 4.0",
        upstream="https://pypi.org/project/nucleardatapy/1.0.2/",
        source_work="Gogny Mondal Viñas Centelles De 2020; GMRS Grams Margueron Reddy Somasundaram 2022; Negele Vautherin 1973; Pearson et al BSk 2018",
        original_source_bytes_preserved=True,original_datasets={},total=0,
        source_tier="UNALTERED_UPSTREAM_RAW_NOT_FULL_SURFACE_THERMODYNAMIC_EOS",
        physical_derivative_models_unverified=True)
    with zipfile.ZipFile(io.BytesIO(binary)) as zp:
        paths=sorted(x for x in zp.namelist() if
                     x.startswith("nucleardatapy/data/crust/") and x.endswith(".dat"))
        # 23 named source models plus the separately supplied D1M* crust+core
        # extension, which is NOT an independent 24th nuclear functional.
        if len(paths)!=24 or "nucleardatapy/data/crust/2020-MVCD-D1MS-full.dat" not in paths:
            raise ValueError("expected 23 native sources plus D1MS-full extension, got "+str(len(paths)))
        licensefile=zp.read("nucleardatapy-1.0.2.dist-info/licenses/LICENSE")
        (DEST/"SOURCE_LICENSE_CC_BY_NC_ND_4_0.txt").write_bytes(licensefile)
        for rel in paths:
            raw=zp.read(rel)
            original_name=Path(rel).name
            (DEST/original_name).write_bytes(raw)
            parsed=parse_numeric_rows(raw)
            if not parsed:raise ValueError("no numerical values in "+rel)
            first=np.asarray([v[0] for v in parsed])
            columns=sorted(set(len(row) for row in parsed))
            family=("Gogny_MVCD" if "2020-MVCD" in rel else
                    "GMRS_2022" if "2022-GMRS" in rel else
                    "BSk_2018" if "2018-PCPF" in rel else
                    "Negele_Vautherin_1973" if "1973-Negele" in rel else "UNKNOWN")
            report=dict(path_in_wheel=rel,filename=original_name,
                sha256=sha(raw),bytes=len(raw),numeric_row_count=len(parsed),
                observed_numeric_column_counts=columns,
                first_column_min=float(min(first)),first_column_max=float(max(first)),
                first_column_monotone=bool(np.all(np.diff(first)>0)),
                raw_bytes_unchanged=True,producer=family,
                physical_model_recomputed=False,source_data_only=True,
                full_source_outer_to_inner_matched=False)
            rows["original_datasets"][original_name]=report
            print("NUDA_ORIGINAL_INTAKE",original_name,"rows",len(parsed),
                  "columns",columns,"sha256",report["sha256"][:12],flush=True)
    rows["total"]=len(rows["original_datasets"])
    rows["nominal_independent_source_models"]=rows["total"]-1
    rows["additional_nonindependent_extension"]="2020-MVCD-D1MS-full.dat"
    rows["model_groups"]={g:sum(v["producer"]==g for v in rows["original_datasets"].values())
                          for g in sorted(set(v["producer"] for v in rows["original_datasets"].values()))}
    rows["status"]="23_UNMODIFIED_NAMED_SOURCES_PLUS_ONE_D1MS_FULL_EXTENSION_NOT_MICROSCOPIC_REGENERATIONS"
    save(DEST/"NUDA_SOURCE_AUDIT.json",rows)
    print("NUDA_AUDIT",rows["total"],rows["model_groups"],flush=True)

if __name__=="__main__":main()
