#!/usr/bin/env python3
"""Acquire pinned primary HFB22/AME2016 data and construct the outer-crust mass input."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
SOURCES={
 "hfb22-dat":("https://www.astro.ulb.ac.be/bruslib/nucdata/hfb22-dat","a5484cd56cf93a694843a6f8e62cd3452248a4cd0a1d7e3bd9a99a489bfd8a69"),
 "mass16.txt":("https://www-nds.iaea.org/amdc/ame2016/mass16.txt","2167f57a2a98331e4649b2dd2b658a9006ed4fba1975729ebfe52a42b4b9218a"),
}
# Welker etal PRL119192502(2017), TableI; keV converted to MeV.
# 78Cu: Penning trap result. MRTOF alternative -44.819MeV is sensitivity-tested.
CU={75:-54.47004,76:-50.98155,77:-48.8628,78:-44.772,79:-42.408}

def build(source: Path, output: Path):
    provenance={}
    for name,(url,sha) in SOURCES.items():
        path=source/name
        if not path.exists():
            source.mkdir(parents=True,exist_ok=True)
            with urllib.request.urlopen(url,timeout=60) as response:
                path.write_bytes(response.read())
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest!=sha: raise ValueError(f"source identity mismatch: {name}")
        provenance[name]={"url":url,"sha256":digest,"bytes":path.stat().st_size}
    nuclei={}
    for line in (source/"hfb22-dat").read_text().splitlines():
        fields=line.split()
        if len(fields)>9 and fields[0].isdigit():
            z,a=int(fields[0]),int(fields[1])
            if (z,a) in nuclei: raise ValueError("duplicate HFB22 nuclide")
            nuclei[z,a]=dict(Z=z,A=a,mass_excess_mev=float(fields[9]),
                            rch_fm=float(fields[4]),source="HFB22")
    measured=0
    for line in (source/"mass16.txt").read_text().splitlines()[39:]:
        try:
            z,a=int(line[9:14]),int(line[14:19])
            raw=line[28:42]
            if "#" in raw or z<=0: continue
            excess=float(raw)/1000
        except ValueError:
            continue
        # Rch=0 is a conservative lower bound on the positive finite-size
        # correction for nuclei outside HFB22. A selected such nucleus refuses.
        row=nuclei.setdefault((z,a),dict(Z=z,A=a,rch_fm=0.0))
        row.update(mass_excess_mev=excess,source="AME2016_measured")
        measured+=1
    for a,excess in CU.items():
        nuclei[29,a].update(mass_excess_mev=excess,source="Welker2017_TableI")
    output.mkdir(parents=True,exist_ok=True)
    table=output/"hfb22_ame2016_cu.csv"
    with table.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=["Z","A","mass_excess_mev","rch_fm","source"])
        writer.writeheader()
        writer.writerows(nuclei[key] for key in sorted(nuclei))
    provenance.update(
        catalogue_sha256=hashlib.sha256(table.read_bytes()).hexdigest(),
        nuclei=len(nuclei),ame_measured_read=measured,
        zero_radius_lower_bound_candidates=sum(r["rch_fm"]==0 for r in nuclei.values()),
        hfb22_mass_units="atomic mass excess MeV; author file rounded0.01MeV",
        radius="HFB22 rms charge radius fm, including for experimentally replaced masses",
        cu={"url":"https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.119.192502",
            "table":"I","mass_excess_mev":CU,
            "78Cu_choice":"Penning trap -44.772MeV; MRTOF alternative -44.819MeV",
            "79Cu_uncertainty_mev":0.105},
        mass_override_priority=["Welker2017_TableI","AME2016_measured_no_hash","HFB22"],
        interpretation="nuclear masses and radii only; no CompOSE EOS rows used")
    (output/"provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    return provenance

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources",type=Path,default=ROOT/"build/outer_crust_sources")
    parser.add_argument("--output",type=Path,default=ROOT/"data/outer_crust")
    args=parser.parse_args()
    print(json.dumps(build(args.sources,args.output),indent=2))
