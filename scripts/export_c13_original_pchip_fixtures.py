#!/usr/bin/env python3
"""Frontend export of C12 source-locked stellar background PCHIP oracle.

Uses the already versioned immutable original-M2/M3 one-sided C12 input
bytes. All expected PCHIP derivatives/values are independently evaluated
with SciPy only to verify a standalone native-Fortran implementation.
"""
from __future__ import annotations
import gzip,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
SOURCE_DIR=ROOT/"data/cowling_fortran_source_exact_root_c12_20261121"
OUT=ROOT/"data/cowling_native_pchip_c13_20261121"
EXPECTED={
  "M2":"eee499eb50de4300f9987e05bf46a4fee3ddcbfb2e9b9d99deb9ac8a423c35e5",
  "M3":"f6a7abb566ae68b8ac5f838d410771823634a98024007ca02264bd6e44a0d3b1",
}

def run():
    records={}
    for key in ("M2","M3"):
        gz=SOURCE_DIR/f"ORIGINAL_FORTRAN_ROOT_{key}.dat.gz"
        source=gzip.decompress(gz.read_bytes())
        sha=hashlib.sha256(source).hexdigest()
        if sha!=EXPECTED[key]:raise RuntimeError("C12 original-input source drift "+key)
        lines=source.decode().splitlines()
        no,nu,jump=map(int,lines[0].split())
        off=3+no+nu
        nh,nl=map(int,lines[off].split())
        start=off+3
        records[key]=dict(source_uncompressed_sha256=sha,source_c12_gzip_sha256=hashlib.sha256(gz.read_bytes()).hexdigest(),phases={})
        for name,count in (("high",nh),("low",nl)):
            part=np.asarray([list(map(float,line.split())) for line in lines[start:start+count]],dtype=np.float64)
            start+=count
            assert part.shape==(count,11)
            original_x=part[:,0]
            assert np.all(np.diff(original_x)>0)
            indices=np.unique(np.round(np.linspace(0,count-1,257)).astype(int))
            x=original_x[indices]
            y=part[indices,1:6]
            oracle=PchipInterpolator(x,y,axis=0,extrapolate=False)
            slopes=oracle(x,nu=1)
            queries=np.linspace(x[0],x[-1],97)
            qvalues=oracle(queries)
            f=OUT/f"SOURCE_PCHIP_C13_{key}_{name}.dat"
            with f.open("w") as stream:
                stream.write(f"{len(x)} {len(queries)}\n")
                np.savetxt(stream,np.column_stack((x,y,slopes)),fmt="%.17e")
                np.savetxt(stream,np.column_stack((queries,qvalues)),fmt="%.17e")
            raw=f.read_bytes()
            zipped=gzip.compress(raw,compresslevel=9,mtime=0)
            z=f.with_suffix(".dat.gz")
            z.write_bytes(zipped)
            f.unlink()
            records[key]["phases"][name]=dict(original_grid_points=count,
               sampled_grid_points=len(x),
               source_phase_first_s=float(original_x[0]),
               source_phase_last_s=float(original_x[-1]),
               oracle_sha256=hashlib.sha256(raw).hexdigest(),
               fixture_gzip_sha256=hashlib.sha256(zipped).hexdigest())
            print("ORIGINAL_PCHIP_FIXTURE",key,name,"SOURCE_SAMPLES",count,
                  "ORACLE_SAMPLES",len(x),"SHA",hashlib.sha256(raw).hexdigest(),flush=True)
    (OUT/"ORIGINAL_PCHIP_FIXTURE_SOURCES.json").write_text(
        json.dumps(records,indent=2,sort_keys=True,allow_nan=False)+"\n")

if __name__=="__main__":run()
