#!/usr/bin/env python3
"""Export original source EOS two-phase TOV *input* arrays only.

No Python stellar structure ODE is evaluated. Fortran computes
PCHIP thermodynamics and complete one-sided TOV solutions.
The original compressed EOS contract and stellar center identity
are separately hash-checked. G and c are retained via the source
system's geom-km conversion constants.
"""
from __future__ import annotations
import hashlib,io,json,sys,zipfile,gzip
from pathlib import Path
from math import pi
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import _one_phase_grid
from compose_eos.stellar_background import MEV_FM3_TO_GEOM_KM2,MSUN_GEOM_KM

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_c15_tov_convergence_20261121"
CASES={
 "M2":("SPG_M2","maxwell",.30385128926615257,1.3999999857354968,14.420419623651872),
 "M3":("SPG_M3","css",.192510767556007,.7000000086373105,14.59794179069283)
}

def run():
 raw=SOURCE.read_bytes()
 if hashlib.sha256(raw).hexdigest()!=SHA:raise RuntimeError("EOS source SHA drift")
 with zipfile.ZipFile(io.BytesIO(raw)) as outer:
  with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as data:
   for short,(model,kind,nc,ref_m,ref_r) in CASES.items():
    rep=data.read(f"FRESH_5x2_WHOLE_STARS/{model}/{kind}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(rep)))
    r0=.001
    hc=float(eos.enthalpy_at_density(nc))
    center=eos.state_at_enthalpy(hc,"high")
    pc=center["p"]*MEV_FM3_TO_GEOM_KM2
    ec=center["epsilon"]*MEV_FM3_TO_GEOM_KM2
    h0=hc-2*pi*(ec/3+pc)*r0*r0
    m0=4*pi*ec*r0**3/3
    surface_h=3e-8
    trans=eos.transitions[0]["h"]
    phases=[]
    for hi,lo in ((h0,trans),(trans,surface_h)):
      grid=_one_phase_grid(h0,surface_h,hi,lo,eos.break_h,
                           2049,extras=lo==surface_h)
      mid=(hi+lo)*.5
      y=np.empty((len(grid),2))
      for i,s in enumerate(grid):
        h=hi if i==0 else lo if i==len(grid)-1 else h0-s*s
        st=eos.state_at_enthalpy(h,eos.side_for_segment(h,mid))
        y[i]=[st["p"],st["epsilon"]]
      y*=MEV_FM3_TO_GEOM_KM2
      phases.append((grid,y))
    controls=[
       (3e-10,3e-11,min(np.sqrt(h0-surface_h)/200,.004/(2*np.sqrt(h0-surface_h))),10),
       (0.,3e-11,min(np.sqrt(h0-surface_h)/200,.004/(2*np.sqrt(h0-surface_h))),10),
       (1e-10,3e-11,min(np.sqrt(h0-surface_h)/200,.002/(2*np.sqrt(h0-surface_h))),22),
       (3e-10,3e-11,min(np.sqrt(h0-surface_h)/200,.004/(2*np.sqrt(h0-surface_h))),22)
    ]
    dest=OUT/f"C15_ORIGINAL_TOV_SOURCE_{short}.dat"
    with dest.open("w") as f:
      f.write(f"{len(phases[0][0])} {len(phases[1][0])} {len(controls)}\n")
      f.write(f"{r0:.17e} {m0:.17e} {ref_m*MSUN_GEOM_KM:.17e} {ref_r:.17e}\n")
      for grid,y in phases:
         np.savetxt(f,np.column_stack((grid,y)),fmt="%.17e")
      for rtol,atol,dsmax,mf in controls:
        f.write(f"{rtol:.17e} {atol:.17e} {dsmax:.17e} {mf}\n")
    full=dest.read_bytes()
    checksum=hashlib.sha256(full).hexdigest()
    zipped=gzip.compress(full,compresslevel=9,mtime=0)
    zg=dest.with_suffix(".dat.gz")
    zg.write_bytes(zipped)
    dest.unlink()
    print("C15_SOURCE_TOV_INPUT",short,
      "CONTRACT_SHA",hashlib.sha256(rep).hexdigest(),
      "UNCOMPRESSED_SHA",checksum,
      "GZIP_SHA",hashlib.sha256(zipped).hexdigest(),
      "PHASE_GRID",len(phases[0][0]),len(phases[1][0]),
      "INPUT",zg,flush=True)

if __name__=="__main__":
 run()
