#!/usr/bin/env python3
"""Frontend-only export of immutable original EOS contract and reference probes.

No Python scientific solver is implemented here. Numerical sound-speed
laws are evaluated at runtime by the separate Fortran implementation.
"""
from __future__ import annotations
import hashlib,io,json,zipfile,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS

SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_fortran_source_cs2_20261120"

def run():
    orig=SOURCE.read_bytes()
    if hashlib.sha256(orig).hexdigest()!=SHA:raise RuntimeError("EOS source SHA drift")
    with zipfile.ZipFile(io.BytesIO(orig)) as container:
        inner=container.read("05_SPG5_PK1_10_WHOLE_STARS.zip")
    with zipfile.ZipFile(io.BytesIO(inner)) as container:
        for model,join in (("SPG_M2","maxwell"),("SPG_M3","css")):
            source_key=f"FRESH_5x2_WHOLE_STARS/{model}/{join}/REPRESENTATION.json"
            blob=container.read(source_key)
            eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(blob)))
            b=eos.barotrope
            out=OUT/f"FORTRAN_EOS_{model}_{join}.dat"
            limits=[0.,eos.h_surface_anchor,eos.h_outer_end,
                    eos.h_bridge_end,eos.transitions[0]["h"],eos.h_max]
            probes=[]
            for region,(a,c) in enumerate(zip(limits[:-1],limits[1:]),start=1):
                for h in np.linspace(a,c,31)[1:-1]:
                    probes.append((float(h),2,eos.state_at_enthalpy(float(h))["cs2_over_c2"],min(4,region)))
            for h in eos.break_h[::max(1,len(eos.break_h)//15)]:
                for change in (-3e-10,3e-10):
                    q=float(h+change)
                    if q<=0 or q>=eos.h_max or eos._transition(q) is not None:continue
                    reg=1 if q<=eos.h_surface_anchor else 2 if q<=eos.h_outer_end else 3 if q<=eos.h_bridge_end else 4
                    probes.append((q,2,eos.state_at_enthalpy(q)["cs2_over_c2"],reg))
            ht=eos.transitions[0]["h"]
            for side in (1,2):
                probes.append((ht,side,eos.state_at_enthalpy(ht,"low" if side==1 else "high")["cs2_over_c2"],4))
            outer=np.column_stack([b.outer.nb,b.outer.p,b.outer.F,
                                   np.r_[b.outer.gamma,b.outer.gamma[-1]]])
            upper=np.column_stack([b.upper_n,b.upper_p,b.upper_epsilon])
            bridge=[b.bridge.n_left,b.bridge.p_left,b.bridge.epsilon_left,
                    b.bridge.n_right,b.bridge.p_right,b.bridge.epsilon_right,b.bridge.kappa]
            surface=[b.surface.rest_energy_MeV,b.surface.n0,b.surface.p0,b.surface.gamma]
            with out.open("w") as f:
                f.write(f"{len(outer)} {len(upper)} {int(next(iter(b.sharp_interfaces)))}\n")
                f.write(" ".join(f"{v:.17e}" for v in surface)+"\n")
                f.write(" ".join(f"{v:.17e}" for v in bridge)+"\n")
                np.savetxt(f,outer,fmt="%.17e")
                np.savetxt(f,upper,fmt="%.17e")
                f.write(f"{len(probes)}\n")
                for h,side,c,region in probes:
                    f.write(f"{h:.17e} {side:d} {c:.17e} {region:d}\n")
            print(model,join,"contract_sha256",hashlib.sha256(blob).hexdigest(),
                  "probes",len(probes),"input_sha256",hashlib.sha256(out.read_bytes()).hexdigest(),
                  "file",out,flush=True)
if __name__=="__main__":run()
