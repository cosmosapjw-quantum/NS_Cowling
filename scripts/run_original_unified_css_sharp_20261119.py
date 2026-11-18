#!/usr/bin/env python3
"""Official full GPPVA unified EOS, minimal first-law corrected split, CSS sharp native stellar control.

Source-n, P and raw bytes are immutable; corrections to epsilon and mu, explicit
outer-to-inner causal gap, and new CSS interface are labelled derived. No Maxwell
coexistence is invented where source mu bands do not overlap.
"""
from __future__ import annotations
import hashlib,json,sys,time,traceback
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos import read_compose,match_crust_core
from compose_eos.whole_crust import WholeCrustBarotrope,PiecewisePolytrope
from compose_eos.maxwell_v2 import _segment_eosdata
from run_official_original_unified_controls_20261119 import _view,OUT as PROVIDER_OUT
from build_original_fsu2_same_rmf_three_layer_20261119 import mode_result

DEST=ROOT/"data/cowling_maxwell_native_20261119/official_unified_css_sharp"
SOURCE=PROVIDER_OUT/"ORIGINAL_CONTROLS.json"

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def one(model,control):
    source_dir=PROVIDER_OUT/model/"raw"
    if not control.get("independent_thermo_matches_zip_member"):
        raise ValueError("original EOS thermo source independent readback invalid")
    if control.get("provider_zip_checksum_matches"):
        checksum="provider ZIP hash matched"
    else:
        checksum="provider reported ZIP sha256 does NOT match actual zip; member thermo independently verified"
    read_dir=ROOT/control["source_projection"]["projection_path"] if control.get("source_projection") else source_dir
    data=read_compose(read_dir)
    phase=np.asarray(data.Iphase)
    low_ids=np.flatnonzero(phase==7)
    high_ids=np.flatnonzero(phase==6)
    assert len(low_ids)>8 and len(high_ids)>8 and low_ids[-1]+1==high_ids[0]
    i,j=int(low_ids[0]),int(high_ids[0])
    end=int(np.searchsorted(data.nb,.49,side="right"))
    original_upper=data.nb[i:end].copy()
    corrected=data.copy()
    f=PiecewisePolytrope(data.nb[i:end],data.p[i:end],
                       data.eps[i]/data.nb[i])
    corrected.Q3[i:end]=f.mu/corrected.mn-1
    corrected.Q6[i:end]=f.F/corrected.mn-1
    corrected.Q7[i:end]=corrected.Q6[i:end]
    low=_view(corrected,i,j)
    high=_view(corrected,j,end)
    matched,report=match_crust_core(low,high,
         "OFFICIAL_"+model+"_FIRST_LAW_INNER",
         "OFFICIAL_"+model+"_FIRST_LAW_CORE",
         "css",40,drop_composition=True)
    if not report["consistency"]["overall"]["pass_"]:
        raise ValueError("native CSS consistency reported false")
    jump=int(report["css_transition_low_index"])
    if not (0<=jump<matched.n-1 and matched.nb[jump+1]>matched.nb[jump]
          and np.isclose(matched.p[jump],matched.p[jump+1],rtol=2e-11)):
        raise ValueError("native CSS pressure or density jump inconsistent")
    original_contract_path=PROVIDER_OUT/model/"thermodynamically_corrected_smooth"/"REPRESENTATION.json"
    smooth=json.loads(original_contract_path.read_text())
    # Only replace the core-containing upper table by a pressure/chemical-potential
    # native CSS doublet. The outer and low-density connector are IDENTICAL.
    upper=dict(n=matched.nb.tolist(),p=matched.p.tolist(),
               epsilon=matched.eps.tolist(),sharp_interfaces=[jump])
    whole_contract=dict(surface=smooth["surface"],outer=smooth["outer"],
                        lower_bridge=smooth["lower_bridge"],upper=upper)
    if not np.array_equal(matched.nb[:i-i+1],low.nb[:1]):
        raise RuntimeError("first original inner density not preserved")
    whole=WholeCrustBarotrope.from_contract(whole_contract)
    outdir=DEST/model
    outdir.mkdir(parents=True,exist_ok=True)
    save(outdir/"REPRESENTATION.json",whole_contract)
    native=mode_result(whole,"OFFICIAL_"+model,"corrected_original_css",.30)
    raw_smooth=control["corrected_smooth_native"]
    smooth_modes={x["label"]:x["frequency_hz"] for x in raw_smooth["modes"]}
    delta_modes={conversion:{
           y["label"]:float(y["frequency_hz"]-smooth_modes[y["label"]])
           for y in details["modes"] if y["label"] in smooth_modes}
           for conversion,details in native["modes"].items()}
    row=dict(model=model,status="CSS_THREE_LAYER_NATIVE_STELLAR_CHECKED",
       original_unified_zip_sha256=control["provider_eos_zip_sha256"],
       archive_checksum_issue=checksum,
       immutable_original_thermo_sha256=control["source_file_hashes"]["eos.thermo"],
       source_direct_inner_nP_preserved_before_native_css=True,
       source_projection=control.get("source_projection"),
       first_law_reconstructed_epsilon_mu=True,
       original_inner_core_support=[float(original_upper[0]),float(original_upper[-1])],
       raw_upper_modified_physical_laws="energy/chemical-potential correction; source n,P at original knots unchanged unless native CSS truncation near transition",
       corresponding_smooth_source= str(original_contract_path.relative_to(ROOT)),
       corresponding_smooth_hash=digest(original_contract_path),
       interface=dict(index=jump,pressure=float(matched.p[jump]),
         density=[float(matched.nb[jump]),float(matched.nb[jump+1])],
         energy_jump=float(matched.eps[jump+1]-matched.eps[jump])),
       native_fortran= native,
       delta_vs_smooth_same_source=dict(
          mass_Msun=native["native_tov"]["mass_msun"]-raw_smooth["mass_msun"],
          radius_km=native["native_tov"]["R_km"]-raw_smooth["radius_km"],
          frequencies_hz=delta_modes),
       whole_contract_sha256=digest(outdir/"REPRESENTATION.json"))
    return row

def main():
    c=json.loads(SOURCE.read_text())
    report=dict(schema="ns_cowling.official_unified_css_sharp.v1",
         source_control_path=str(SOURCE.relative_to(ROOT)),
         source_control_sha256=digest(SOURCE),
         physics="original source n,P + thermodynamic correction of F,mu + declared causal outer-inner connector + CSS finite jump",
         original_unified_provider_binary_id="NOT_VERIFIED_ZIP_CHECKSUM",
         cases={},status="RUNNING")
    for model in ("FSU2","TM1e","DDME2"):
        try:
            row=one(model,c["models"][model])
        except Exception as e:
            row=dict(status="BLOCKED",type=type(e).__name__,message=str(e),
                     traceback=traceback.format_exc()[-1800:])
        report["cases"][model]=row
        print("OFFICIAL_CSS_SHARP",model,row["status"],
              row.get("native_fortran",{}).get("native_tov"),flush=True)
        save(DEST/"PARTIAL.json",report)
    report["status"]="COMPLETED_WITH_EXPLICIT_CASE_STATUSES"
    save(DEST/"RESULT.json",report)

if __name__=="__main__":main()
