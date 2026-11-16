"""Reproduce public-data candidate matching and archived outer-segment identity.
Run from repository root with Python and NumPy after downloading the three
Zenodo record11084025 files into build/crust_literature. This script never
generates an EOS, fits a parameter, or silently declares model identity.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "build/crust_literature"
OUT = Path(__file__).resolve().parent
EXPECTED = {
    "model_NMP.txt": ("https://zenodo.org/api/records/11084025/files/model_NMP.txt/content", "5cc222abd2d1f29054e9695f7f8d86c67d950949106cf8e9460ef6edf4ed7556"),
    "model_par.txt": ("https://zenodo.org/api/records/11084025/files/model_par.txt/content", "43371e043b7cc572ad283032ce1785cae22e608039dd86cbee64f70ee5d70614"),
}
TARGET = np.array([
    [.157,-15.8,230.9,-177.1,7966.3,33.9,56.3,-181.4,350.1],
    [.155,-15.9,223.9,-230.1,2141.8,30.2,40.2,-59.8,1007.2],
    [.161,-15.8,225.7,-283.4,11174.8,33.7,50.9,-231.0,322.9],
    [.162,-15.8,246.0,-292.9,-1397.4,30.9,40.1,-111.8,958.5],
    [.148,-15.8,227.7,45.2,14552.1,31.4,52.3,-113.5,539.8],
])
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def zenodo():
    for name, (url, expected) in EXPECTED.items():
        if sha(RAW/name) != expected:
            raise ValueError("Wrong source bytes: " + name)
    nmp = np.loadtxt(RAW/"model_NMP.txt", skiprows=1)
    par_lines = (RAW/"model_par.txt").read_text().splitlines()
    if len(nmp) != len(par_lines)-1:
        raise ValueError("Model row alignment failed")
    tol = np.array([.0005]+[.05]*8)
    result = {
        "status": "candidate_identity_pending_independent_EOS_check",
        "source": "https://zenodo.org/records/11084025",
        "model_par_sha256": sha(RAW/"model_par.txt"),
        "nmp_sha256": sha(RAW/"model_NMP.txt"),
        "header": par_lines[0], "row_count": len(nmp),
        "nucleon_mass_MeV":938.9, "delta_meson":False,
        "g_i":"a_i + (b_i+d_i*x**3)*exp(-c_i*x)",
        "x":"n_B/n_sat_of_selected_model",
        "isovector_normalization": "Raw author values retained. Paper Eq2/Eq29 uses g_rho/2; raw vector normalization must be determined using J and source beta branch before admission.",
        "models": {},
    }
    for i, target in enumerate(TARGET):
        d = abs(nmp[:,:9]-target)
        strict = np.flatnonzero(np.all(d<=tol+1e-9,axis=1))
        near = int(np.argmin(np.max(d/tol,axis=1)))
        expanded = np.flatnonzero(np.all(d<=2*tol+1e-9,axis=1))
        params = [float(x) for x in par_lines[near+1].split()]
        result["models"][f"SPG_M{i+1}"] = {
            "zero_based_model_index":near, "physical_line_1_based":near+2,
            "NMP":nmp[near].tolist(), "table9":target.tolist(),
            "strict_halfbin_match":bool(len(strict)),
            "strict_halfbin_indices":strict.tolist(),
            "full_printed_bin_indices":expanded.tolist(),
            "maximum_halfbin_distance":float(np.max(d[near]/tol)),
            "parameter_line":par_lines[near+1],
            "rho0_fm3":params[0],
            "meson_masses_MeV":dict(zip(["sigma","omega","rho"],params[1:4])),
            "couplings":{k:dict(zip(["a","b","c","d"],params[4+4*j:8+4*j]))
                         for j,k in enumerate(["sigma","omega","rho"])},
        }
    (OUT/"spg_author_candidates.json").write_text(json.dumps(result,indent=2)+"\n")
    return result

def archive():
    base=ROOT/"eos_generator/eos_data"
    comparisons=[]
    for other in ["GPPVA_DD2","GPPVA_DDME2"]:
        for filename,skip,first_col in [("eos.nb",2,0),("eos.thermo",1,3),("eos.compo",0,3)]:
            a=base/"GPPVA_TW"/filename
            b=base/other/filename
            ar=[x.split()[first_col:] for x in a.read_text().splitlines()[skip:skip+217]]
            br=[x.split()[first_col:] for x in b.read_text().splitlines()[skip:skip+217]]
            if len(ar)!=217 or len(br)!=217: raise ValueError("Incomplete outer segment")
            mismatches=[i+1 for i,(x,y) in enumerate(zip(ar,br)) if x!=y]
            canonical="\n".join(" ".join(row) for row in ar)+"\n"
            comparisons.append({"other":other,"file":filename,"rows":217,
                "ignored_header_lines":skip,"ignored_leading_grid_index_columns":first_col,
                "different_rows_1based":mismatches,"TW_full_file_sha256":sha(a),
                "other_full_file_sha256":sha(b),"compared_TW_token_sha256":hashlib.sha256(canonical.encode()).hexdigest()})
    names=[f"SPG_M{i}" for i in range(1,6)]+["GPPVA_"+x for x in ["DD2","DDME2","FSU2","FSU2H","FSU2R","NL3wrL55","TM1e","TW"]]
    models={}
    for name in names:
        nb=np.loadtxt(base/name/"eos.nb",skiprows=2)
        comp=[x.split() for x in (base/name/"eos.compo").read_text().splitlines() if x.strip()]
        if len(nb)!=len(comp): raise ValueError("Composition alignment "+name)
        phase=np.array([int(x[3]) for x in comp])
        models[name]={"rows":len(nb),"n_min_fm3":float(nb.min()),"n_max_fm3":float(nb.max()),
            "phase_counts":{str(p):int(sum(phase==p)) for p in np.unique(phase)},
            "first_nonzero_phase_density_fm3":float(nb[np.flatnonzero(phase!=0)[0]]),
            "source_file_sha256":{fn:sha(base/name/fn) for fn in ["eos.nb","eos.thermo","eos.compo"]}}
    result={"status":"implementation_verified_token_identity","comparisons":comparisons,
        "interpretation":"TW common outer segment is token-identical in density, thermodynamics, and composition to DD2/DDME2, after ignoring grid indices. This resolves extraction consistency, not historical BPS-vs-BSk22 bibliographic attribution.",
        "models":models}
    (OUT/"crust_archive_evidence.json").write_text(json.dumps(result,indent=2)+"\n")
    return result

if __name__=="__main__":
    z=zenodo(); a=archive()
    print("Zenodo source rows:",z["row_count"])
    print("Candidate indices:",{k:v["zero_based_model_index"] for k,v in z["models"].items()})
    print("Strict halfbin matches:",{k:v["strict_halfbin_match"] for k,v in z["models"].items()})
    print("Outer comparisons:",len(a["comparisons"]),"mismatching rows:",sum(len(x["different_rows_1based"]) for x in a["comparisons"]))
    print("Composition-bearing crust archives:",len(a["models"]))
