"""Independent validation of species-restored finite seam candidates."""
import json,hashlib
from pathlib import Path
import numpy as np
from run_rest_mass_seam import restored_state,corrected_cell,bsk_state,DIFF
from compose_eos.candidate_gdfm_cld import reduced_cell,solve_two_phase
HERE=Path(__file__).resolve().parent
experiment=json.loads((HERE/"rest_mass_seam_experiment.json").read_text())
out={"scope":"finite first-law, beta re-solve, Hessian and common-tangent checks of explicitly distinct restored candidates","models":{}}
for model,d in experiment["models"].items():
    root=d["matches"]["energy_consistent"]["roots"][0];lo,hi=root["n_bsk"],root["n_spg"];p=root["pressure"];mu=root["spg_state"]["mu_B"]
    tangent=[]
    for n in np.linspace(lo,hi,7):
        b=bsk_state(float(n));s=restored_state(model,float(n))
        tangent.append({"n_B":float(n),"bsk_above_tangent":b["epsilon"]-mu*n+p,"restored_spg_above_tangent":s["epsilon"]-mu*n+p})
    firstlaw=[]
    for n in (.001,.002,hi):
        s=restored_state(model,float(n));errors=[]
        for relative in (2e-4,1e-4):
            h=n*relative
            ee=[restored_state(model,float(n+k*h))["epsilon"] for k in (-2,-1,1,2)]
            derivative=(ee[0]-8*ee[1]+8*ee[2]-ee[3])/(12*h)
            errors.append(derivative-s["mu_B"])
        firstlaw.append({"n_B":n,"derivative_errors_MeV":errors})
    z=np.array(root["spg_state"]["coordinates"]);n=hi
    corrected,g,_=corrected_cell(model,n,z);old_energy,oldg,_=reduced_cell(model,n,z)
    h=np.array([2e-5,2e-6,2e-6]);matrix=[]
    for j in range(3):
        q=np.zeros(3);q[j]=h[j]
        matrix.append((corrected_cell(model,n,z+q)[1]-corrected_cell(model,n,z-q)[1])/(2*h[j]))
    mat=np.asarray(matrix).T;mat=(mat+mat.T)/2;scale=np.array([1.,.1,.1])
    eig=np.linalg.eigvalsh(mat*scale[:,None]*scale[None,:])
    original,_=solve_two_phase(model,n)
    state=root["spg_state"]
    record={"common_tangent_probes":tangent,"first_law_checks":firstlaw,
        "min_tangent_gap_MeV_fm3":min(x[k] for x in tangent for k in ("bsk_above_tangent","restored_spg_above_tangent")),
        "max_firstlaw_error_MeV":max(abs(x) for r in firstlaw for x in r["derivative_errors_MeV"]),
        "scaled_hessian_eigenvalues":eig.tolist(),"corrected_gradient":g.tolist(),
        "old_uncorrected_composition_gradient_at_new_solution":float(oldg[2]),
        "required_linear_composition_gradient":float(DIFF*z[0]),
        "original_Ye_at_same_density":original["composition"]["Ye"],"restored_Ye":state["Y_e"],
        "gas_proton_insertion_margin_MeV":state["gas_proton_insertion_margin_MeV"]}
    assert record["min_tangent_gap_MeV_fm3"]>=-2e-9
    assert record["max_firstlaw_error_MeV"]<2e-5
    assert eig[0]>0 and np.max(abs(g))<1e-6
    assert abs(oldg[2]+DIFF*z[0])<1e-6
    assert abs(record["restored_Ye"]-record["original_Ye_at_same_density"])>1e-5
    record["status"]="PASS_SCOPED";out["models"][model]=record
    print(model,record["min_tangent_gap_MeV_fm3"],record["max_firstlaw_error_MeV"],eig[0],
        "deltaYe",record["restored_Ye"]-record["original_Ye_at_same_density"],flush=True)
out["experiment_sha256"]=hashlib.sha256((HERE/"rest_mass_seam_experiment.json").read_bytes()).hexdigest()
(HERE/"rest_mass_seam_validation.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
