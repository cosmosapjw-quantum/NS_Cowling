from pathlib import Path
import dataclasses, json
import numpy as np
from compose_eos import read_compose
import compose_eos.skyrme_uniform as sk
from compose_eos.cold_leptons import free_fermion
import compose_eos.cold_leptons as cold
ROOT=Path.cwd()
out=ROOT/"build/skyrme_research/pressure_diagnosis"
out.mkdir(exist_ok=True)
results=[]
for model in sk.PARAMETERS:
    source=read_compose(ROOT/"eos_generator/eos_data"/model)
    mask=source.Iphase==3
    ns=source.nb[mask]; ps=source.p[mask]; es=source.eps[mask]
    table,states=sk.compose_table(model,ns)
    diff=(table.p-ps)/ps
    def runs(flag):
        selected=[]
        for i,ok in enumerate(flag):
            if ok:
                if selected and selected[-1]["end_index"]==i-1:
                    selected[-1].update(end_index=i,n_max_fm3=float(ns[i]))
                else:selected.append(dict(start_index=i,end_index=i,n_min_fm3=float(ns[i]),n_max_fm3=float(ns[i])))
        return selected
    worst=int(np.argmax(np.abs(diff)))
    indexes=sorted(set([0,worst,len(ns)//2,len(ns)-1]))
    pointcases=[]
    for i in indexes:
        n=ns[i]
        original=states[i]["pressure_mev_fm3"]
        variants={}
        oldk=sk.KINETIC_COEFFICIENT
        sk.KINETIC_COEFFICIENT=20.734
        variants["kinetic_coefficient_20_734"]=sk.beta_equilibrium(model,n)["pressure_mev_fm3"]
        sk.KINETIC_COEFFICIENT=oldk
        mn,mp=sk.MN,sk.MP
        sk.MN=sk.MP=939.0
        variants["common_rest_mass_939"]=sk.beta_equilibrium(model,n)["pressure_mev_fm3"]
        sk.MN,sk.MP=mn,mp
        variants["omit_muons"]=sk.beta_equilibrium(model,n,include_muons=False)["pressure_mev_fm3"]
        oldh=cold.HBARC
        cold.HBARC=197.32891
        variants["lepton_hbarc_197_32891"]=sk.beta_equilibrium(model,n)["pressure_mev_fm3"]
        cold.HBARC=oldh
        oldmu=sk.MMU
        sk.MMU=105.658
        variants["rounded_muon_mass_105_658"]=sk.beta_equilibrium(model,n)["pressure_mev_fm3"]
        sk.MMU=oldmu
        p=sk.PARAMETERS[model]
        if model=="RG_SLY4":
            sk.PARAMETERS[model]=dataclasses.replace(p,t0=-2488.913,t1=486.818,t2=-546.395)
            variants["higher_precision_SLY4_coefficients_sensitivity_only"]=sk.beta_equilibrium(model,n)["pressure_mev_fm3"]
            sk.PARAMETERS[model]=p
        rawidx=np.flatnonzero(mask)[i]
        nuclear=sk.nuclear_state(model,n*source.Yn[rawidx],n*source.Yp[rawidx])
        variants["fixed_source_composition"]=nuclear["pressure_mev_fm3"]+free_fermion(n*source.Ye[rawidx],sk.ME)[1]+free_fermion(n*source.Ymu[rawidx],sk.MMU)[1]
        pointcases.append({"n_fm3":float(n),"source_pressure":float(ps[i]),"baseline_pressure":original,
                          "baseline_relative_difference":float(diff[i]),
                          "variants":[{"name":key,"pressure":value,"relative_difference_to_source":(value-ps[i])/ps[i],"relative_change_from_baseline":(value-original)/original} for key,value in variants.items()]})
    p_from_source=np.gradient(es/ns,ns,edge_order=2)*ns**2
    p_from_model_grid=np.gradient(table.eps/ns,ns,edge_order=2)*ns**2
    result={"model_id":model,"source_homogeneous_support_fm3":[float(ns[0]),float(ns[-1])],
            "pressure_within_1percent_sampled_ranges":runs(np.abs(diff)<.01),
            "pressure_within_1percent_rows":int(np.count_nonzero(np.abs(diff)<.01)),"rows":len(ns),
            "point_sensitivities":pointcases,
            "source_energy_grid_derivative_pressure_max_rel_difference":float(np.max(np.abs((p_from_source-ps)/ps))),
            "model_energy_same_grid_derivative_pressure_max_rel_difference":float(np.max(np.abs((p_from_model_grid-table.p)/table.p))),
            "derivative_pressure_first_four":[{"n":float(ns[i]),"source_P":float(ps[i]),"source_E_grid_derivative_P":float(p_from_source[i]),"model_P":float(table.p[i]),"model_E_grid_derivative_P":float(p_from_model_grid[i])} for i in range(4)],
            "attribution":"Sensitivity calculations only; provider pressure-construction method remains unresolved."}
    results.append(result)
    print(json.dumps(result,indent=2))
(out/"diagnosis.json").write_text(json.dumps(results,indent=2)+"\n")
