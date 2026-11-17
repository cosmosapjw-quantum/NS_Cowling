from pathlib import Path
import json
from compose_eos.skyrme_uniform import PARAMETERS,beta_equilibrium,nuclear_state,sound_speed_ratio,ME,MMU
from compose_eos.cold_leptons import free_fermion
rows=[]
for model in PARAMETERS:
    for n in [.08,.16,.5,1.0,1.4]:
        eq=beta_equilibrium(model,n)
        xp,xe,xmu=eq["Y_p"],eq["Y_e"],eq["Y_mu"]
        def frozen_pressure(density):
            matter=nuclear_state(model,density*(1-xp),density*xp)
            return matter["pressure_mev_fm3"]+free_fermion(density*xe,ME)[1]+free_fermion(density*xmu,MMU)[1]
        h=n*1e-4
        frozen=(frozen_pressure(n+h)-frozen_pressure(n-h))/(2*h*eq["mu_n_mev"])
        equilibrium=sound_speed_ratio(model,n)
        row={"model_id":model,"n_fm3":n,"equilibrium_cs2_over_c2":equilibrium,"frozen_composition_cs2_over_c2":frozen,
             "relaxation_difference":frozen-equilibrium}
        assert frozen>=equilibrium-1e-7
        rows.append(row)
        print(json.dumps(row))
Path("build/skyrme_research/frozen_composition.json").write_text(json.dumps(rows,indent=2)+"\n")
