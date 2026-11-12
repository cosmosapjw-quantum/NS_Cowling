import sys, json
sys.path.insert(0,'python_frontend')
import numpy as np, types, os
from pathlib import Path
from compose_eos import read_compose, match_crust_core
CR=['SPG_M1','SPG_M2','SPG_M3','SPG_M4','SPG_M5','GPPVA_DD2','GPPVA_DDME2','GPPVA_FSU2','GPPVA_FSU2H','GPPVA_FSU2R','GPPVA_NL3wrL55','GPPVA_TM1e','GPPVA_TW']
N=Path('eos_generator/eos_data_normalized'); F=Path('data/eos_original_normalized')
COLS={'nb':'n_b_fm3','Q1':'q1_p_over_nb_mev','Q3':'q3_mu_b_over_mn_minus_1','Q4':'q4_mu_q_over_mn','Q6':'q6_f_over_nb_mn_minus_1','Q7':'q7_e_over_nb_mn_minus_1'}
out=[]
for co in ['BBB_BHF2','GM_GM1']:
  for cr in CR:
    d,i=match_crust_core(read_compose(N/cr),read_compose(N/co),cr,co,method='legacy_npphb',n_bridge=40,drop_composition=True)
    fz=np.genfromtxt(F/f'{cr}__{co}__maxwell.csv',delimiter=',',names=True)
    mr=max(float(np.max(np.abs(getattr(d,k)-fz[c])/np.maximum(np.abs(fz[c]),1e-300))) for k,c in COLS.items())
    out.append({'case':f'{cr}__{co}__maxwell','status':'PARITY_MATCH' if all(np.allclose(getattr(d,k),fz[c],rtol=8e-7,atol=2e-9) for k,c in COLS.items()) else 'PARITY_MISMATCH','rows':d.n,'max_rel_diff':mr,'first_law_residual_same_rule':i['first_law_relative_residual'],'cs2_min_bridge':i['cs2_min_bridge'],'cs2_max_bridge':i['cs2_max_bridge'],'pressure_span_bridge':i['pressure_span_bridge'],'dP':i['dP']})
json.dump(out,open('/home/cosmosapjw/ns_cowling_wt/scratch/explore_npphb.json','w'),indent=1)
print(len(out),'worst max_rel',max(o['max_rel_diff'] for o in out),'statuses',{o['status'] for o in out})
print('pressure_span/dP range',min(o['pressure_span_bridge']/o['dP'] for o in out),max(o['pressure_span_bridge']/o['dP'] for o in out))
