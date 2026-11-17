import sys,json
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'scripts'))
from continue_tf_20261117 import *
base=json.loads((OUT/'N256_d3_droplet_n1_warm_fine_r090_restart/STATE.json').read_text());rows=[]
h=1e-6
for k in [1,2,4]:
 s=run(base,256,base['n_b']+k*h,f'N256_d3_derivative_h{k}');rows.append(s)
 if not s.get('accepted_stationary') and s.get('optimizer_mode')==8:
  s=run(s,256,s['n_b'],s['id']+'_restart');rows.append(s)
 dump(OUT/'DERIVATIVE_STATES.json',rows)
values={k:next(s for s in rows[::-1] if s['id'].split('_restart')[0]==f'N256_d3_derivative_h{k}') for k in [1,2,4]}
e0=base['epsilon']/base['n_b'];e={k:s['epsilon']/s['n_b'] for k,s in values.items()}
p1=base['n_b']**2*(-3*e0+4*e[1]-e[2])/(2*h);p2=base['n_b']**2*(-3*e0+4*e[2]-e[4])/(4*h)
dump(OUT/'THERMODYNAMIC_CHECK.json',{'definition':'P=n_B^2 d(epsilon/n_B)/dn_B; all displaced profiles reoptimized native, same N256 and branch','base_id':base['id'],'h_fm3':h,'P_Euler':base['pressure'],'P_forward_h':p1,'P_forward_2h':p2,'h_difference':p1-p2,'difference_from_Euler_h':p1-base['pressure'],'relative_difference_from_Euler_h':(p1-base['pressure'])/base['pressure'],'all_displaced_stationary':all(s.get('accepted_stationary') for s in values.values()),'state_ids':[s['id'] for s in values.values()],'scope':'Finite difference consistency; residual/mesh error distinct from derivative truncation.'})
