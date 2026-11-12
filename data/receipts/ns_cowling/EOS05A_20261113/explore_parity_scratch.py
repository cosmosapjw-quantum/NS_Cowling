import sys, json, time
from pathlib import Path
import numpy as np
W=Path('/home/cosmosapjw/ns_cowling_wt/eos05')
from compose_eos import read_compose, match_crust_core
CR=['SPG_M1','SPG_M2','SPG_M3','SPG_M4','SPG_M5','GPPVA_DD2','GPPVA_DDME2','GPPVA_FSU2','GPPVA_FSU2H','GPPVA_FSU2R','GPPVA_NL3wrL55','GPPVA_TM1e','GPPVA_TW']
CO=['BBB_BHF2','GM_GM1','RG_SLY4']
M=sys.argv[1:] or ['cubic_F','cubic_P','quintic_F','quintic_P']
N=W/'eos_generator/eos_data_normalized'; F=W/'data/eos_original_normalized'
COLS={'nb':'n_b_fm3','Q1':'q1_p_over_nb_mev','Q3':'q3_mu_b_over_mn_minus_1','Q4':'q4_mu_q_over_mn','Q6':'q6_f_over_nb_mn_minus_1','Q7':'q7_e_over_nb_mn_minus_1'}
cache={}
def rd(n):
    if n not in cache: cache[n]=read_compose(N/n)
    return cache[n]
rows=[];t0=time.time()
for co in CO:
  for cr in CR:
    for m in M:
      key=f'{cr}__{co}__{m}'
      try:
        d,info=match_crust_core(rd(cr),rd(co),cr,co,method=m,n_bridge=40,drop_composition=True)
      except Exception as e:
        rows.append({'case':key,'native':'REJECT','msg':str(e)[:160]});continue
      fz=np.genfromtxt(F/f'{key}.csv',delimiter=',',names=True)
      r={'case':key,'native':'OK','n_native':d.n,'n_frozen':len(fz),'spline_mode':info['spline_mode'],'i1':info['i1'],'i2':info['i2'],'fl':info['first_law_relative_residual']}
      if d.n==len(fz):
        i1,i2=info['i1'],info['i2']; nb_=40
        for k,c in COLS.items():
            a=getattr(d,k);b=fz[c]
            rel=np.abs(a-b)/np.maximum(np.abs(b),1e-300)
            ok=np.isclose(a,b,rtol=8e-7,atol=2e-9)
            r[k+'_bad']=int((~ok).sum()); r[k+'_maxrel']=float(rel.max())
            if k=='Q4':
                br=slice(i1+1,i1+1+nb_); mask=np.ones(len(a),bool); mask[br]=False
                r['Q4_nonbridge_bad']=int((~ok[mask]).sum()); r['Q4_bridge_bad']=int((~ok[br]).sum())
      rows.append(r)
json.dump(rows,open('/home/cosmosapjw/ns_cowling_wt/scratch/explore_%s.json'%('_'.join(M)),'w'),indent=1)
print('cases',len(rows),'secs',round(time.time()-t0,1))
rej=[r for r in rows if r['native']!='OK']; print('rejected',len(rej)); [print(' ',r) for r in rej[:5]]
ok=[r for r in rows if r['native']=='OK']
print('size mismatch',sum(r['n_native']!=r['n_frozen'] for r in ok))
for k in COLS: print(k,'cases with bad pts',sum(1 for r in ok if r.get(k+'_bad',0)>0),'worst maxrel',max((r.get(k+'_maxrel',0) for r in ok)))
print('Q4 nonbridge bad cases',sum(1 for r in ok if r.get('Q4_nonbridge_bad',0)>0),'bridge bad cases',sum(1 for r in ok if r.get('Q4_bridge_bad',0)>0))
from collections import Counter
print(Counter(r['spline_mode'] for r in ok))
