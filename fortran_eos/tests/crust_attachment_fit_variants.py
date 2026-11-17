#!/usr/bin/env python3
"""Quantify released-author-code vs printed-paper chemical fit; no production changes."""
import os,sys,pathlib,json,math
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.bsk22_analytic import state
p=ROOT/'data/research_blockers_20261116/lower'
rows=np.loadtxt(p/'reference_integral_rows.csv',delimiter=',');out=[]
for row in rows:
 n=row[0];s=state(n,thermodynamic_mode='source_fit');nf=n*s['Y_n_free'];mn=939.5654
 code=mn+.211*mn*nf**(2/3)/(1+9.41*math.sqrt(nf)+7.62*nf-259*nf*nf)
 paper=mn+197.3269804**2*(3*math.pi**2*nf)**(2/3)/(2*mn)/(1+9.41*math.sqrt(n)+7.62*n+259*n*n)
 out.append({'n':n,'n_free':nf,'Y_n_free':s['Y_n_free'],'mu_author_code':code,'mu_printed_paper':paper,'difference_MeV':code-paper,'source_euler':s['mu_b_euler'],'integrated_euler':row[3]})
res={'definition':'Same source neutron rest mass939.5654MeV; printed C18/HTML92 uses n,+p3n²; released code INCRST uses n_free,-p3n_free². HBARC197.3269804 for paper kinetic numerator.','max_abs_code_paper_mu_MeV':max(abs(r['difference_MeV']) for r in out),'max_abs_paper_mu_minus_source_euler_MeV':max(abs(r['mu_printed_paper']-r['source_euler']) for r in out),'max_abs_code_mu_minus_source_euler_MeV':max(abs(r['mu_author_code']-r['source_euler']) for r in out),'rows':out}
(p/'CHEMICAL_FIT_VARIANTS.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({k:v for k,v in res.items() if k!='rows'},indent=2))
