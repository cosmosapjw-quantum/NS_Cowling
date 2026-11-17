"""Bounded sampling check: exactly two refinements of one counted matrix cell."""
import json, os, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'scripts'))
import expand_supported_pairs_20261117 as m
os.environ.update(COMPOSE_EOS_FORTRAN_EXE='/home/cosmosapjw/ns_cowling_wt/core_support_20261117/build/compose_eos_fortran/compose_eos_fortran',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
BASE=ROOT/'data/tf_pair_expansion_20261117/matching'
r=json.load(open(BASE/'campaign_v1/RECEIPT.json'))
c=max((c for c in r['cells'] if c['method'] in m.SMOOTH),key=lambda c:c['result']['report']['consistency']['C8']['max_rel'])
sources={s['id']:s for s in r['sources']}
for side in ('crust','core'):m.revalidate(ROOT,sources[c[f'{side}_source_id']])
results={40:c['result']};paths={40:BASE/'campaign_v1/pairs'/c['id']}
for n in (80,160):
 cell=dict(c,n_bridge=n);dest=BASE/'sampling_refinement_v1'/f'{c["id"]}__n{n}'
 result=m.historical.run_smooth(ROOT/sources[c['crust_source_id']]['compose_path'],ROOT/sources[c['core_source_id']]['compose_path'],dest,cell)
 assert result['status']=='CONSTRUCTED_VERIFIED',result
 m.historical.save(dest/'RESULT.json',result);results[n]=result;paths[n]=dest
summary=dict(schema='ns-cowling.bridge-sampling-sensitivity.v1',cell_id=c['id'],new_matrix_coverage_cells=0,new_native_sampling_runs=2,
 independent_method_admission=False,cohort_receipt_unchanged=True,runs={},common_retained_nodes={},common_polynomial_windows=[])
arrays={}
for n,result in results.items():
 report=result['report'];v=report['consistency'];a=np.genfromtxt(paths[n]/'work/readback_internal/eos_data.csv',delimiter=',',names=True)
 arrays[n]=a
 summary['runs'][n]=dict(n_bridge=n,status=result['status'],rows=result['rows'],C7=v['C7']['max_rel'],C8=v['C8']['max_rel'],C9=v['C9']['max_rel'],first_law_relative_residual=report['first_law_relative_residual'],
 endpoints={k:report[k] for k in ('i1','i2','n1','n2','P1','P2','eps1','eps2','mu1','mu2')},spline_mode=report['spline_mode'],native_consistency=result['verification']['consistency_native']['overall'])
 # Reproduce native nonuniform three-point C8 diagnostic and categorical mask.
 x=a['nb'];p=x*a['Q1'];f=939.*(1+a['Q7']);d=np.empty(len(x))
 for i in range(1,len(x)-1):
  h1=x[i]-x[i-1];h2=x[i+1]-x[i]
  d[i]=-h2/(h1*(h1+h2))*f[i-1]+(h2-h1)/(h1*h2)*f[i]+h1/(h2*(h1+h2))*f[i+1]
 mask=np.ones(len(x),dtype=bool);mask[[0,-1]]=False
 for i in range(len(x)-1):
  if a['Iphase'][i]!=a['Iphase'][i+1] or np.log(x[i+1]/x[i])>.25:mask[i:i+2]=False
 residual=np.abs(p-x*x*d)/np.maximum(np.abs(p),np.abs(x*x*d));residual[~mask]=-1;i=int(np.argmax(residual))
 summary['runs'][n]['reproduced_C8_max']=float(residual[i]);summary['runs'][n]['C8_peak']=dict(row_zero_based=i,n_b_fm3=float(x[i]),u_log=float(np.log(x[i]/report['n1'])/np.log(report['n2']/report['n1'])),phase_id=int(a['Iphase'][i]),stencil_density_fm3=list(map(float,x[i-1:i+2])))
 assert abs(residual[i]-v['C8']['max_rel'])<1e-9
base=arrays[40]
for n in (80,160):
 a=arrays[n];common,ia,ib=np.intersect1d(base['nb'],a['nb'],return_indices=True)
 P0=base['nb'][ia]*base['Q1'][ia];P1=a['nb'][ib]*a['Q1'][ib];E0=base['nb'][ia]*939.*(1+base['Q7'][ia]);E1=a['nb'][ib]*939.*(1+a['Q7'][ib])
 summary['common_retained_nodes'][n]=dict(nodes=len(common),pressure_max_abs=float(np.max(abs(P1-P0))),epsilon_max_abs=float(np.max(abs(E1-E0))),all_shared_columns_exact=all(np.array_equal(base[name][ia],a[name][ib]) for name in base.dtype.names))
# On windows that do not cross spline knots, q=P/n is cubic and F=eps/n
# quartic in u=log(n/n1)/log(n2/n1). Fit each sampling independently, reject
# nonpolynomial windows, compare physical P/eps on the same nine u points.
report=results[40]['report'];dx=np.log(report['n2']/report['n1'])
for lo,hi in ((.05,.25),(.25,.45),(.45,.65),(.65,.85)):
 fits={};errors=[]
 for n,a in arrays.items():
  u=np.log(a['nb']/report['n1'])/dx;sel=(u>=lo)&(u<=hi)
  if sel.sum()<5:break
  q=a['Q1'][sel];f=939.*(1+a['Q7'][sel])-report['eps1']/report['n1']
  cq=np.polynomial.polynomial.polyfit(u[sel],q,3);cf=np.polynomial.polynomial.polyfit(u[sel],f,4)
  errors.append(max(np.max(abs(np.polynomial.polynomial.polyval(u[sel],cq)-q)),np.max(abs(np.polynomial.polynomial.polyval(u[sel],cf)-f))))
  fits[n]=(cq,cf)
 if len(fits)!=3 or max(errors)>1e-9:continue
 u=np.linspace(lo,hi,9);nb=report['n1']*np.exp(dx*u);evals={n:(nb*np.polynomial.polynomial.polyval(u,cq),nb*(np.polynomial.polynomial.polyval(u,cf)+report['eps1']/report['n1'])) for n,(cq,cf) in fits.items()}
 summary['common_polynomial_windows'].append(dict(u_interval=[lo,hi],max_fit_absolute_residual=float(max(errors)),common_physical_nodes=9,pressure_max_abs_difference=float(max(np.max(abs(evals[n][0]-evals[40][0])) for n in (80,160))),epsilon_max_abs_difference=float(max(np.max(abs(evals[n][1]-evals[40][1])) for n in (80,160)))))
assert summary['common_polynomial_windows']
summary['endpoints_exact']=all(summary['runs'][n]['endpoints']==summary['runs'][40]['endpoints'] for n in (80,160))
summary['interpretation']='one fixed source-bound analytic bridge, sampling changed only; native C7-C9 are three-point finite-difference diagnostics, not algebraic first-law identities; finite refinement is not a continuum convergence proof'
summary['status']='COMPLETE_BOUNDED_SAMPLING_CHECK'
m.historical.save(BASE/'sampling_refinement_v1/REFINEMENT.json',summary)
print(json.dumps(summary,indent=2))
