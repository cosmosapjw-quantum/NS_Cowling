#!/usr/bin/env python3
"""Independent dimensional, seam, composition and local-curvature audit."""
from pathlib import Path
import sys,json,gzip,hashlib,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.api import read_compose
from compose_eos.candidate_gdfm_cld import reduced_cell,SurfaceParameters
from compose_eos.crust_assembly import rebase_header,select_rows,certify_interface
BASE=ROOT/'data/research_finish_20261116/crust_assembly'
def corrected_gradient(model,n,z,surface,rest_difference):
    gradient=reduced_cell(model,n,z,surface=surface)[1]
    return gradient+np.array([rest_difference*z[2],0.,rest_difference*z[0]])
def load(p):return json.loads(p.read_text())
def zipped(p):return json.loads(gzip.decompress(p.read_bytes()))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
out=[]
for kind in ['assembled_tf','restored_spg']:
 for directory in sorted((BASE/kind).glob('*')):
  if not (directory/'ASSEMBLY.json').exists():continue
  receipt=load(directory/'RECEIPT.json');manifest=load(directory/'ASSEMBLY.json')
  data=read_compose(directory/'compose')
  assert manifest['row_count']==data.n and np.all(np.diff(data.nb)>0)
  seams=[]
  for seam in manifest['interfaces']:
   i,j=seam['left_row_zero_based'],seam['right_row_zero_based']
   assert j==i+1 and seam['interpolation_across_density_gap']=='FORBIDDEN'
   assert data.nb[j]>data.nb[i]
   pe=abs(data.p[i]/data.p[j]-1);me=abs(data.mub[i]-data.mub[j])
   assert pe<3e-7 and me<2e-6,(directory,pe,me)
   seams.append(dict(rows=[i,j],n_left=data.nb[i],n_right=data.nb[j],relative_pressure_error=pe,mu_error_MeV=me))
   # Nonfinite thresholds must not bypass admission.
   left=select_rows(data,[i]);right=select_rows(data,[j])
   for key in ['pressure_rtol','mu_atol_MeV','euler_atol_MeV']:
    for bad in [float('nan'),float('inf')]:
     try:certify_interface(left,right,**{key:bad})
     except ValueError:pass
     else:raise AssertionError('invalid tolerance admitted')
  changed=rebase_header(data,mn=1000.,mp=999.)
  normal_errors={k:float(max(abs(getattr(data,k)-getattr(changed,k)))) for k in ['p','eps','mub','muq','mul','F']}
  assert max(normal_errors.values())<3e-10
  known=data.C3_applicable & data.C4_applicable
  baryon=float(max(abs(data.Yn[known]+data.Yp[known]+data.Yh[known]-1))) if known.any() else None
  charge=float(max(abs(data.Yp[known]+data.Yh[known]*data.Zav[known]/data.Aav[known]-data.Ye[known]-data.Ymu[known]))) if known.any() else None
  assert baryon<1e-12 and charge<1e-12
  row=dict(kind=kind,model=directory.name,rows=data.n,seams=seams,normalization_physical_errors=normal_errors,
   known_composition_rows=int(known.sum()),baryon_error=baryon,charge_error=charge,
   bsk_composition_mu_residual_range_MeV=receipt['bsk_composition_mu_residual_range_MeV'])
  if kind=='restored_spg':
   inner=zipped(directory/'inner_rows.json.gz');eigs=[];bounds=[];margins=[]
   for state in inner:
    n=state['n_B'];z=np.array(state['coordinates']);cols=[]
    surface=SurfaceParameters(**state['surface_parameters'])
    diff=state['rest_energies_MeV']['mp']-state['rest_energies_MeV']['mn']
    for j,h in enumerate([2e-5,2e-6,2e-6]):
     a=z.copy();b=z.copy();a[j]+=h;b[j]-=h
     cols.append((corrected_gradient(directory.name,n,a,surface,diff)-corrected_gradient(directory.name,n,b,surface,diff))/(2*h))
    H=np.array(cols).T;scale=np.array([1.,.1,.1])
    eigs.append(float(np.linalg.eigvalsh((H+H.T)*.5*scale[:,None]*scale[None,:])[0]))
    low=np.array([1e-7,max(.035,n*(1+1e-6)),.001]);high=np.array([1.,.25,.499])
    bounds.append(float(np.min(np.minimum(z-low,high-z)/(high-low))))
    margins.append(state['muon_margin_MeV'])
   assert min(eigs)>-1e-6 and min(bounds)>1e-6 and min(margins)>0
   row.update(inner_rows=len(inner),all_inner_minimum_scaled_hessian_eigenvalue=min(eigs),
    min_fractional_distance_from_artificial_bounds=min(bounds),min_muon_insertion_margin=min(margins),
    min_proton_insertion_margin=min(s['gas_proton_insertion_margin_MeV'] for s in inner))
  else:
   comp=load(directory/'COMPOSITION.json');profile=zipped(directory/'inner_source/selected_states.json.gz')
   assert comp['row_count']==data.n and comp['row_segments'][-1]['stop']==data.n
   idx=comp['row_segments'][-1]['start'];assert len(profile)==data.n-idx
   assert not data.Yn_present[idx:].any() and not data.Yp_present[idx:].any()
   assert data.Ye_present[idx:].all() and data.Ymu_present[idx:].all()
   errors=[]
   for i,s in enumerate(profile):
    w=np.array(s['profile']['weights']);p=s['profile']
    yn=w@np.array(p['nn'])/s['n_b'];yp=w@np.array(p['np'])/s['n_b'];ye=w@np.array(p['ne'])/s['n_b']
    errors.append(max(abs(yn+yp-1),abs(yp-ye),abs(ye-data.Ye[idx+i])))
   assert max(errors)<1e-9
   row.update(tf_profile_rows=len(profile),max_integrated_composition_error=max(errors))
  out.append(row)
result=dict(scope='Nine finite partial-crust tables; exact normalization and finite endpoint checks, not thermochemical consistency of BSk fits or a complete crust',
 tables=out,table_count=len(out),row_count=sum(r['rows'] for r in out),
 module_sha256={p:sha(ROOT/p) for p in ['python_frontend/compose_eos/crust_assembly.py','python_frontend/compose_eos/bsk22_extended_representation.py']},
 runner_sha256=sha(Path(__file__)))
assert len(out)==9
Path(__file__).with_name('CRUST_TABLE_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
