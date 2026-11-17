#!/usr/bin/env python3
"""Direct common-density prescription response and independent continuation checks."""
from __future__ import annotations
import argparse,concurrent.futures,gzip,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.candidate_gdfm_cld import solve_two_phase
from compose_eos.spg_crust_generator import sequence
from reoptimize_spg_surfaces import prescriptions,LABELS,BRACKETS,sha,write_json
from generate_spg_crust import compare_grid

def read(p):return json.loads(Path(p).read_text())
def rows(p):return json.load(gzip.open(p/'all_attempts.json.gz','rt'))['fine']
def observables(row):
    b=row['best']
    return {'pressure_MeV_fm3':b['pressure_MeV_fm3'],'epsilon_MeV_fm3':b['energy'],
            'mu_B_MeV':b['mu_B_MeV'],**b['composition']}
def onset(states):
    for i,r in enumerate(states):
        if r['best'].get('gas_protons')=='VARIATIONAL_POSITIVE':
            return [states[max(0,i-1)]['n_B'],r['n_B']]
    return None

def calculate(task):
    model,label,grid,outroot=task
    p=prescriptions(model)[label];out=Path(outroot)/model/label;out.mkdir(parents=True)
    states=sequence(model,grid,surface=p)
    with gzip.open(out/'common_grid_attempts.json.gz','wt') as f:json.dump(states,f,allow_nan=False)
    checks=[]
    for i in (len(grid)//4,len(grid)//2,3*len(grid)//4,7*len(grid)//8,len(grid)-2):
        state=states[i]['best'];n=float(grid[i]);h=n*1e-5
        low,lowtries=solve_two_phase(model,n-h,surface=p)
        high,hightries=solve_two_phase(model,n+h,surface=p)
        mu=(high['energy']-low['energy'])/(2*h)
        checks.append({'n_B':n,'h':h,'center':state,'low':low,'high':high,
            'low_attempts':lowtries,'high_attempts':hightries,'mu_B_finite_difference':mu,
            'mu_B_analytic':state['mu_B_MeV'],'relative_mu_error':abs(mu/state['mu_B_MeV']-1),
            'absolute_pressure_error_MeV_fm3':abs(n*mu-state['energy']-state['pressure_MeV_fm3']),
            'stationary':all(r['status']=='STATIONARY_CANDIDATE' for r in (low,state,high))})
    with gzip.open(out/'first_law_attempts.json.gz','wt') as f:json.dump(checks,f,allow_nan=False)
    return {'model':model,'label':label,'states':states,'first_law':[
        {k:v for k,v in c.items() if k not in ('center','low','high','low_attempts','high_attempts')} for c in checks]}

def compare_states(a,b):
    metrics={}
    for key in ('pressure_MeV_fm3','epsilon_MeV_fm3','Ye','Yn','Yp','Yh','Aav','Zav'):
        x=np.array([observables(r)[key] for r in a]);y=np.array([observables(r)[key] for r in b])
        delta=y-x;mask=abs(x)>1e-14;idx=int(np.argmax(abs(delta)))
        metric={'max_absolute':float(max(abs(delta))),'n_at_max_absolute':a[idx]['n_B'],
                'signed_range':[float(min(delta)),float(max(delta))]}
        if mask.any():
            rel=np.zeros_like(x);rel[mask]=abs(delta[mask]/x[mask]);idx=int(np.argmax(rel))
            metric.update(max_relative_nonzero_baseline=float(rel[idx]),n_at_max_relative=a[idx]['n_B'],
                          relative_denominator_min=1e-14)
        metrics[key]=metric
    return metrics

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('coarse',type=Path);ap.add_argument('fine',type=Path)
    ap.add_argument('output',type=Path);ap.add_argument('--points',type=int,default=129);ap.add_argument('--workers',type=int,default=4)
    a=ap.parse_args()
    if a.points<17 or (a.points-1)%8:ap.error('points must be8k+1, at least17')
    a.output.mkdir(parents=True,exist_ok=False)
    cs=read(a.coarse/'SUMMARY.json');fs=read(a.fine/'SUMMARY.json')
    if cs['status']!='COMPLETED_SCOPED' or fs['status']!='COMPLETED_SCOPED':raise RuntimeError('completed campaigns required')
    if cs['input_sha256']!=fs['input_sha256'] or cs['source_sha256']!=fs['source_sha256']:
        raise RuntimeError('input/source identity differs between grid campaigns')
    for name,digest in fs['source_sha256'].items():
        if sha(ROOT/name)!=digest:raise RuntimeError('live source changed: '+name)
    source_hashes={**fs['source_sha256'],str(Path(__file__).relative_to(ROOT)):sha(__file__)}
    result={'schema':'spg-reoptimized-prescription-comparison-v1','source_sha256':source_hashes,
            'input_sha256':fs['input_sha256'],'model_results':{},
            'claim':'Actual fixed-p=3 prescription response, not original calibration recovery or statistical confidence band',
            'units':{'n_B':'fm^-3','pressure':'MeV fm^-3','epsilon':'MeV fm^-3','Aav':'baryon count','Zav':'proton count'}}
    tasks=[]
    for model in BRACKETS:
        rec={label:read(a.fine/model/label/'RECEIPT.json') for label in LABELS}
        lower=max(r['drip']['n_B'] for r in rec.values());upper=min(r['transition']['n_cluster_fm3'] for r in rec.values())
        q=(a.points-1)//4;edge=.9*upper
        grid=np.r_[np.geomspace(lower,edge,3*q+1),np.linspace(edge,upper,q+1)[1:]].tolist()
        result['model_results'][model]={'common_support_fm3':[lower,upper],'common_points':len(grid),'variants':{}}
        tasks.extend((model,label,grid,str(a.output)) for label in LABELS)
        base=rec[LABELS[0]]
        for label,r in rec.items():
            coarse=rows(a.coarse/model/label);fine=rows(a.fine/model/label)
            if not np.allclose([x['n_B'] for x in coarse],[x['n_B'] for x in fine[::2]],rtol=2e-14,atol=0):
                raise RuntimeError('grid not nested '+model+label)
            result['model_results'][model]['variants'][label]={
                'surface_parameters':r['surface_parameters'],'drip_n_B':r['drip']['n_B'],
                'drip_shift_vs_baseline':r['drip']['n_B']-base['drip']['n_B'],
                'transition':{k:r['transition'][k] for k in ('n_cluster_fm3','n_uniform_fm3','pressure_MeV_fm3','mu_B_MeV')},
                'transition_cluster_shift_vs_baseline':r['transition']['n_cluster_fm3']-base['transition']['n_cluster_fm3'],
                'gas_proton_onset_bracket_fm3':onset(fine),
                'independent_coarse_to_fine_interpolation':compare_grid(coarse,fine),
                'same_density_independent_continuation':compare_states(coarse,fine[::2]),
                'native_physical_checks':{k:r['native_exports'][0]['native_checks'][k] for k in ('C6','C10','C11')},
                'max_projected_gradient':r['max_projected_gradient']}
    actual={}
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
        for data in pool.map(calculate,tasks):
            model=data['model'];label=data['label'];actual[(model,label)]=data['states']
            result['model_results'][model]['variants'][label]['first_law']=data['first_law']
            print(model,label,'common_grid_complete',flush=True)
    for model in BRACKETS:
        baseline=actual[(model,LABELS[0])]
        for label in LABELS:
            states=actual[(model,label)]
            result['model_results'][model]['variants'][label]['common_density_response_vs_baseline']=compare_states(baseline,states)
            result['model_results'][model]['variants'][label]['common_stationary_rows']=sum(
                r['best']['status']=='STATIONARY_CANDIDATE' for r in states)
    result['source_unchanged_during_execution']=source_hashes=={name:sha(ROOT/name) for name in source_hashes}
    result['all_common_rows_stationary']=all(r['best']['status']=='STATIONARY_CANDIDATE' for s in actual.values() for r in s)
    result['first_law_max_relative_mu_error']=max(c['relative_mu_error'] for m in result['model_results'].values()
                                               for v in m['variants'].values() for c in v['first_law'])
    result['all_first_law_stencils_stationary']=all(c['stationary'] for m in result['model_results'].values()
                                                   for v in m['variants'].values() for c in v['first_law'])
    result['first_law_relative_mu_tolerance']=1e-7
    result['first_law_pass']=(result['all_first_law_stencils_stationary'] and
        result['first_law_max_relative_mu_error']<=result['first_law_relative_mu_tolerance'])
    result['output_sha256']={str(p.relative_to(a.output)):sha(p) for p in a.output.glob('*/*/*.gz')}
    write_json(a.output/'COMPARISON.json',result)
    return 0 if (result['all_common_rows_stationary'] and result['source_unchanged_during_execution']
                 and result['first_law_pass']) else 1
if __name__=='__main__':raise SystemExit(main())
