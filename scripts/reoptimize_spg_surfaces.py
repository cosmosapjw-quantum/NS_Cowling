#!/usr/bin/env python3
"""Re-equilibrate five SPG candidates under explicit p=3 surface prescriptions."""
from __future__ import annotations
# Historical research oracle; production physics lives in the Fortran backend.
# Imports remain available to pinned-campaign audit and regression helpers.
import sys as _oracle_sys
if __name__ == '__main__':
    if '--historical-python-oracle' not in _oracle_sys.argv:
        raise SystemExit('This is a historical Python research oracle. Use compose-eos generate for production; add --historical-python-oracle only to reproduce the explicitly archived Python campaign.')
    _oracle_sys.argv.remove('--historical-python-oracle')
import argparse,concurrent.futures,gzip,hashlib,json,os,sys,time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.candidate_gdfm_cld import SurfaceParameters,default_surface_parameters,solve_two_phase
from compose_eos.spg_crust_generator import neutron_drip,homogeneous_endpoint,sequence,to_eos_data
from compose_eos.gdfm_uniform import beta_state
from compose_eos.physical_domain import physical_domain_report
from compose_eos.api import write_compose,read_compose,check_consistency
from generate_spg_crust import BRACKETS,compare_grid

VARIANT_PATH=ROOT/'data/research_completion_20261116/theory_review/surface_fit_variants_321.json'
BASE_PATH=ROOT/'data/research_generator_20261116/candidate_surface/candidate_surface_fit.json'
LABELS=('AME2408_candidate','geometric_area__nine','geometric_area__eight_measured',
        'fixed_symmetric_area__nine','fixed_symmetric_area__eight_measured')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write_json(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def prescriptions(model):
    raw=VARIANT_PATH.read_bytes();data=json.loads(raw)['models'][model]['variants']
    return {'AME2408_candidate':default_surface_parameters(model),
        **{label:SurfaceParameters(model,float(fit['sigma0_MeV_fm2']),float(fit['b']),
            label,hashlib.sha256(raw).hexdigest()) for label,fit in data.items()}}

def resolve_endpoint(model,surface):
    attempts=[]
    base=BRACKETS[model]
    for bracket in (base,(max(.025,base[0]-.01),min(.119,base[1]+.01))):
        try:
            result=homogeneous_endpoint(model,bracket,surface=surface)
            attempts.append({'bracket':bracket,'status':'CONVERGED'})
            return result,attempts
        except (ValueError,RuntimeError) as e:
            attempts.append({'bracket':bracket,'status':'FAILED','error':str(e)})
    # Retain the actual stationary crossing survey instead of extrapolating.
    survey=[]
    for n in np.linspace(max(.025,base[0]-.015),min(.119,base[1]+.02),18):
        try:
            cell,tries=solve_two_phase(model,float(n),surface=surface)
            survey.append({'n_B':float(n),'status':cell['status'],
                'delta_epsilon':cell['energy']-beta_state(model,float(n))['epsilon'],
                'max_projected_gradient':cell['max_projected_gradient']})
        except (ValueError,RuntimeError) as e:survey.append({'n_B':float(n),'status':'ERROR','error':str(e)})
    attempts.append({'stationary_crossing_survey':survey})
    for left,right in zip(survey,survey[1:]):
        if (left['status']==right['status']=='STATIONARY_CANDIDATE' and
            left['delta_epsilon']*right['delta_epsilon']<0):
            bracket=(left['n_B'],right['n_B'])
            try:
                result=homogeneous_endpoint(model,bracket,surface=surface)
                attempts.append({'bracket':bracket,'status':'SURVEY_CONVERGED'})
                return result,attempts
            except (ValueError,RuntimeError) as e:attempts.append({'bracket':bracket,'error':str(e)})
    raise RuntimeError(json.dumps({'endpoint_attempts':attempts}))

def export_native(out,rows):
    data=to_eos_data(rows);write_compose(out,data);restored=read_compose(out)
    errors={}
    for f in ('nb','Q1','Q2','Q3','Q4','Q5','Q6','Q7','Yn','Yp','Ye','Ymu','Yh','Aav','Zav'):
        errors[f]=float(np.max(abs(getattr(data,f)-getattr(restored,f))))
        if not np.allclose(getattr(data,f),getattr(restored,f),rtol=3e-14,atol=3e-15):
            raise RuntimeError('native readback mismatch '+f)
    if not np.array_equal(data.Iphase,restored.Iphase):raise RuntimeError('native phase mismatch')
    for f in ('Yn_present','Yp_present','Ye_present','Ymu_present','Yh_present',
              'Abar_present','Zbar_present','C3_applicable','C4_applicable'):
        if not getattr(restored,f).all():raise RuntimeError('native composition mask '+f)
    checks=check_consistency(restored,verbose=False)
    # Functional/constraint implementation failures differ from model support.
    # Acausality/compressibility findings stay density-resolved diagnostics.
    for key in ('C1','C2','C3','C4','C5'):
        if not checks[key]['pass_']:raise RuntimeError('implementation constraint '+key)
    return {'rows':len(rows),'native_readback_max_abs':errors,'native_checks':checks,
        'physical_support':physical_domain_report(restored),
        'density_interval_fm3':[float(data.nb[0]),float(data.nb[-1])],
        'composition_ranges':{f:[float(np.min(getattr(data,f))),float(np.max(getattr(data,f)))]
                            for f in ('Yn','Yp','Ye','Yh','Aav','Zav')}}


def case(task):
    model,label,points,outroot,endpoint_only=task
    started=time.monotonic();out=Path(outroot)/model/label;out.mkdir(parents=True)
    p=prescriptions(model)[label]
    record={'model':model,'label':label,'surface_parameters':asdict(p),'status':'RUNNING'}
    write_json(out/'RECEIPT.json',record)
    try:
        drip=neutron_drip(model,surface=p)
        endpoint,endpoint_attempts=resolve_endpoint(model,p)
        record.update(drip=drip,transition=endpoint,endpoint_search=endpoint_attempts)
        if endpoint_only:
            record.update(status='ENDPOINTS_RESOLVED',runtime_seconds=time.monotonic()-started)
            write_json(out/'RECEIPT.json',record);return record
        q=(points-1)//4;edge=.9*endpoint['n_cluster_fm3']
        grid=np.r_[np.geomspace(drip['n_B'],edge,3*q+1),
                   np.linspace(edge,endpoint['n_cluster_fm3'],q+1)[1:]]
        fine=sequence(model,grid,surface=p)
        failures=[{'n_B':r['n_B'],'violations':r['best']['violations']} for r in fine
                  if r['best']['status']!='STATIONARY_CANDIDATE']
        payload={'fine':fine,'drip':drip,'transition':endpoint,'surface_parameters':asdict(p)}
        with gzip.open(out/'all_attempts.json.gz','wt') as stream:json.dump(payload,stream,allow_nan=False)
        record['failures']=failures
        # Preserve and export consecutive stationary runs; never bridge HOLDs.
        blocks=[];block=[]
        for row in fine:
            if row['best']['status']=='STATIONARY_CANDIDATE':block.append(row)
            elif block:blocks.append(block);block=[]
        if block:blocks.append(block)
        exports=[]
        for i,block in enumerate(blocks):
            if len(block)<3:
                exports.append({'rows':len(block),'status':'TOO_SHORT_FOR_NATIVE_TABLE',
                                'density_interval_fm3':[block[0]['n_B'],block[-1]['n_B']]})
                continue
            name='compose' if not failures else f'compose_segment_{i:03d}'
            exports.append({'path':name,'status':'NATIVE_READBACK_VERIFIED',**export_native(out/name,block)})
        record.update(status='REOPTIMIZED_NATIVE_COMPOSITION_VERIFIED' if not failures else 'PARTIAL_STATIONARY_SEGMENTS',
            rows=len(fine),stationary_rows=sum(len(b) for b in blocks),native_exports=exports,
            max_projected_gradient=max(r['best']['max_projected_gradient'] for r in fine),
            positive_gas_rows=sum(r['best'].get('gas_protons')=='VARIATIONAL_POSITIVE' for r in fine))
        if not failures:
            record['density_refinement']={'quarter_to_half':compare_grid(fine[::4],fine[::2]),
                'half_to_fine':compare_grid(fine[::2],fine),'row_counts':[len(fine[::4]),len(fine[::2]),len(fine)],
                'scope':'Nested interpolation accuracy using independently optimized fine-grid states; not separate optimizer-continuation reruns.'}
        record['output_sha256']={str(x.relative_to(out)):sha(x) for x in out.glob('compose*/eos.*')}
    except Exception as e:
        record.update(status='HOLD',error_type=type(e).__name__,error=str(e))
    record['runtime_seconds']=time.monotonic()-started
    write_json(out/'RECEIPT.json',record)
    return record

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    parser.add_argument('--points',type=int,default=129);parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--models',nargs='+',choices=BRACKETS,default=list(BRACKETS))
    parser.add_argument('--variants',nargs='+',choices=LABELS,default=list(LABELS))
    parser.add_argument('--endpoint-only',action='store_true');a=parser.parse_args()
    if a.points<17 or (a.points-1)%8:parser.error('points must be8k+1, at least17')
    a.output.mkdir(parents=True,exist_ok=False)
    sources=[Path(__file__),ROOT/'scripts/generate_spg_crust.py',
        *[ROOT/'python_frontend/compose_eos'/x for x in ('candidate_gdfm_cld.py','spg_crust_generator.py',
            'gdfm_uniform.py','scurto_interface.py','cold_leptons.py','eos_data.py','api.py','physical_domain.py')]]
    hashes={str(p.relative_to(ROOT)):sha(p) for p in sources}
    summary={'status':'RUNNING','base_commit':'573a31485bc9af30915e01fe74f259e4fa5ac54f',
        'claim':'Reoptimized p=3 candidate prescription ensemble; not original-fit recovery or a statistical uncertainty band',
        'source_sha256':hashes,'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in (BASE_PATH,VARIANT_PATH)},
        'points':a.points,'cases':[]}
    write_json(a.output/'SUMMARY.json',summary)
    tasks=[(m,v,a.points,str(a.output),a.endpoint_only) for m in a.models for v in a.variants]
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
        for result in pool.map(case,tasks):
            summary['cases'].append({k:result[k] for k in ('model','label','status','runtime_seconds')})
            write_json(a.output/'SUMMARY.json',summary)
            print(json.dumps(summary['cases'][-1]),flush=True)
    summary['source_unchanged_during_execution']=hashes=={str(p.relative_to(ROOT)):sha(p) for p in sources}
    summary['status']='COMPLETED_SCOPED' if all(r['status'] in ('REOPTIMIZED_NATIVE_COMPOSITION_VERIFIED','ENDPOINTS_RESOLVED') for r in summary['cases']) else 'PARTIAL_WITH_HOLD'
    write_json(a.output/'SUMMARY.json',summary)
    return 0 if summary['status']=='COMPLETED_SCOPED' and summary['source_unchanged_during_execution'] else 1
if __name__=='__main__':raise SystemExit(main())
