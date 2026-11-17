#!/usr/bin/env python3
"""Generate candidate SPG inner crusts on refined grids with native composition."""
from __future__ import annotations
import argparse,gzip,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.spg_crust_generator import neutron_drip,homogeneous_endpoint,sequence,to_eos_data
from compose_eos.api import write_compose,read_compose,check_consistency
from compose_eos.physical_domain import physical_domain_report
from compose_eos.candidate_gdfm_cld import MASS, solve_two_phase
BRACKETS={'SPG_M1':(.07,.075),'SPG_M2':(.085,.09),'SPG_M3':(.075,.08),'SPG_M4':(.09,.1),'SPG_M5':(.07,.075)}


def compare_grid(coarse,fine):
    nc=np.array([r['n_B'] for r in coarse]);nf=np.array([r['n_B'] for r in fine])
    result={}
    for name in ['pressure_MeV_fm3','energy','mu_B_MeV']:
        vc=np.array([r['best'][name] for r in coarse]);vf=np.array([r['best'][name] for r in fine])
        pred=PchipInterpolator(nc,vc)(nf);delta=pred-vf
        result[name]={'max_absolute':float(max(abs(delta))),
            'max_relative':float(max(abs(delta)/np.maximum(abs(vf),1e-20)))}
    for name in ['Ye','Yn','Yp','Yh','Aav','Zav']:
        vc=np.array([r['best']['composition'][name] for r in coarse]);vf=np.array([r['best']['composition'][name] for r in fine])
        result[name]={'max_absolute':float(max(abs(PchipInterpolator(nc,vc)(nf)-vf)))}
    return result


def write_json(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output',type=Path);ap.add_argument('--points',type=int,default=257)
    ap.add_argument('--models',nargs='+',choices=list(BRACKETS),default=list(BRACKETS))
    args=ap.parse_args()
    if args.points<17 or (args.points-1)%8:ap.error('points must be 8k+1 and at least 17 for nested grids')
    args.output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();summary={'status':'RUNNING','models':{},'points':args.points,
        'comparison_policy':'CompOSE differences are diagnostics; implementation acceptance uses variational equations, conservation, numerical refinement and native readback',
        'scope':'Spherical two-phase GDFM CLD, calibrated candidate surface, electron gas; standalone inner crust from candidate neutron drip to same-functional coexistence endpoint',
        'source_calibrated':False,'global_minimum_status':'FINITE_STARTS_AND_LOCAL_CURVATURE_ONLY'}
    write_json(args.output/'SUMMARY.json',summary)
    for model in args.models:
        out=args.output/model;out.mkdir()
        drip=neutron_drip(model);transition=homogeneous_endpoint(model,BRACKETS[model])
        def grid(points):
            q=(points-1)//4;edge=.9*transition['n_cluster_fm3']
            return np.r_[np.geomspace(drip['n_B'],edge,3*q+1),np.linspace(edge,transition['n_cluster_fm3'],q+1)[1:]]
        coarse_n=grid((args.points+1)//2)
        fine_n=grid(args.points)
        coarse=sequence(model,coarse_n);fine=sequence(model,fine_n)
        violations=[{'n_B':r['n_B'],'violations':r['best']['violations']} for r in coarse+fine if r['best']['status']!='STATIONARY_CANDIDATE']
        positive_start_checks=[]
        for density in (.001,.01,.04,.95*transition['n_cluster_fm3']):
            candidate,attempts=solve_two_phase(model,density,force_proton_search=True)
            positive_start_checks.append({'n_B':density,'best':candidate,'attempts':attempts})
            if candidate['status']!='STATIONARY_CANDIDATE':
                violations.append({'n_B':density,'violations':candidate['violations'],'scope':'positive_gas_probe'})
        payload={'drip':drip,'transition':transition,'coarse':coarse,'fine':fine,
                 'independent_positive_gas_starts':positive_start_checks}
        with gzip.open(out/'all_attempts.json.gz','wt') as stream:json.dump(payload,stream,allow_nan=False)
        if violations:
            write_json(out/'FAILURES.json',violations);raise RuntimeError(f'{model}: {len(violations)} unresolved rows')
        data=to_eos_data(fine);write_compose(out/'compose',data);restored=read_compose(out/'compose')
        for field in ('nb','Q1','Q2','Q3','Q4','Q5','Q6','Q7','Yn','Yp','Ye','Ymu','Yh','Aav','Zav'):
            if not np.allclose(getattr(data,field),getattr(restored,field),rtol=3e-14,atol=3e-15):raise RuntimeError('native readback mismatch '+field)
        if not np.array_equal(data.Iphase,restored.Iphase):raise RuntimeError('phase changed')
        for field in ('Yn_present','Yp_present','Ye_present','Ymu_present','Yh_present','Abar_present','Zbar_present','C3_applicable','C4_applicable'):
            if not getattr(restored,field).all():raise RuntimeError('composition mask lost '+field)
        checks=check_consistency(restored,verbose=False)
        for key in ('C1','C2','C3','C4','C5','C6','C10','C11'):
            if not checks[key]['pass_']:raise RuntimeError('generated physical check failed '+key)
        source=read_compose(ROOT/'eos_generator/eos_data'/model)
        mask=(source.nb>=fine_n[0])&(source.nb<=fine_n[-1])&(source.Iphase==7)
        compare={}
        for key,values,ref in [('pressure',data.p,source.p),('epsilon',data.eps,source.eps),('Ye',data.Ye,source.Ye)]:
            diff=PchipInterpolator(fine_n,values)(source.nb[mask])-ref[mask]
            compare[key]={'max_absolute':float(max(abs(diff))),'max_relative':float(max(abs(diff)/np.maximum(abs(ref[mask]),1e-20)))}
        domain=physical_domain_report(restored)
        receipt={'model':model,'status':'GENERATED_STATIONARY_NATIVE_COMPOSITION_READBACK_VERIFIED',
            'rows':len(fine),'coarse_rows':len(coarse),'drip':drip,'transition':transition,
            'physical_diagnostics':checks,'sampled_domain_diagnostics':domain,
            'density_interpolation_refinement':{'half_grid_to_coarse':compare_grid(coarse[::2],coarse),'coarse_to_fine':compare_grid(coarse,fine),'row_counts':[len(coarse[::2]),len(coarse),len(fine)]},
            'source_comparison':compare,'source_comparison_rows':int(sum(mask)),
            'source_comparison_scope':'overlap with provider phase7 only; differences do not govern acceptance',
            'density_grid':'3/4 logarithmic through0.9endpoint,1/4linear toclusterendpoint',
            'positive_gas_multistart_probe_count':len(positive_start_checks),
            'max_projected_gradient':max(r['best']['max_projected_gradient'] for r in fine),
            'minimum_local_scaled_hessian_eigenvalue':min(r['best']['minimum_scaled_hessian_eigenvalue'] for r in fine if r['best']['minimum_scaled_hessian_eigenvalue'] is not None),
            'positive_proton_gas_rows':sum(r['best'].get('gas_protons')=='VARIATIONAL_POSITIVE' for r in fine),
            'composition_semantics':{'Yh':'cluster baryon fraction','cluster_number_fraction':'Yh/Aav','Aav':'continuous CLD baryon count','Zav':'continuous CLD proton count'},
            'input_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'data/research_generator_20261116/spg_author_candidates.json',ROOT/'data/research_generator_20261116/candidate_surface/candidate_surface_fit.json']},
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'python_frontend/compose_eos/candidate_gdfm_cld.py',ROOT/'python_frontend/compose_eos/spg_crust_generator.py',ROOT/'python_frontend/compose_eos/gdfm_uniform.py',ROOT/'python_frontend/compose_eos/cold_leptons.py',ROOT/'python_frontend/compose_eos/gm1_uniform.py',ROOT/'python_frontend/compose_eos/physical_domain.py',Path(__file__)]}}
        with gzip.open(out/'receipt.json.gz','wt') as stream:json.dump(receipt,stream,allow_nan=False)
        summary['models'][model]={k:v for k,v in receipt.items() if k in ('status','rows','coarse_rows','source_comparison','max_projected_gradient','minimum_local_scaled_hessian_eigenvalue','positive_proton_gas_rows','density_interpolation_refinement')}
        summary['models'][model]['density_support_fm3']=[float(fine_n[0]),float(fine_n[-1])]
        write_json(args.output/'SUMMARY.json',summary)
        print(json.dumps({'model':model,'rows':len(fine),'status':receipt['status'],'support':summary['models'][model]['density_support_fm3'],'positive_gas_rows':receipt['positive_proton_gas_rows']}),flush=True)
    summary['status']='COMPLETED_SCOPED';summary['runtime_seconds']=time.monotonic()-started
    write_json(args.output/'SUMMARY.json',summary)

if __name__=='__main__':main()
