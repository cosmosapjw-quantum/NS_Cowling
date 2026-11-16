#!/usr/bin/env python3
"""Execute spatial TF geometry/mesh campaigns and export sparse CompOSE data."""
from __future__ import annotations
# Historical research oracle; production physics lives in the Fortran backend.
# Imports remain available to pinned-campaign audit and regression helpers.
import sys as _oracle_sys
if __name__ == '__main__':
    if '--historical-python-oracle' not in _oracle_sys.argv:
        raise SystemExit('This is a historical Python research oracle. Use compose-eos generate for production; add --historical-python-oracle only to reproduce the explicitly archived Python campaign.')
    _oracle_sys.argv.remove('--historical-python-oracle')
import argparse,csv,gzip,hashlib,json,os,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from compose_eos.spatial_tf_equilibrium import solve_screened_cell as solve_cell
from compose_eos import gppva_rmf as rmf
from compose_eos.api import write_compose,read_compose,check_consistency
from compose_eos.eos_data import EOSData
from compose_eos.physical_domain import physical_domain_report
from compose_eos.tf_selection import select_geometry_candidates
from compose_eos.tf_morphology import geometry_label,profile_morphology,TF_PHASE_CODES
from compose_eos.tf_search import weak_cluster_seed

MODELS=('GPPVA_DD2','GPPVA_DDME2','GPPVA_TW','GPPVA_FSU2','GPPVA_FSU2R','GPPVA_FSU2H','GPPVA_NL3wrL55','GPPVA_TM1e')
GEOMETRIES=((3,'droplet'),(2,'droplet'),(1,'droplet'),(2,'hole'),(3,'hole'),(3,'uniform'))
DEFAULT_DENSITIES=(.0002,.0005,.001,.002,.005,.01,.02,.03,.04,.05,.06,.07,.08,.09,.10)

def serial(value):
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,np.ndarray):return value.tolist()
    raise TypeError(type(value).__name__)

def dump(path,value):
    text=json.dumps(value,indent=2,allow_nan=False,default=serial)+'\n'
    path=Path(path)
    if path.name.endswith('.gz'):
        with gzip.GzipFile(filename=str(path),mode='wb',mtime=0) as f:f.write(text.encode())
    else:path.write_text(text)

def physical_geometry(row):
    return geometry_label(row)

def export(model,states,out):
    n=np.array([s['n_b'] for s in states]);nrows=len(n)
    if nrows<3:return {'status':'INSUFFICIENT_STATIONARY_ROWS','row_count':nrows}
    m=rmf.model(model);zero=np.zeros(nrows);true=np.ones(nrows,bool)
    energy=np.array([s['epsilon'] for s in states]);pressure=np.array([s['pressure'] for s in states])
    mub=np.array([s['mu_B'] for s in states]);muq=np.array([s['mu_Q_lepton'] for s in states])
    yn=np.array([s['Y_n'] for s in states]);yp=np.array([s['Y_p'] for s in states]);ye=np.array([s['Y_e'] for s in states])
    geometry=[physical_geometry(s) for s in states];uniform=np.array([k=='uniform' for k in geometry])
    phase=np.array([TF_PHASE_CODES[k] for k in geometry])
    qenergy=energy/(n*m['mn'])-1
    data=EOSData(nb=n,mn=m['mn'],mp=m['mp'],Q1=pressure/n,Q2=zero,Q3=mub/m['mn']-1,
        Q4=-muq/m['mn'],Q5=zero,Q6=qenergy,Q7=qenergy.copy(),q_extra=zero,Iphase=phase,
        Yn=np.where(uniform,yn,0.),Yp=np.where(uniform,yp,0.),Ye=ye,Ymu=zero,
        Aav=np.where(uniform,1.,0.),Zav=zero,Yh=zero,composition_file_present=True,
        composition_available=True,beta_equilibrated=True,row_status=np.where(uniform,2,1),
        Yn_present=uniform,Yp_present=uniform,Ye_present=true,Ymu_present=true,
        Abar_present=uniform,Zbar_present=uniform,Yh_present=uniform,
        C3_applicable=uniform,C4_applicable=uniform)
    write_compose(out/'compose',data);back=read_compose(out/'compose')
    for key in ('nb','Q1','Q3','Q4','Q6','Q7','Ye','Ymu'):
        if not np.allclose(getattr(back,key),getattr(data,key),rtol=2e-13,atol=2e-14):raise ValueError('native readback changed '+key)
    for key in ('Iphase','Yn_present','Yp_present','Ye_present','Ymu_present','Yh_present','C3_applicable','C4_applicable'):
        if not np.array_equal(getattr(back,key),getattr(data,key)):raise ValueError('native readback changed '+key)
    composition={'composition_representation':'spatial_TF_species_densities',
        'total_nucleon_fractions':'Integrals of local neutron and proton density, including diffuse clustered matter.',
        'compose_semantics':'Ye/Ymu present in every row. Free-neutron/free-proton/discrete-cluster masks unavailable in spatial TF rows; full free-nucleon fields only for uniform rows.',
        'profiles':'selected_states.json.gz','phase_codes':TF_PHASE_CODES,
        'row_count':nrows,'rows':[{'n_b':s['n_b'],'Y_n_total':s['Y_n'],'Y_p_total':s['Y_p'],'Y_e':s['Y_e'],'Y_mu':0.,'geometry':g,'morphology':profile_morphology(s),'radius_fm':s['radius_fm'],
            'baryon_constraint_error':s['baryon_constraint_error'],'charge_constraint_error':s['charge_constraint_error']} for s,g in zip(states,geometry)]}
    dump(out/'composition.json',composition)
    dump(out/'selected_states.json.gz',states)
    files=[*sorted((out/'compose').glob('eos.*')),out/'composition.json',out/'selected_states.json.gz']
    if (out/'INTERVAL_SUPPORT.json').is_file():files.append(out/'INTERVAL_SUPPORT.json')
    output_hashes={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return {'output_sha256':output_hashes,'status':'EXPORTED_NATIVE_READBACK_VERIFIED','row_count':nrows,'native_checks':check_consistency(back,verbose=False),
            'physical_domain':physical_domain_report(back),'composition_representation':'spatial_TF_species_densities',
            'free_cluster_decomposition':'NOT_AVAILABLE in nonuniform TF; physically computed total composition is in companion composition.json and profiles'}

def export_consecutive(model,states,densities,out):
    """Export consecutive supported runs; retain unresolved gaps explicitly.

    Phase codes describe finite sampled candidates. Even within a complete run,
    a phase change does not authorize continuous interpolation or certify a
    common-tangent transition.
    """
    lookup={float(s['n_b']):s for s in states}
    blocks=[];block=[];gaps=[]
    for density in densities:
        state=lookup.get(float(density))
        if state is None:
            if block:blocks.append(block);block=[]
            gaps.append(float(density))
        else:block.append(state)
    if block:blocks.append(block)
    phase_changes=[]
    for block in blocks:
        for left,right in zip(block,block[1:]):
            a,b=physical_geometry(left),physical_geometry(right)
            if a!=b:
                phase_changes.append({'density_bracket_fm3':[left['n_b'],right['n_b']],
                                      'phases':[a,b]})
    policy={'unresolved_requested_densities_fm3':gaps,
            'phase_change_brackets':phase_changes,
            'interpolation_policy':'No interpolation across unresolved density gaps or phase-change brackets. Same-phase interpolation still needs a density-refinement error estimate.'}
    dump(out/'INTERVAL_SUPPORT.json',policy)
    if not gaps:
        return {**export(model,states,out),'interval_support':policy}
    exports=[]
    for i,block in enumerate(blocks):
        dest=out/f'segment_{i:03d}'
        dest.mkdir()
        dump(dest/'selected_states.json.gz',block)
        exports.append({'path':str(dest.relative_to(out)),
                        'density_interval_fm3':[block[0]['n_b'],block[-1]['n_b']],
                        **export(model,block,dest)})
    return {'status':'PARTIAL_SEGMENTS_NATIVE_READBACK_VERIFIED' if
            any(e['status']=='EXPORTED_NATIVE_READBACK_VERIFIED' for e in exports)
            else 'INSUFFICIENT_STATIONARY_ROWS',
            'row_count':len(states),'segments':exports,'interval_support':policy}

def recover_unresolved_families(model,n,points,attempts):
    """Finite weak-profile retries; original failures remain in attempts.

    Every retry re-solves the unchanged physical functional and frees the
    radius. A collapse to uniform matter resolves that seed-family attempt,
    not the existence/absence of every possible metastable branch.
    """
    uniform=next((r for r in attempts if r.get('seed')=='uniform' and
        'profile' in r and select_geometry_candidates([r],[(r['dimension'],r['seed'])],
                                                     n_b=n)['selected'] is not None),None)
    if uniform is None:return
    for dimension,seed in GEOMETRIES:
        if any(select_geometry_candidates([r],[(dimension,seed)],n_b=n)['selected'] is not None
               for r in attempts if (r.get('dimension'),r.get('seed'))==(dimension,seed)):
            continue
        resolved=False
        for amplitude in (.03,.3):
            for radius in (15.,30.,60.):
                signed=-amplitude if seed=='hole' else amplitude
                metadata={'n_b':n,'points':points,'dimension':dimension,'seed':seed,
                          'initial_kind':'WEAK_CLUSTER_RECOVERY','initial_amplitude':signed,
                          'initial_radius_fm':radius}
                try:
                    guess=weak_cluster_seed(uniform,dimension,signed,radius)
                    initial={k:np.asarray(guess['profile'][k],float) for k in ('nn','np','ne')}
                    initial['radius']=radius
                    result=solve_cell(model,n,points=points,dimension=dimension,seed=seed,
                        radius=radius,radius_bounds=(5.,180.),maxiter=900,ftol=1e-13,initial=initial)
                    result.update(initial_kind=metadata['initial_kind'],
                                  initial_amplitude=signed,initial_radius_fm=radius)
                    result['geometry']=physical_geometry(result);attempts.append(result)
                    resolved=select_geometry_candidates([result],[(dimension,seed)],n_b=n)['selected'] is not None
                except Exception as exc:
                    attempts.append(dict(metadata,execution='FAILED',exception=type(exc).__name__,message=str(exc)))
                if resolved:break
            if resolved:break

def run_model(model,densities,resolutions,outroot,warm_start=None):
    start=time.monotonic();out=Path(outroot)/model;out.mkdir()
    stages=[];previous={};selected_by_mesh=[]
    if warm_start:
        oldpath=Path(warm_start)/model/f'attempts_N{resolutions[0]:03}.json.gz'
        if oldpath.exists():
            for row in json.loads(gzip.decompress(oldpath.read_bytes())):
                for a in row['attempts']:
                    if 'profile' in a:previous[(float(row['n_b']),a['dimension'],a['seed'])]=a
    for points in resolutions:
        records=[];selected=[]
        for n in densities:
            attempts=[]
            for dimension,seed in GEOMETRIES:
                key=(float(n),dimension,seed)
                initial=None
                if key in previous:
                    old=previous[key];p=old['profile'];x=(np.arange(points)+.5)/points
                    vals={k:np.interp(x,p['x'],p[k]) for k in ('nn','np','ne')}
                    # Conservative correction after profile interpolation.
                    edges=np.linspace(0,1,points+1);w=np.diff(edges**dimension)
                    fac=n/(w@(vals['nn']+vals['np']));vals['nn']*=fac;vals['np']*=fac
                    vals['ne']*=float(w@vals['np'])/max(float(w@vals['ne']),1e-30)
                    initial={**vals,'radius':old['radius_fm']}
                try:
                    r=solve_cell(model,float(n),points=points,dimension=dimension,seed=seed,
                        radius=45. if n<.003 else 30.,radius_bounds=(5.,180.),
                        maxiter=900,ftol=1e-13,initial=initial)
                    r['geometry']=physical_geometry(r);attempts.append(r);previous[key]=r
                except Exception as exc:
                    attempts.append({'n_b':float(n),'points':points,'dimension':dimension,'seed':seed,
                        'execution':'FAILED','exception':type(exc).__name__,'message':str(exc)})
            recover_unresolved_families(model,float(n),points,attempts)
            decision=select_geometry_candidates(attempts,GEOMETRIES,n_b=float(n))
            best=decision['best'];selection=decision['status']
            if decision['selected'] is not None:selected.append(decision['selected'])
            records.append({'n_b':float(n),'selection':selection,'selected_geometry':best['geometry'] if best else None,
                'selection_diagnostics':{k:v for k,v in decision.items() if k not in ('best','selected')},'attempts':attempts})
            dump(out/f'progress_N{points:03}.json',{'status':'RUNNING','completed':len(records),'total':len(densities),'last_density':float(n)})
        dump(out/f'attempts_N{points:03}.json.gz',records)
        summary={'points':points,'requested_points':len(densities),'selected_stationary_points':len(selected),
            'selections':[{'n_b':r['n_b'],'selection':r['selection'],'geometry':r['selected_geometry']} for r in records]}
        dump(out/f'progress_N{points:03}.json',{'status':'COMPLETE',**summary})
        stages.append(summary);selected_by_mesh.append(selected)
    final=selected_by_mesh[-1];dump(out/'selected_states.json.gz',final)
    mesh_comparisons=[]
    for coarse,fine in zip(selected_by_mesh,selected_by_mesh[1:]):
        lookup={r['n_b']:r for r in coarse}
        for f in fine:
            if f['n_b'] not in lookup:continue
            c=lookup[f['n_b']]
            mesh_comparisons.append({'n_b':f['n_b'],'coarse':c['points'],'fine':f['points'],
                'delta_energy_per_baryon_MeV':f['energy_per_baryon']-c['energy_per_baryon'],
                'relative_pressure_change':(f['pressure']-c['pressure'])/max(abs(f['pressure']),1e-30),
                'delta_Yp':f['Y_p']-c['Y_p'],'coarse_geometry':c['geometry'],'fine_geometry':f['geometry']})
    export_result=export_consecutive(model,final,densities,out)
    report={'model_id':model,'status':'EXECUTED','seconds':time.monotonic()-start,'stages':stages,
        'mesh_comparison':mesh_comparisons,'export':export_result,
        'global_minimum_status':'Finite geometry/seed/radius search only; unresolved families and feasible lower unconverged attempts retained.',
        'source_matching':'Diagnostic only; no source-table agreement threshold used for physical acceptance.',
        'electron_model':'local zero-temperature screened npe; muon threshold checked',
        'radius_search_fm':[5.,180.]}
    dump(out/'RECEIPT.json',report)
    return {'model':model,'rows':len(final),'requested_rows':len(densities),'seconds':report['seconds'],'export_status':export_result['status']}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);parser.add_argument('--models',nargs='+',default=list(MODELS))
    parser.add_argument('--densities',nargs='+',type=float,default=DEFAULT_DENSITIES)
    parser.add_argument('--resolutions',nargs='+',type=int,default=[48,96]);parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--warm-start',type=Path,help='Previous campaign profiles used only as optimization initial guesses')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    sources=['python_frontend/compose_eos/spatial_tf.py','python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/gppva_rmf.py','python_frontend/compose_eos/cold_leptons.py','scripts/build_gppva_tf.py','python_frontend/compose_eos/tf_selection.py','python_frontend/compose_eos/tf_morphology.py','python_frontend/compose_eos/tf_search.py']
    initial_sha={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources}
    dump(args.output/'INPUT.json',{'models':args.models,'densities_fm3':args.densities,'resolutions':args.resolutions,'warm_start':str(args.warm_start) if args.warm_start else None,'source_sha256':initial_sha})
    if any(m not in MODELS for m in args.models):raise ValueError('unknown model')
    if any(np.diff(args.densities)<=0) or any(np.diff(args.resolutions)<=0):raise ValueError('increasing density/resolution required')
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(run_model,m,args.densities,args.resolutions,str(args.output),args.warm_start):m for m in args.models}
        for job in as_completed(jobs):
            try:r=job.result()
            except Exception as e:r={'model':jobs[job],'status':'FAILED','exception':type(e).__name__,'message':str(e)}
            results.append(r);print(json.dumps(r),flush=True)
    final_sha={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources}
    if final_sha!=initial_sha:raise RuntimeError('campaign source files changed during execution')
    dump(args.output/'CAMPAIGN.json',{'models':results,'densities_fm3':args.densities,'resolutions':args.resolutions,
        'source_sha256':initial_sha,'source_unchanged_during_execution':True})
    return 1 if any(r.get('rows',0)!=len(args.densities) or r.get('status')=='FAILED' for r in results) else 0
if __name__=='__main__':raise SystemExit(main())
