#!/usr/bin/env python3
"""Derive finite observed phase changes from actual corrected stationary profiles."""
from __future__ import annotations
import argparse,gzip,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from refine_gppva_geometry import FAMILIES,DEPENDENCIES,branch_minima,transition_evidence,digest
from compose_eos.tf_selection import select_geometry_candidates
from compose_eos.tf_morphology import geometry_label
from build_gppva_tf import dump

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('primary',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--extension',nargs='*',type=Path,default=[])
    a=p.parse_args();inputs={};by_model={};origins={}
    for base in [a.primary,*a.extension]:
        for path in sorted(base.glob('*/attempts_N*.json.gz')):
            inputs[str(path)]=digest(path)
            rows=json.loads(gzip.decompress(path.read_bytes()))
            model=path.parent.name;lookup=by_model.setdefault(model,{})
            for index,r in enumerate(rows):
                key=(round(r['n_b'],12),r['points'])
                origin={'file':str(path),'record_index':index}
                if key in lookup:
                    if lookup[key]['attempts']!=r['attempts']:
                        raise ValueError('conflicting repeated anchor attempts')
                    origins[(model,*key)].append(origin)
                else:
                    lookup[key]=r;origins[(model,*key)]=[origin]
    reports=[]
    for model,lookup in sorted(by_model.items()):
        meshes=sorted({m for _,m in lookup});records=list(lookup.values())
        windows=[];sequences={}
        for mesh in meshes:
            stage=sorted([r for r in records if r['points']==mesh],key=lambda r:r['n_b'])
            sequence=[]
            for r in stage:
                decision=select_geometry_candidates(r['attempts'],FAMILIES.values(),n_b=r['n_b'])
                q=decision['selected']
                item={'n_b':r['n_b'],'status':decision['status'],'geometry':geometry_label(q) if q else None,
                    'evidence':origins[(model,round(r['n_b'],12),mesh)],
                    'best_attempt_index':decision.get('best_attempt_index')}
                if q:item.update({k:q.get(k) for k in ['energy_per_baryon','pressure','mu_B','Y_n','Y_p','Y_e','Y_mu',
                    'max_chemical_KKT_residual_MeV','radius_fm']})
                sequence.append(item)
            sequences[str(mesh)]=sequence
            for left,right in zip(sequence,sequence[1:]):
                if left['geometry'] and right['geometry'] and left['geometry']!=right['geometry']:
                    window={'left':left['n_b'],'right':right['n_b'],
                        'geometries':[left['geometry'],right['geometry']]}
                    if window not in windows:windows.append(window)
        windows.sort(key=lambda w:(w['left'],w['right'],w['geometries']))
        comparisons=transition_evidence(windows,records,meshes)
        for w in comparisons:
            w['definition']='Change in selected descriptive profile among adjacent evaluated density nodes.'
            w['optimizer_uncertainty']='Stationarity criteria are checked, but no bound on residual energy error is established; mesh and roundoff flags do not certify an energy-order sign.'
            for sample in w['samples']:
                for mesh,value in sample['mesh_values'].items():
                    record=lookup[(round(sample['n_b'],12),int(mesh))]
                    branches=branch_minima(record['attempts'])
                    value['stationarity_diagnostics']=[{k:branches[g].get(k) for k in
                        ['max_chemical_KKT_residual_MeV','radius_stationarity_MeV',
                         'scalar_pde_residual','vector_pde_residual','electron_pde_residual',
                         'baryon_constraint_error','charge_constraint_error','radius_fm','dimension']}
                        for g in w['geometries']]
            w['simple_physical_morphologies']=all(g in ('droplet','rod','slab','tube','bubble','uniform')
                for g in w['geometries'])
        reports.append({'model_id':model,'sequences':sequences,'observed_changes':comparisons})
    sources=list(dict.fromkeys(DEPENDENCIES+[str(Path(__file__).resolve().relative_to(ROOT))]))
    dump(a.output,{'input_sha256':inputs,'source_sha256':{n:digest(ROOT/n) for n in sources},
        'unique_mesh_density_nodes':sum(len(x) for x in by_model.values()),'models':reports,
        'scope':'Finite observed profile ordering only. Historical requested pairs are retained separately. Near-uniform descriptive changes and branch disappearance are not invented sign roots; no common-tangent coexistence or continuum/global-minimum certificate.'})
    print(json.dumps({'models':len(reports),'nodes':sum(len(x) for x in by_model.values()),
        'observed_changes':sum(len(r['observed_changes']) for r in reports)}))
if __name__=='__main__':main()
