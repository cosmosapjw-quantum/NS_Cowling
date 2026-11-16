#!/usr/bin/env python3
"""Reclassify retained geometry solutions and export with separate provenance."""
from __future__ import annotations
import argparse,gzip,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.tf_morphology import geometry_label,profile_morphology
from compose_eos.tf_selection import select_geometry_candidates
from refine_gppva_geometry import FAMILIES,digest,transition_evidence
from build_gppva_tf import dump,export_consecutive

def classified(row):
    return dict(row,geometry_at_solve=row.get('geometry'),geometry=geometry_label(row),
                morphology=profile_morphology(row))

def additional_attempts(directories):
    lookup={};hashes={}
    for directory in directories:
        for path in sorted(Path(directory).glob('*.json.gz')):
            result=json.loads(gzip.decompress(path.read_bytes()))
            if not isinstance(result,dict) or 'attempts' not in result:continue
            sha=digest(path);hashes[str(path)]=sha
            key=(result['model_id'],round(result['n_b'],12),result['points'])
            for index,row in enumerate(result['attempts']):
                lookup.setdefault(key,[]).append(dict(row,supplemental_attempt_source=str(path),supplemental_attempt_index=index,supplemental_attempt_source_sha256=sha))
    return lookup,hashes

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('campaign',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--additional',nargs='*',type=Path,default=[])
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);supplement,inputs=additional_attempts(a.additional);reports=[]
    sources={n:digest(ROOT/n) for n in ['scripts/export_gppva_geometry.py','scripts/build_gppva_tf.py',
        'scripts/refine_gppva_geometry.py','python_frontend/compose_eos/tf_morphology.py','python_frontend/compose_eos/tf_selection.py']}
    for path in sorted(a.campaign.glob('*/SUMMARY.json')):
        model=path.parent.name;original=json.loads(path.read_text());inputs[str(path)]=digest(path)
        records=[];exports=[]
        for points in original['meshes']:
            rawpath=path.parent/f'attempts_N{points:03}.json.gz';inputs[str(rawpath)]=digest(rawpath)
            stage=json.loads(gzip.decompress(rawpath.read_bytes()))
            for record in stage:
                record['_base_attempt_count']=len(record['attempts'])
                record['attempts'].extend(supplement.get((model,round(record['n_b'],12),points),[]))
            records.extend(stage)
            rows=[];evidence=[]
            for record_index,r in enumerate(stage):
                decision=select_geometry_candidates(r['attempts'],FAMILIES.values(),n_b=r['n_b'])
                r['status_at_original_run']=r['status'];r['status']=decision['status']
                row=decision['selected']
                evidence.append({'n_b':r['n_b'],'points':points,'base_attempt_file':str(rawpath),'base_record_index':record_index,
                    'supplemental_attempts':[{'file':x['supplemental_attempt_source'],'attempt_index':x['supplemental_attempt_index']} for x in r['attempts'][r['_base_attempt_count']:] if 'supplemental_attempt_source' in x],
                    'selection':{k:v for k,v in decision.items() if k not in ('best','selected')}})
                if row:rows.append(classified(row))
            out=a.output/model/f'N{points:03}';out.mkdir(parents=True)
            dump(out/'selected_states.json.gz',rows)
            dump(out/'EVIDENCE.json',{'rows':evidence,'order':'Base record attempts first, then supplemental_attempts in the listed order.'})
            result=export_consecutive(model,rows,[r['n_b'] for r in stage],out)
            receipt={'model_id':model,'points':points,'export':result,'selection_evidence':'EVIDENCE.json','input_sha256':{str(rawpath):inputs[str(rawpath)]},
                'supplemental_attempt_sources':{str(d):digest(d/'INPUT.json') for d in a.additional},
                'selection_status':[{'n_b':r['n_b'],'original':r['status_at_original_run'],'combined':r['status']} for r in stage],
                'export_source_sha256':sources,'physical_solve_provenance':str(a.campaign/'INPUT.json'),
                'scope':'Unchanged stationary solutions with corrected descriptive geometry and sparse composition masks; no physical energy competitor removed.'}
            dump(out/'RECEIPT.json',receipt);exports.append(receipt)
        transitions=transition_evidence(original['windows'],records,original['meshes'])
        for window in transitions:
            for sample in window['samples']:
                fine=sample['mesh_values'].get(str(original['meshes'][-1]))
                shift=sample['adjacent_mesh_deltaE_change_MeV_per_baryon']
                deltaP=fine['pressures'][1]-fine['pressures'][0] if fine else None
                sample['linearized_mesh_density_shift_fm3']=(abs(sample['n_b']**2*shift/deltaP) if deltaP and shift is not None else None)
                sample['linearization_scope']='Local envelope-slope estimate only; neither an error bound nor an evaluated transition root.'
        report={'model_id':model,'exports':exports,
            'transitions':transitions,
            'scope':'Named simple-morphology comparisons omit noncanonical labels only from those named pairs; all candidate energies remain in raw records.'}
        dump(a.output/model/'SUMMARY.json',report);reports.append(report)
    dump(a.output/'RECEIPT.json',{'input_sha256':inputs,'export_source_sha256':sources,'models':reports,
        'raw_campaign':str(a.campaign),'raw_campaign_unchanged':all(digest(Path(n))==h for n,h in inputs.items())})
    print(json.dumps({'models':len(reports),'exports':sum(len(r['exports']) for r in reports)}))
if __name__=='__main__':main()
