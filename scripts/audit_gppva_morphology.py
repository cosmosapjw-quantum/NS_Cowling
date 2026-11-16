#!/usr/bin/env python3
"""Audit retained TF density morphology without discarding lower competitors."""
from __future__ import annotations
import argparse,gzip,json,hashlib
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.tf_morphology import profile_morphology,geometry_label
from compose_eos.tf_selection import select_geometry_candidates
from refine_gppva_geometry import FAMILIES,digest
from build_gppva_tf import dump,physical_geometry

def contrasts(row):
    p=row['profile'];w=np.array(p['weights']);out={}
    for key in ['nn','np','ne']:
        x=np.asarray(p[key]);mean=float(w@x);span=float(np.ptp(x))
        out[key]={'mean_fm3':mean,'absolute_contrast_fm3':span,'relative_contrast':span/max(abs(mean),1e-30)}
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('campaign',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();records=[];inputs={}
    for path in sorted(a.campaign.glob('*/attempts_N*.json.gz')):
        raw=path.read_bytes();inputs[str(path)]=hashlib.sha256(raw).hexdigest()
        for r in json.loads(gzip.decompress(raw)):
            decision=select_geometry_candidates(r['attempts'],FAMILIES.values(),n_b=r['n_b'])
            selected=decision['selected']
            for index,(attempt,diagnostic) in enumerate(zip(r['attempts'],decision['attempt_diagnostics'])):
                if not diagnostic['stationary_supported']:continue
                records.append({'model_id':path.parent.name,'n_b':r['n_b'],'points':r['points'],
                    'attempt_index':index,'requested_family':attempt['requested_family'],
                    'retained_geometry_label':attempt.get('geometry'),'current_geometry_label':geometry_label(attempt),
                    'selected':attempt is selected,'energy_per_baryon':attempt['energy_per_baryon'],
                    'morphology':profile_morphology(attempt),'species_contrasts':contrasts(attempt)})
    selected=[r for r in records if r['selected']]
    report={'input_sha256':inputs,'code_sha256':{n:digest(ROOT/n) for n in [
        'scripts/audit_gppva_morphology.py','python_frontend/compose_eos/tf_morphology.py','scripts/build_gppva_tf.py']},
        'stationary_attempt_count':len(records),'selected_count':len(selected),
        'selected_noncanonical':[r for r in selected if 'NONCANONICAL' in r['morphology']['classification'] or r['morphology']['classification']=='MULTIPLE_RADIAL_DENSE_INTERVALS'],
        'selected_legacy_uniform_with_species_contrast_above_1e_minus4':[r for r in selected if r['retained_geometry_label']=='uniform' and max(v['relative_contrast'] for v in r['species_contrasts'].values())>1e-4],
        'scope':'Descriptive shape and composition audit; all original energies remain physical competitors.',
        'records':records}
    dump(a.output,report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('records','input_sha256','code_sha256')}))
if __name__=='__main__':main()
