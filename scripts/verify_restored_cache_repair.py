#!/usr/bin/env python3
"""Replay all archived restored-SPG inner rows after immutable-cache repair."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import gzip,json,hashlib
import numpy as np
from compose_eos.restored_spg_candidate import RestoredSPGCandidate
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/research_finish_20261116/crust_assembly'
def one(model):
    path=BASE/f'restored_spg/{model}/inner_rows.json.gz'
    old=json.loads(gzip.decompress(path.read_bytes()));c=RestoredSPGCandidate(model)
    fields=['epsilon','pressure','mu_B','energy_per_baryon','Y_e','Y_n_free',
            'coordinates','gradient','gas_proton_insertion_margin_MeV','muon_margin_MeV']
    maximum=0.;minimum=float('inf')
    for previous in old:
        current=c.state(previous['n_B'])
        for key in fields:
            difference=np.max(np.abs(np.asarray(current[key])-np.asarray(previous[key])))
            maximum=max(maximum,float(difference))
            if difference!=0:raise RuntimeError(model+' changed '+key+': '+str(difference))
        minimum=min(minimum,current['minimum_scaled_hessian_eigenvalue'])
        if current['artificial_search_bound_active']:raise RuntimeError('active artificial boundary')
    return dict(model=model,rows=len(old),all_physical_fields_bit_identical=True,max_absolute_difference=maximum,
      minimum_scaled_hessian_eigenvalue=minimum,input_rows_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=5) as pool:
        results=list(pool.map(one,[f'SPG_M{i}' for i in range(1,6)]))
    p=ROOT/'python_frontend/compose_eos/restored_spg_candidate.py'
    old=json.loads((BASE/'restored_spg/CAMPAIGN.json').read_text())['source_sha256'][str(p.relative_to(ROOT))]
    report=dict(status='ALL_325_ARCHIVED_PHYSICAL_ROWS_BIT_IDENTICAL_AFTER_IMMUTABLE_CACHE_REPAIR',
      old_code_sha256=old,current_code_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),results=results)
    (BASE/'cache_repair_replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
