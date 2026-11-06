"""Fixed TOV background transport with its executed Fortran provenance.

This is not a mode calibration, a cryptographic signature, or an EOS inference tool.
"""
from __future__ import annotations
import csv,json,math,shutil
from pathlib import Path
from tools.atomic_io import read_json,write_json,digest
from tools.run_case import read_completed_run,RunRecordMismatch,RunNotComplete
FIELDS={'q':'q','r_m':'r','m_kg':'m','p_pa':'p','energy_j_m3':'eps',
        'n_b_m3':'n','mu_j':'mu','chi':'chi','nu':'nu','lambda':'lambda','dnu_dr':'dnu'}
POLICY='finite_cutoff_schwarzschild_time'

def profile_rows(result):
    n=result.get('profile_n')
    if type(n) is not int or not 2<=n<=99999:raise ValueError('Missing/bad profile count')
    try:rows=[{key:result[f'{tag}_{i:05d}'] for key,tag in FIELDS.items()} for i in range(1,n+1)]
    except KeyError as e:raise ValueError('Full profile was not emitted') from e
    validate_rows(rows)
    return rows

def validate_rows(rows):
    if len(rows)<2:raise ValueError('Insufficient background')
    for i,r in enumerate(rows):
        if set(r)!=set(FIELDS) or not all(type(v) in (int,float) and math.isfinite(v) for v in r.values()):
            raise ValueError('Nonfinite or malformed background row')
        if min(r[k] for k in ['r_m','m_kg','p_pa','energy_j_m3','n_b_m3','mu_j','chi','dnu_dr'])<=0 or r['chi']>1:
            raise ValueError('Invalid smooth barotropic background')
        if i and not (r['r_m']>rows[i-1]['r_m'] and r['m_kg']>=rows[i-1]['m_kg'] and r['p_pa']<rows[i-1]['p_pa']):
            raise ValueError('Nonmonotone background segment')

def write_background(result,destination):
    d=Path(destination)
    if d.exists() or d.is_symlink():raise FileExistsError(d)
    if result.get('status')!='ok' or not result.get('result_valid') or result.get('eos_kind')!='polytrope':
        raise ValueError('Only successful physical polytrope backgrounds are exportable here')
    if result.get('surface_policy')!=POLICY:raise ValueError('Unexpected surface')
    try:
        src=Path(result['_evidence']['run_receipt']).parent
        bound=read_completed_run(src,result['request_token'],result['source_sha256'])
    except (KeyError,RunRecordMismatch,RunNotComplete) as e:raise ValueError('Unverified source run') from e
    if any(result.get(k)!=v for k,v in bound.items()):raise ValueError('Result was modified after the run')
    rows=profile_rows(bound);s=read_json(src/'run.json')
    d.mkdir(parents=True,exist_ok=False)
    with (d/'background.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(FIELDS));w.writeheader();w.writerows(rows)
    for name in ['run.json','input.nml','result.json','process/stdout.txt','process/stderr.txt','process/command.json']:
        p=d/'provenance'/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src/name,p)
    meta={'schema_version':1,'status':'complete','kind':'fixed_background','units':'SI',
          'closure':'barotropic_equilibrium_derivative','surface_policy':POLICY,
          'coordinate':'log_pressure','eos_kind':'polytrope','exact_vacuum_surface':False,
          'validated_scope':'TOV_ONLY','mode_calibration':False,'row_count':len(rows),
          'source_sha256':s['source_sha256'],'spec_sha256':s['spec_sha256'],
          'request_token':s['request_token'],'input_sha256':s['input_sha256'],
          'constants_source_sha256':s['constants_source_sha256'],
          'csv_sha256':digest(d/'background.csv'),'pc_pa':bound['pc_pa'],'p_cut_pa':bound['p_cut_pa'],
          'r0_m':bound['r0_m'],'radius_m':bound['radius_m'],'mass_kg':bound['mass_kg'],
          'limits':'finite positive-pressure cutoff and time anchor; no modes or exact vacuum surface'}
    write_json(d/'metadata.json',meta)
    return meta

def read_background(directory,expected_request_token=None):
    d=Path(directory)
    try:
        if d.is_symlink() or (d/'metadata.json').is_symlink() or (d/'background.csv').is_symlink():raise ValueError('Background symlink')
        m=read_json(d/'metadata.json')
        if m.get('schema_version')!=1 or m.get('status')!='complete' or m.get('units')!='SI':raise ValueError('Background contract')
        if m.get('kind')!='fixed_background' or m.get('coordinate')!='log_pressure' or m.get('eos_kind')!='polytrope':raise ValueError('Wrong background interpretation')
        if m.get('exact_vacuum_surface') is not False:raise ValueError('Unsupported vacuum-surface claim')
        if m.get('closure')!='barotropic_equilibrium_derivative' or m.get('surface_policy')!=POLICY:raise ValueError('Wrong physical contract')
        if m.get('mode_calibration') is not False or m.get('validated_scope')!='TOV_ONLY':raise ValueError('Unsupported verification claim')
        if expected_request_token is not None and m['request_token']!=expected_request_token:raise ValueError('Wrong background identity')
        if digest(d/'background.csv')!=m['csv_sha256']:raise ValueError('Background byte mismatch')
        bound=read_completed_run(d/'provenance',m['request_token'],m['source_sha256'])
        run=read_json(d/'provenance/run.json')
        if m['spec_sha256']!=run['spec_sha256']:raise ValueError('Spec identity mismatch')
        if m['input_sha256']!=run['input_sha256'] or m['constants_source_sha256']!=run['constants_source_sha256']:raise ValueError('Input/source mismatch')
        if bound.get('eos_kind')!='polytrope':raise ValueError('Formal benchmark not a physical background')
        with (d/'background.csv').open(newline='') as f:
            reader=csv.DictReader(f)
            if reader.fieldnames!=list(FIELDS):raise ValueError('Background columns')
            rows=[{k:float(v) for k,v in r.items()} for r in reader]
        if len(rows)!=m['row_count'] or rows!=profile_rows(bound):raise ValueError('CSV does not match executed result')
        for key in ['pc_pa','p_cut_pa','r0_m','radius_m','mass_kg']:
            if m[key]!=bound[key]:raise ValueError('Altered background scalars')
        validate_rows(rows)
        return rows,m
    except (OSError,KeyError,TypeError,RunRecordMismatch,RunNotComplete) as e:
        raise ValueError('Invalid background or missing provenance') from e
