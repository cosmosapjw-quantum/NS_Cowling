#!/usr/bin/env python3
"""Stage C: 80 previously uncovered source-bound native matching cells.

No physics or immutable SUPPORT changes. Source tables remain the historical
129-row .08--.5 cores; the B .02--1 audit is cited separately, not transplanted.
The 390-cell merge counts model/method IDs once while preserving source variants.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'python_frontend')]
import generated_pair_campaign as historical
import match_core_crust as pipeline
from compose_eos.native import executable, run_native

CORES = ('GM_GM1', 'RG_SLY4', 'RG_KDE0v1', 'BBB_BHF2')
SMOOTH = ('cubic_F', 'cubic_P', 'quintic_F', 'quintic_P')
OLD = 'data/matched_eos_reference_20261117/campaign/RECEIPT.json'
MAXWELL = 'data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json'
PREP = 'data/maxwell_three_families_20261117/prepared_v2/MANIFEST.json'
UNIVERSE = 'data/current_generator_intake_20261117/MANIFEST_APPROVED_SCOPED.json'
B_DOMAIN = 'data/core_support_20261117/domain_audit/REPORT.json'


def selected_cells():
    cells = []
    for crust in historical.CRUSTS:
        for core in CORES:
            for method in (*SMOOTH, 'css'):
                if core in ('GM_GM1', 'RG_SLY4') and method in ('quintic_P', 'css'):
                    continue
                cells.append(dict(id=f'{crust}__{core}__{method}', crust=crust,core=core,method=method,n_bridge=40,
                                  crust_source_id=f'crust:{crust}',core_source_id=f'core:{core}',
                                  requested=True,eligible=False,attempted=False,admitted=False,reused=False,status='PENDING'))
    return cells


def counts(cells):
    if len({c['id'] for c in cells}) != len(cells):
        raise ValueError('duplicate model/method cell ID')
    return dict(requested=len(cells),eligible=sum(c['eligible'] for c in cells),
                kernel_attempts=sum(c['attempted'] for c in cells),
                constructed_verified=sum(c['status']=='CONSTRUCTED_VERIFIED' for c in cells),
                serialized_readback_verified=sum(c['status'] in ('CONSTRUCTED_VERIFIED','ADMITTED') for c in cells),
                method_admitted=sum(c['admitted'] for c in cells),
                diagnostic_rejected=sum(c['attempted'] and c['status']=='REJECTED' for c in cells),
                unsupported_or_source_hold=sum(not c['eligible'] for c in cells),
                by_status=dict(Counter(c['status'] for c in cells)),
                by_method={m:dict(Counter(c['status'] for c in cells if c.get('method')==m)) for m in (*SMOOTH,'css')})


def source_snapshot(root, family, model):
    # Extend only the explicit runner cohort in this process; the historical
    # function still checks outer/native/request/output/SUPPORT seals in full.
    previous = historical.CORES
    try:
        historical.CORES = CORES
        return historical.source_snapshot(root, family, model)
    finally:
        historical.CORES = previous


def revalidate(root, source):
    if source_snapshot(root,source['family'],source['model']) != source:
        raise ValueError('source identity changed after snapshot')


def verify_output_seal(directory, hashes):
    if not hashes:
        raise ValueError('missing output seal')
    for name, digest in hashes.items():
        if historical.sha(directory/name) != digest:
            raise ValueError(f'reused output changed: {directory/name}')


def verify_maxwell_output(root, record):
    if record.get('admission',{}).get('status') != 'ADMITTED':
        raise ValueError('historical Maxwell independent admission missing/failed')
    if (record.get('readback',{}).get('numerical',{}).get('status') != 'PASS'
            or record.get('readback',{}).get('metadata',{}).get('status') != 'PASS'
            or record.get('serialized_native_verification',{}).get('overall') != 'PASS'):
        raise ValueError('historical Maxwell serialized verification missing/failed')
    verify_output_seal(root/record['output_path'],record['output_sha256'])


def merge_390(root, output, receipt):
    old = historical._json(root / OLD)
    maxwell = historical._json(root / MAXWELL)
    prep = historical._json(root / PREP)
    if historical.sha(root/PREP) != maxwell['prepared_manifest_sha256']:
        raise ValueError('historical Maxwell preparation seal mismatch')
    # Verify native source and prepared payload bytes without rerunning physics.
    for source in prep['source_records']:
        if source['family'] != 'native_generated':
            continue
        for field, location in [('source_table_sha256','source_path'),('internal_sha256','internal_path'),('compose_sha256','compose_path')]:
            for name, digest in source.get(field,{}).items():
                path = root/source[location]/name
                if historical.sha(path) != digest:
                    raise ValueError(f'historical Maxwell source/payload changed: {path}')
        if historical.sha(root/source['source_receipt_path']) != source['source_receipt_sha256']:
            raise ValueError('historical Maxwell receipt seal mismatch')
    universe = historical._json(root/UNIVERSE)['cells']
    mapped = {c['id']:dict(id=c['id'],crust_source_id=c['crust_source_id'],core_source_id=c['core_source_id'],
              method=c['method'],requested=True,eligible=False,kernel_attempt=False,constructed_readback=False,
              method_admitted=False,reused_evidence=False,status='UNSUPPORTED_OR_UNCOVERED',
              reason=c['reason_code'],evidence_variants=[]) for c in universe}
    for c in old['cells']:
        if c['status'] in ('CONSTRUCTED_VERIFIED','ADMITTED','REUSED_ADMITTED'):
            verify_output_seal(root/c['compose_path'],c.get('output_sha256') or c['result']['output_sha256'])
        target=mapped[c['id']]
        variant=dict(path=OLD,sha256=historical.sha(root/OLD),status=c['status'],scope=old['scope'],
                     source_ids=[c['crust_source_id'],c['core_source_id']],actual_matching_core_density_fm3=[.08,.5],
                     original_new_kernel_attempt=c['attempted'],original_reused=c['reused'])
        target['evidence_variants'].append(variant)
        target.update(eligible=c['eligible'],kernel_attempt=c['attempted'] or c['reused'],
                      constructed_readback=c['status'] in ('CONSTRUCTED_VERIFIED','ADMITTED','REUSED_ADMITTED'),
                      method_admitted=c['admitted'],reused_evidence=True,status=c['status'],reason='HISTORICAL_SOURCE_BOUND_EVIDENCE')
    for c in maxwell['cells']:
        if c['family']!='native_generated':
            continue
        target=mapped[c['id']+'__maxwell']
        attempted=any(x.get('kernel_invoked',False) for x in c['candidate_results'])
        generated=[x for x in c['candidate_results'] if x['status']=='GENERATED_VERIFIED_FINITE_DOMAIN']
        if len(generated) != c['generated_count']:
            raise ValueError('historical Maxwell output count mismatch')
        for record in generated:
            verify_maxwell_output(root,record)
        constructed=bool(generated)
        variant=dict(path=MAXWELL,sha256=historical.sha(root/MAXWELL),family=c['family'],status=c['status'],
                     kernel_invocations=c['kernel_invocations'],component_candidates=len(c['candidate_results']),
                     source_ids=[c['crust_source_id'],c['core_source_id']],
                     prepared_manifest=PREP,scope='finite component Maxwell; not full-source stellar EOS')
        target['evidence_variants'].append(variant)
        target.update(eligible=bool(c['candidate_results']),kernel_attempt=target['kernel_attempt'] or attempted,
                      constructed_readback=target['constructed_readback'] or constructed,
                      method_admitted=target['method_admitted'] or constructed,reused_evidence=True,
                      status=c['status'],reason='HISTORICAL_NATIVE_COMPONENT_MAXWELL_EVIDENCE')
    for c in receipt['cells']:
        target=mapped[c['id']]
        target['evidence_variants'].append(dict(path=str((output/'RECEIPT.json').relative_to(root)),status=c['status'],
                   scope=receipt['scope'],actual_matching_core_density_fm3=[.08,.5],new_in_C=True))
        target.update(eligible=c['eligible'],kernel_attempt=c['attempted'],
                      constructed_readback=c['status'] in ('CONSTRUCTED_VERIFIED','ADMITTED'),
                      method_admitted=c['admitted'],reused_evidence=False,status=c['status'],reason='NEW_C_COHORT')
    for cell in mapped.values():
        if cell['evidence_variants']:
            continue
        if cell['crust_source_id'].startswith('crust:SPG_') and cell['core_source_id']=='core:APR_APR':
            cell.update(status='DEFERRED_APR_GIBBS_SUPPORT',reason='retained6_to8_to7 requires boundary-aware construction/consumer support in D')
        else:
            cell.update(status='DEFERRED_TF_PRODUCTION_SUPPORT',reason='finite TF components exist; full-source smooth/CSS eligibility not established; see separate TF lane')
    cells=list(mapped.values())
    assert len(cells)==len({c['id'] for c in cells})==390
    result=dict(schema='ns-cowling.stage-C-matching-union.v1',scope='matching evidence only; TF production eligibility decided separately',
        native_generated_only=True,public_original_or_reconstructed_not_counted=True,
        eligible_scope='eligible for the recorded finite-component/method kernel only; not full-source TF production or full-density certification',
        summary=dict(requested=390,eligible=sum(c['eligible'] for c in cells),
          unique_cell_kernel_attempts=sum(c['kernel_attempt'] for c in cells),
          unique_cell_constructed_readback=sum(c['constructed_readback'] for c in cells),
          unique_cell_method_admitted=sum(c['method_admitted'] for c in cells),
          reused_evidence_cells=sum(c['reused_evidence'] for c in cells),
          deferred_uncovered=sum(not c['evidence_variants'] for c in cells),
          supported_not_yet_attempted_in_selected_C_scope=0,
          permanently_unsupported_claimed=0,
          native_Maxwell_component_kernel_invocations=sum(c['kernel_invocations'] for c in maxwell['cells'] if c['family']=='native_generated'),
          by_status=dict(Counter(c['status'] for c in cells))),cells=cells,
        deduplication='model/method ID counted once; all physical source/domain variants retained under evidence_variants',
        physical_scope='smooth native construction+serialized consistency is not independent method admission; finite Maxwell components do not certify original source intervals',
        inputs_sha256={p:historical.sha(root/p) for p in (OLD,MAXWELL,PREP,UNIVERSE,B_DOMAIN)})
    historical.save(output/'MATCHING_390_MAP.json',result)
    return result


def execute(root,output):
    output.mkdir(parents=True,exist_ok=False)
    cells=selected_cells()
    sources,issues={},{}
    for family,models in [('crust',historical.CRUSTS),('core',CORES)]:
        for model in models:
            sid=f'{family}:{model}'
            try:
                sources[sid]=source_snapshot(root,family,model)
            except Exception as exc:
                issues[sid]=f'{type(exc).__name__}: {exc}'
    deps=sorted((root/'fortran_eos').glob('**/*.f90'))+sorted((root/'python_frontend/compose_eos').glob('*.py'))
    deps += [Path(__file__), root/'scripts/generated_pair_campaign.py', root/'scripts/match_core_crust.py',root/'scripts/match_report_admission.py']
    before={str(p.relative_to(root)):historical.sha(p) for p in deps}
    receipt=dict(schema='ns-cowling.stage-C-supported-pair-expansion.v1',stage='C4_C5',
       base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
       started_utc=datetime.now(timezone.utc).isoformat(),scope='finite native inner-crust+core matching; no whole-star/closure certification',
       source_physics='unchanged native Fortran; explicit thermo-only composition drop',n_bridge=40,workers=1,threads=1,
       skill_context_used=['physmath-research-loop1.3.0 main','POLICY_INTEGRATION','CODING_LOOP C5-C8/C10','superpowers TDD','superpowers verification-before-completion'],
       historical_reuse=dict(core_B_domain_path=B_DOMAIN,sha256=historical.sha(root/B_DOMAIN),
         applies_to=['GM_GM1','RG_SLY4','RG_KDE0v1'],new_core_generation=False,
         output_density_scope='new matches use actual historical .08--.5 source tables; B .02--1 not transplanted',
         BBB_scope='own sealed single-phase source and per-output native consistency/readback; no B domain certification'),
       native_executable=dict(path=str(executable()),sha256=historical.sha(executable())),sources=list(sources.values()),source_issues=issues,
       source_sha256_before=before,cells=cells,status='RUNNING')
    historical.save(output/'MANIFEST.json',receipt)
    historical.save(output/'RECEIPT.json',receipt)
    for cell in cells:
        dest=output/'pairs'/cell['id']
        bad={sid:issues[sid] for sid in (cell['crust_source_id'],cell['core_source_id']) if sid in issues}
        if bad:
            cell.update(status='HOLD_SOURCE',error=bad)
        else:
            cell['eligible']=True
            try:
                cr,co=[sources[cell[f'{side}_source_id']] for side in ('crust','core')]
                revalidate(root,cr); revalidate(root,co)
                if cell['method'] in SMOOTH:
                    result=historical.run_smooth(root/cr['compose_path'],root/co['compose_path'],dest,cell)
                    cell.update(status=result['status'],attempted=result['matching_kernel_invoked'],result=result)
                else:
                    old=pipeline.N_BRIDGE
                    try:
                        pipeline.N_BRIDGE=40
                        result=pipeline.run_pipeline('css',root/cr['compose_path'],root/co['compose_path'],dest)
                    finally:
                        pipeline.N_BRIDGE=old
                    cell.update(status='ADMITTED' if result.get('survived') else 'REJECTED' if result.get('matching_kernel_invoked') else 'NOT_EVALUATED',
                                attempted=bool(result.get('matching_kernel_invoked')),admitted=bool(result.get('survived')),result=result)
                if cell['status'] in ('CONSTRUCTED_VERIFIED','ADMITTED'):
                    cell['compose_path']=str((dest/'compose').relative_to(root))
                cell['method_admission_scope']='independent CSS interface admission' if cell['admitted'] else 'native smooth construction+serialized consistency only' if cell['status']=='CONSTRUCTED_VERIFIED' else 'not admitted'
            except Exception as exc:
                cell.update(status='NOT_EVALUATED',error=f'{type(exc).__name__}: {exc}',failure_category='data_provenance_or_environment')
        historical.save(dest/'PAIR.json',cell)
        receipt['census']=counts(cells)
        historical.save(output/'RECEIPT.json',receipt)
        print(cell['id'],cell['status'],flush=True)
    after={str(p.relative_to(root)):historical.sha(p) for p in deps}
    receipt.update(source_sha256_after=after,sources_unchanged=before==after,
      finished_utc=datetime.now(timezone.utc).isoformat(),status='COMPLETE_SCOPED' if before==after else 'IDENTITY_CHANGED')
    historical.save(output/'RECEIPT.json',receipt)
    merge_390(root,output,receipt)
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--native-exe',type=Path,required=True)
    args=parser.parse_args()
    os.environ.update(COMPOSE_EOS_FORTRAN_EXE=str(args.native_exe.resolve()),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    output=args.out.resolve()
    if not output.is_relative_to(ROOT):
        raise ValueError('output must be in repository')
    result=execute(ROOT,output)
    print(json.dumps(result['census'],indent=2))
    return 0 if result['status']=='COMPLETE_SCOPED' else 1

if __name__=='__main__':
    raise SystemExit(main())
