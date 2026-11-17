#!/usr/bin/env python3
"""Assemble final scientific execution index from completed, hash-bound evidence."""
from pathlib import Path
import datetime, hashlib, json
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'data/research_finish_20261116'

def read(rel): return json.loads((ROOT/rel).read_text())
def digest(rel): return hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()

def main():
    export_rel='data/research_finish_20261116/geometry/exported/RECEIPT.json'
    exported=read(export_rel)
    if len(exported['models'])!=8: raise RuntimeError('Final eight-model export missing')
    geometry=[]
    for model in exported['models']:
        name=model['model_id']
        path='data/research_finish_20261116/geometry/accelerated_final/'+name+'/SUMMARY.json'
        s=read(path)
        geometry.append({'model_id':name,'density_nodes':s['requested_density_nodes'],
                         'fine_selected_nodes':s['fine_selected_nodes'],'fine_held_nodes':s['fine_held_nodes'],
                         'summary':path,'summary_sha256':digest(path)})
    regression_rel='data/research_finish_20261116/integration/resume_verification/REGRESSION.json'
    regression=read(regression_rel)
    if regression['status']!='PASS' or not regression['tested_sources_unchanged']:raise RuntimeError('Regression not source-stable PASS')
    evidence=[export_rel,regression_rel,
      'data/research_finish_20261116/independent_review/REVIEW_RECEIPT.json',
      'docs/research_finish_20261116/HANDOFF_KO.md',
      'data/research_finish_20261116/geometry/RECEIPT.json',
      'data/research_finish_20261116/geometry/nl3_upper_exported/RECEIPT.json',
      'data/research_finish_20261116/geometry/OBSERVED_PHASES.json',
      'data/research_finish_20261116/integration/ENVIRONMENT.json',
      'data/research_finish_20261116/integration/EVIDENCE_REUSE.json',
      'data/research_finish_20261116/MODEL_COVERAGE.json',
      'data/research_finish_20261116/GOAL_DAG.json',
      'data/research_finish_20261116/RESEARCH_LOOP.json',
      'data/research_finish_20261116/crust_band_recheck/RECEIPT.json',
      'docs/research_finish_20261116/REPORT_KO.md',
      'docs/research_finish_20261116/INDEPENDENT_REVIEW.md',
      'docs/research_finish_20261116/REPRODUCE.md',
      'docs/research_finish_20261116/SOURCE_COMPLETION_ADDENDUM_KO.md']
    result={'schema':'ns-cowling.research-finish-execution.v1',
      'base_commit':'573a31485bc9af30915e01fe74f259e4fa5ac54f',
      'assembled_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'project_status':'PARTIALLY_COMPLETE',
      'scope':'Five cores and thirteen original crusts; independently checked explicit representations/candidates with physical composition. No new core-crust or stellar calculations.',
      'policies':{'compose_matching':'diagnostic_only','physical_support':'density_resolved','source_equivalence':'separate_from_executable_route_availability'},
      'core':{'models':5,'prior_rows_reused_with_source_output_identity':807,'new_full_campaign_this_loop':False},
      'spg':{'models':5,'surface_prescriptions_per_model':5,'grid129_selected_rows':3225,'grid257_selected_rows':6425,'actual_common_density_rows':3225,'first_law_stencils':125,'original_calibration_established':False,'source_transition_prescription_identical':False},
      'partial_crust':{'tables':9,'native_rows':6458,'unchanged_TF_models':4,'unchanged_TF_upper_fm3':.03,'changed_rest_SPG_candidates':5,'changed_rest_SPG_upper_fm3':.01,'whole_original_crusts_completed':0},
      'geometry':{'models':geometry,'target_unique_density_nodes':76,'target_mesh_density_nodes':152,'meshes':[192,384],'dense_nodes_retained':112,'band_nodes_dispatched':40,'band_nodes_completed':40,'upper_extension_new_nodes':4,'unique_mesh_density_nodes_including_extension':156,'native_tables':18,'native_rows_including_two_repeated_anchors':158,'final_held_nodes':0,'global_minimum_proven':False,'continuum_transition_proven':False},
      'regression':{'frontend_pipeline_passed':971,'frontend_pipeline_skipped':9,'generator_passed':142,'existing_generator_warning':'Matplotlib Axes3D import','source_stable':True,'receipt':regression_rel},
      'public_cli':{'catalog_routes':18,'actual_executed_families':['GM1 core','APR phase endpoints','SPG candidate','coarse FSU2H spatial TF'],'whole_original_model_completion_implied':False},
      'evidence_sha256':{p:digest(p) for p in evidence},
      'publication':{'branch':'remediation/2026-11-15-m5a-maxwell-eos-contracts','pr':15,'required_mode':'fast-forward; draft/open/unmerged','authority':'Actual commit SHA, push and CI outcomes are observed from the live branch and PR; this receipt does not self-certify its future publication.'}}
    (BASE/'EXECUTION.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('Wrote final execution index:',sum(x['density_nodes'] for x in geometry),'unique density nodes')

if __name__=='__main__':main()
