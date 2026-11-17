from pathlib import Path
import json
root=Path(__file__).resolve().parents[3]
d=root/'data/research_finish_20261116'
prior=json.loads((root/'data/research_completion_20261116/MODEL_COVERAGE.json').read_text())
rows=[]
for old in prior['models']:
 model=old['model_id']
 row={'model_id':model,'component':old['component'],
      'prior_evidence_index':'data/research_completion_20261116/MODEL_COVERAGE.json',
      'source_comparison':'Diagnostic; not an unconditional matching gate'}
 if old['component']=='core':
  row.update(generation_status='SCOPED_CORE_GENERATION_AVAILABLE',
    evidence='data/research_completion_20261116/core/generated/'+model+'/RECEIPT.json',
    evidence_reuse='data/research_finish_20261116/integration/EVIDENCE_REUSE.json',
    representation=old['representation'],rows=old['rows'],
    pointwise_causal_boundary=old['pointwise_causal_boundary'],
    new_full_core_campaign_this_loop=False)
 elif model.startswith('SPG_'):
  p=d/'spg/grid257'/model/'AME2408_candidate/RECEIPT.json'
  r=json.loads(p.read_text())
  row.update(generation_status='EXPLICIT_SPHERICAL_CANDIDATE_REOPTIMIZED',
    original_surface_calibration_equivalence='UNRESOLVED',
    original_surface_calibration_equivalence_established=False,whole_original_crust_complete=False,
    prescriptions=5,rows_per_prescription=257,additional_independent_grids=[129],
    actual_common_density_reoptimizations_per_prescription=129,
    baseline_inner_support_fm3=[r['drip']['n_B'],r['transition']['n_cluster_fm3']],
    endpoint_construction=r['transition']['construction'],
    source_transition_prescription='Same-density inhomogeneous versus beta-homogeneous energy-density comparison, Scurto2024 SecII.2 afterEq25.',
    source_transition_prescription_identical=False,
    fixed_density_energy_crossing_fm3=r['transition']['fixed_density_energy_crossing'],
    common_tangent_uniform_endpoint_fm3=r['transition']['n_uniform_fm3'],
    evidence=str(p.relative_to(root)),composition='Free nucleons, leptons, cluster fraction, A and Z',
    partial_crust_variant={'status':'EXECUTED_CHANGED_REST_ENERGY_FUNCTIONAL',
      'upper_density_fm3':.01,'rows':741,'original_model_equivalence':False,
      'evidence':'data/research_finish_20261116/crust_assembly/restored_spg/'+model+'/RECEIPT.json'},
    remaining=['Original surface parameters/fit protocol and low-Yp source equivalence.',
      'Original-model lower/inner reference consistency and full crust; changed-rest candidate is separate.'])
 else:
  raw=d/'geometry/accelerated_final'/model/'SUMMARY.json'
  exported=d/'geometry/exported'/model/'SUMMARY.json'
  row.update(generation_status='SPATIAL_TF_ROUTE_AVAILABLE_REFINEMENT_IN_PROGRESS' if not exported.exists() else 'SPATIAL_TF_REFINEMENT_EXPORTED',
    whole_original_crust_complete=False,
    prior_sampled_inner_support_fm3=[.002,.1],
    composition='Actual spatial n/p/e profiles and totals; free/discrete-cluster native masks unavailable for nonuniform or ambiguous rows.',
    radial_meshes=[192,384],raw_refinement_evidence=str(raw.relative_to(root)),
    corrected_export_evidence=str(exported.relative_to(root)),
    remaining=['Finite geometry/mesh ordering is not a continuum transition or global minimum.',
               'Complete source/reference-consistent crust including supported upper branch and lower piece.'])
  if raw.exists():
   r=json.loads(raw.read_text())
   row['refinement_summary']={k:r.get(k) for k in ('requested_density_nodes','fine_selected_nodes','fine_held_nodes')}
  if model=='GPPVA_NL3wrL55':
   row['upper_extension']={'new_mesh_density_nodes':4,'repeated_anchor_rows':2,
     'densities_fm3':[.0875,.09],'meshes':[192,384],
     'selected_morphology_bracket_fm3':[.085,.0875],
     'interpretation':'Tube selected at.085; uniform at.0875/.09 on both meshes. Finite profile ordering, not a certified energy-crossing root.',
     'evidence':'data/research_finish_20261116/geometry/nl3_upper_exported/RECEIPT.json'}
  partial=d/'crust_assembly/assembled_tf'/model/'RECEIPT.json'
  if partial.exists():
   r=json.loads(partial.read_text())
   row['partial_crust']={'status':'EXECUTED_WITH_ANALYTIC_LOWER_APPROXIMATION',
     'upper_density_fm3':r['inner_cutoff_fm3'],'rows':r['row_count'],
     'evidence':str(partial.relative_to(root)),'whole_crust_complete':False}
  else:
   row['partial_crust']={'status':'UNCHANGED_FUNCTIONAL_INTERFACE_UNRESOLVED_IN_SAMPLED_SUPPORT',
     'evidence':'data/research_finish_20261116/crust_assembly/extended_original_gap_scan.json',
     'scope':'No sampled zero is not a global impossibility proof.'}
 rows.append(row)
result={'schema':'ns-cowling.independent-generation-coverage.v3',
 'base_commit':'573a31485bc9af30915e01fe74f259e4fa5ac54f','status':'PARTIALLY_COMPLETE',
 'objective':'Original five core and thirteen crust models independently regenerable with composition.',
 'declared_executable_route_count':18,'route_breakdown':{'core':5,'candidate_SPG':5,'spatial_TF':8},
 'count_scope':'Includes declared reconstructions/candidates; not eighteen exact original microscopic models or thirteen complete original crusts.',
 'policy':{'CompOSE_agreement':'diagnostic_only','acausality':'density_resolved_support','failed_numerics':'retain_affected_density_and_family','source_equivalence':'separate_from_candidate_execution'},
 'common_lower_crust':{'outer_evidence':'data/research_completion_20261116/outer_crust/generated/RECEIPT.json',
   'analytic_BSk_evidence':'data/research_completion_20261116/bsk22_reconstruction/generated/RECEIPT.json',
   'full_ETFSI_solved':False,'composition_thermodynamic_mu_residual_MeV':[-.0215,.0255],
   'partial_table_scope':'Four unchanged TF and five changed-rest SPG partial tables; local interface certificates for declared approximations.'},
 'models':rows}
assert len(rows)==18 and len({r['model_id'] for r in rows})==18
(d/'MODEL_COVERAGE.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print('Wrote18-model coverage; geometry exported models:',sum((d/'geometry/exported'/r['model_id']/'SUMMARY.json').exists() for r in rows))
