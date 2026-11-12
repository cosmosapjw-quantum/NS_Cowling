"""M7: Native manifest-driven 390 batch orchestration - manifest generator.

Generates the exact pair/method manifest for the EOS390 regeneration:
  13 crust_models x 5 core_models x 6 methods = 390 unique tasks.

Each task is a (crust, core, method) triple with:
  - a unique task_id (deterministic hash of the triple)
  - provenance (source crust/core EOS paths, method physics reference)
  - an isolated output folder path
  - a status slot (PENDING / RUNNING / PASS / FAIL / SKIPPED_BLOCKED)

The manifest is consumed by the batch orchestrator (m7_batch_orchestrator.py)
which dispatches each task to the calculation backend (Fortran in production,
Python oracle in this sandbox) and writes a per-task machine receipt.

This is the Python reproduction of the M7 acceptance criterion "Fortran
consumes one exact pair/method manifest" - the manifest format and
cardinality (exactly 390) are the production contract.
"""
from __future__ import annotations
import json
import hashlib
import os
from pathlib import Path
from itertools import product


# The 390 matrix (from the DAG objective.matrix)
CRUST_MODELS = [
    'SPG_M1', 'SPG_M2', 'SPG_M3', 'SPG_M4', 'SPG_M5',
    'GPPVA_DD2', 'GPPVA_DDME2', 'GPPVA_FSU2',
    'GPPVA_FSU2H', 'GPPVA_FSU2R', 'GPPVA_NL3wrL55',
    'GPPVA_TM1e', 'GPPVA_TW',
]
CORE_MODELS = ['APR_APR', 'BBB_BHF2', 'GM_GM1', 'RG_SLY4', 'RG_KDE0v1']
METHODS = ['maxwell', 'css', 'cubic_F', 'cubic_P', 'quintic_F', 'quintic_P']

# Physics reference for each matching method (provenance)
METHOD_PHYSICS = {
    'maxwell': {
        'kind': 'sharp_first_order',
        'description': 'Sharp Maxwell construction: equal P and mu_B coexistence with density jump',
        'reference': 'Glendenning 1991 PhysRevD.46.1274 (relativistic two-phase); Maxwell 1875',
        'acceptance_node': 'M5',
    },
    'css': {
        'kind': 'sharp_first_order_high_density_branch',
        'description': 'Constant-Sound-Speed high-density branch with sharp interface jump',
        'reference': 'Alford, Han & Prakash 2013 PRD 88 083013 (arXiv:1302.4732); 2015 PRD 92 083002 (arXiv:1501.07902)',
        'acceptance_node': 'M5A',
    },
    'cubic_F': {
        'kind': 'smooth_crossover',
        'description': 'Cubic Hermite interpolation in free energy F = eps/n_B (continuous, no jump)',
        'reference': 'Legacy NPPHB Hermite bridge (project-defined; not a literature attribution per legacy_method_semantics)',
        'acceptance_node': 'E390-02C',
    },
    'cubic_P': {
        'kind': 'smooth_crossover',
        'description': 'Cubic Hermite interpolation in pressure P (continuous, no jump)',
        'reference': 'Legacy NPPHB Hermite bridge (project-defined)',
        'acceptance_node': 'E390-02C',
    },
    'quintic_F': {
        'kind': 'smooth_crossover',
        'description': 'Quintic Hermite interpolation in free energy F (continuous, no jump, C2 smooth)',
        'reference': 'Legacy NPPHB Hermite bridge (project-defined)',
        'acceptance_node': 'E390-02C',
    },
    'quintic_P': {
        'kind': 'smooth_crossover',
        'description': 'Quintic Hermite interpolation in pressure P (continuous, no jump, C2 smooth)',
        'reference': 'Legacy NPPHB Hermite bridge (project-defined)',
        'acceptance_node': 'E390-02C',
    },
}


def task_id(crust: str, core: str, method: str) -> str:
    """Deterministic 12-char task id from the (crust, core, method) triple.
    Uses SHA-256 of the canonical string, truncated. The id is unique per triple
    and stable across runs (so a re-run resumes by id)."""
    canon = f'{crust}|{core}|{method}'
    return hashlib.sha256(canon.encode('utf-8')).hexdigest()[:12]


def generate_manifest(*, output_root: str, crust_eos_root: str = None,
                      core_eos_root: str = None,
                      available_crust: list = None,
                      available_core: list = None) -> dict:
    """Generate the 390-task manifest.

    output_root:    root directory under which each task gets an isolated folder
    crust_eos_root: where the crust EOS CompOSE tables live (per-model subdirs)
    core_eos_root:  where the core EOS CompOSE tables live (per-model subdirs)
    available_crust: subset of CRUST_MODELS actually present (others BLOCKED_MISSING_SOURCE)
    available_core:  subset of CORE_MODELS actually present (others BLOCKED_MISSING_SOURCE)

    Returns the manifest dict. The manifest has exactly 390 tasks; unavailable
    sources are marked BLOCKED_MISSING_SOURCE but still occupy their slot
    (no silent skipping - M7 acceptance: "Failures remain failures/diagnostic-only").
    """
    ac = set(available_crust or CRUST_MODELS)
    aco = set(available_core or CORE_MODELS)
    tasks = []
    seen_ids = set()
    seen_triples = set()
    # Deterministic ordering: crust (outer), core (middle), method (inner)
    for crust, core, method in product(CRUST_MODELS, CORE_MODELS, METHODS):
        triple = (crust, core, method)
        triple_str = f'{crust}|{core}|{method}'
        if triple_str in seen_triples:
            raise AssertionError(f'duplicate triple detected: {triple} - manifest must be unique')
        seen_triples.add(triple_str)
        tid = task_id(crust, core, method)
        if tid in seen_ids:
            raise AssertionError(f'task id collision: {tid} for {triple} - hash must be unique')
        seen_ids.add(tid)
        # Status: PENDING if both sources available, else BLOCKED_MISSING_SOURCE
        blocked = []
        if crust not in ac:
            blocked.append(f'crust:{crust}')
        if core not in aco:
            blocked.append(f'core:{core}')
        status = 'PENDING' if not blocked else 'BLOCKED_MISSING_SOURCE'
        task = {
            'task_id': tid,
            'crust_model': crust,
            'core_model': core,
            'method': method,
            'method_kind': METHOD_PHYSICS[method]['kind'],
            'method_reference': METHOD_PHYSICS[method]['reference'],
            'acceptance_node': METHOD_PHYSICS[method]['acceptance_node'],
            'crust_eos_path': f'{crust_eos_root}/{crust}' if crust_eos_root else f'inputs/crust/{crust}',
            'core_eos_path': f'{core_eos_root}/{core}' if core_eos_root else f'inputs/core/{core}',
            'output_folder': f'{output_root}/{tid}_{crust}_{core}_{method}',
            'status': status,
            'blocked_reason': ';'.join(blocked) if blocked else None,
            'provenance': {
                'matrix_cardinality': '13 crust x 5 core x 6 method = 390',
                'task_index_global': len(tasks),
                'canonical_triple': triple_str,
            },
        }
        tasks.append(task)
    # Cardinality check (M7 acceptance: "Exactly 390 unique tasks are scheduled")
    if len(tasks) != 390:
        raise AssertionError(f'cardinality violation: generated {len(tasks)} tasks, expected 390')
    if len(seen_triples) != 390:
        raise AssertionError(f'uniqueness violation: {len(seen_triples)} unique triples, expected 390')
    if len(seen_ids) != 390:
        raise AssertionError(f'id uniqueness violation: {len(seen_ids)} unique ids, expected 390')
    # Sanity: cross product cardinality
    if len(CRUST_MODELS) * len(CORE_MODELS) * len(METHODS) != 390:
        raise AssertionError('matrix cardinality drift: 13x5x6 != 390 - update CRUST/CORE/METHOD')
    manifest = {
        'schema': 'ns-cowling.eos390-m7-batch-manifest',
        'schema_version': '1.0.0',
        'matrix': {
            'crust_models': CRUST_MODELS,
            'core_models': CORE_MODELS,
            'methods': METHODS,
            'cardinality_check': f'{len(CRUST_MODELS)} x {len(CORE_MODELS)} x {len(METHODS)} = {len(CRUST_MODELS)*len(CORE_MODELS)*len(METHODS)}',
            'expected_pairs': 390,
            'actual_tasks': len(tasks),
        },
        'available_sources': {
            'crust': sorted(ac),
            'core': sorted(aco),
            'crust_coverage': f'{len(ac & set(CRUST_MODELS))}/13',
            'core_coverage': f'{len(aco & set(CORE_MODELS))}/5',
        },
        'tasks': tasks,
        'task_count': len(tasks),
        'unique_triple_count': len(seen_triples),
        'unique_id_count': len(seen_ids),
        'status_summary': {
            'PENDING': sum(1 for t in tasks if t['status'] == 'PENDING'),
            'BLOCKED_MISSING_SOURCE': sum(1 for t in tasks if t['status'] == 'BLOCKED_MISSING_SOURCE'),
        },
    }
    return manifest


def write_manifest(manifest: dict, path: str) -> None:
    """Write the manifest atomically (tmp + rename)."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, allow_nan=False)
        f.write('\n')
    os.replace(tmp, p)


if __name__ == '__main__':
    import sys
    # Default: full 390 manifest with the 6 crust models we actually have
    available_crust = ['SPG_M1', 'SPG_M3', 'SPG_M5', 'GPPVA_DD2', 'GPPVA_DDME2', 'GPPVA_TW']
    # Core models: synthetic polytropes stand in for the real APR/BBB/GM/RG models
    # (the real cores would come from the M0 replacement ZIP)
    available_core = ['APR_APR', 'BBB_BHF2', 'GM_GM1', 'RG_SLY4', 'RG_KDE0v1']
    out_root = '/home/z/my-project/download/m7_batch/outputs'
    manifest = generate_manifest(
        output_root=out_root,
        crust_eos_root='/home/z/my-project/download/m4_spliced',
        core_eos_root='/home/z/my-project/download/m7_batch/inputs_core_synth',
        available_crust=available_crust,
        available_core=available_core,
    )
    out_path = '/home/z/my-project/download/m7_batch/M7_BATCH_MANIFEST.json'
    write_manifest(manifest, out_path)
    print(f'manifest written: {out_path}')
    print(f'task_count: {manifest["task_count"]}')
    print(f'unique_triples: {manifest["unique_triple_count"]}')
    print(f'unique_ids: {manifest["unique_id_count"]}')
    print(f'status_summary: {manifest["status_summary"]}')
    print(f'crust coverage: {manifest["available_sources"]["crust_coverage"]}')
    print(f'core coverage: {manifest["available_sources"]["core_coverage"]}')
