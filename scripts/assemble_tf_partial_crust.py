#!/usr/bin/env python3
"""Export finite-grid partial crusts through n=.03 using actual TF seam roots.

The high-density pasta and crust endpoint are intentionally separate pieces.
The original TF bulk reference is unchanged; analytic BSk22 is a declared
pressure-integrated approximation, whose composition chemical residual persists.
"""
from pathlib import Path
import json,gzip,hashlib
import numpy as np
from compose_eos.bsk22_thermodynamic_reconstruction import BoundaryAnchoredBSk22
from compose_eos.crust_assembly import assemble_crust
from compose_eos.api import read_compose,write_compose
from compose_eos.tf_composition import read_tf_composition
from build_gppva_tf import export,dump

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/research_finish_20261116/crust_assembly'
def main():
    outerpath=ROOT/'data/research_completion_20261116/outer_crust/generated'
    d=json.loads((outerpath/'RECEIPT.json').read_text())['drip']
    bsk=BoundaryAnchoredBSk22(outer_density=d['n_b'],outer_pressure=d['pressure'],
        outer_mu_b=d['g'],outer_epsilon=d['epsilon'])
    outer=read_compose(outerpath/'compose')
    code_paths=['scripts/assemble_tf_partial_crust.py','scripts/build_gppva_tf.py','python_frontend/compose_eos/crust_assembly.py','python_frontend/compose_eos/bsk22_thermodynamic_reconstruction.py','python_frontend/compose_eos/bsk22_analytic.py','python_frontend/compose_eos/api.py','python_frontend/compose_eos/native.py','python_frontend/compose_eos/eos_data.py','python_frontend/compose_eos/tf_composition.py']
    code_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in code_paths}
    results=[]
    for model in ('GPPVA_DD2','GPPVA_TW','GPPVA_FSU2R','GPPVA_FSU2H'):
        out=BASE/'assembled_tf'/model;out.mkdir(parents=True,exist_ok=True)
        ipath=BASE/f'tf_interfaces/{model}/INTERFACE.json.gz'
        seam=json.loads(gzip.decompress(ipath.read_bytes()))
        if seam['status']!='CANDIDATE_COMMON_TANGENT_REOPTIMIZED':
            raise RuntimeError('unresolved seam '+model)
        endpoint=seam['tf_endpoint'];n=endpoint['n_b']
        oldpath=ROOT/f'data/research_completion_20261116/gppva_tf/final/{model}/selected_states.json.gz'
        old=json.loads(gzip.decompress(oldpath.read_bytes()))
        states=[endpoint]+[r for r in old if n<r['n_b']<=.03]
        if any(r.get('geometry','droplet')!='droplet' for r in states):
            raise RuntimeError('upper branch not spherical')
        innerdir=out/'inner_source';innerdir.mkdir(exist_ok=True)
        dump(innerdir/'selected_states.json.gz',states)
        innerreceipt=export(model,states,innerdir)
        inner=read_compose(innerdir/'compose')
        composition=read_tf_composition(innerdir)
        lower,bsrows=bsk.eos_data(np.geomspace(bsk.n0,seam['bsk_endpoint']['n_b'],129))
        merged,manifest=assemble_crust([outer,lower,inner],
            labels=['independent_outer','analytic_BSk22','unchanged_spherical_TF_partial'],
            mn=939.56542052,mp=938.27208816,pressure_rtol=1e-7,mu_atol_MeV=1e-6)
        write_compose(out/'compose',merged);back=read_compose(out/'compose')
        for key in ('nb','p','eps','mub','Ye','Yn','Yp','Yh','Aav','Zav'):
            if not np.allclose(getattr(back,key),getattr(merged,key),rtol=2e-12,atol=2e-12):
                raise RuntimeError('readback changed '+key)
        for key in ('Iphase','Yn_present','Yp_present','Ye_present','Yh_present','row_status'):
            if not np.array_equal(getattr(back,key),getattr(merged,key)):
                raise RuntimeError('readback changed '+key)
        manifest.update(required_companions=['COMPOSITION.json','inner_source/composition.json','inner_source/selected_states.json.gz'],
            physical_scope='PARTIAL_CRUST_THROUGH_0.03_NOT_COMPLETE_CRUST',
            phase_interpolation='No interpolation across either declared density jump; upper branch is finite mixed-mesh spherical sample.')
        dump(out/'ASSEMBLY.json',manifest)
        # Explicit global row mapping prevents TF total fractions from being
        # misread as free-neutron/proton fractions in CompOSE.
        comp={'row_count':back.n,'total_composition_available_every_row':True,
          'free_cluster_decomposition':'unavailable in nonuniform spatial TF; masks preserved',
          'row_segments':[
            {'start':0,'stop':outer.n,'representation':'discrete_nuclear_cell','carrier':'compose/eos.compo'},
            {'start':outer.n,'stop':outer.n+lower.n,'representation':'analytic_cluster_and_gas','carrier':'compose/eos.compo'},
            {'start':outer.n+lower.n,'stop':back.n,'representation':'spatial_TF_total_species',
             'carrier':'inner_source/composition.json','profiles':'inner_source/selected_states.json.gz'}]}
        dump(out/'COMPOSITION.json',comp)
        receipt=dict(model=model,status='PARTIAL_CRUST_NATIVE_READBACK_AND_PROFILE_READER_VERIFIED',
            row_count=back.n,whole_crust_complete=False,inner_cutoff_fm3=.03,
            source_functional='unchanged RMF functional/reference',
            approximation='independent outer + pressure-integrated analytic BSk22 + spherical TF finite-grid states',
            bsk_composition_mu_residual_range_MeV=[min(r['source_composition_mu_n_residual_MeV'] for r in bsrows),
                max(r['source_composition_mu_n_residual_MeV'] for r in bsrows)],
            seam=seam['status'],mesh_convergence_claim=False,global_minimum_claim=False,
            tf_rows=inner.n,tf_mesh_points=sorted(set(s['points'] for s in states)),
            assembly_code_sha256=code_hashes,
            seam_execution_code_sha256=json.loads((BASE/'tf_interfaces/RECEIPT.json').read_text())['code_sha256'],
            source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ipath,oldpath,outerpath/'RECEIPT.json',*sorted((outerpath/'compose').glob('*'))] if p.is_file()},
            inner_export=innerreceipt)
        receipt['output_sha256']={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file() and p.name!='RECEIPT.json'}
        dump(out/'RECEIPT.json',receipt)
        results.append({k:v for k,v in receipt.items() if k not in ('inner_export','source_sha256')})
        print(model,back.n,'readback and TF composition reader verified',flush=True)
    dump(BASE/'assembled_tf/CAMPAIGN.json',{'results':results})
if __name__=='__main__':main()
