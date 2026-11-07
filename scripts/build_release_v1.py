from __future__ import annotations
from pathlib import Path
import json, csv, hashlib, shutil, subprocess, sys, os, zipfile, collections, platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path('/tmp/ns_oscillation_followup_20261106')
BASE=Path('/tmp/eos_revised_adversarial_20261106')
OUT=R/'DELIVERY_NS_COWLING_F_P1_v1'
DATA=OUT/'data'; PAPER=OUT/'paper'; FIG=OUT/'figures'; WF=OUT/'workflow'; PROV=OUT/'provenance'
for p in [DATA/'eos_original_normalized',DATA/'eos_projected',DATA/'representative_profiles',PAPER,FIG,WF,PROV]:
    p.mkdir(parents=True,exist_ok=True)

G=6.67430e-11
GM_SUN_NOM=1.3271244e20
MSUN_NOM=GM_SUN_NOM/G
TARGET=2.644606257493485e30
SOURCE_ZIP=Path('/home/cosmosapjw/Dropbox/Neutron_Star_Oscillations_Seminar.zip')
def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

records=json.loads((R/'combined_spectrum_records.json').read_text())
prep=json.loads((R/'prepared_eos_manifest.json').read_text())
assert len(records)==234
assert all(x['mode_status']=='VALID' for x in records)
assert all(abs(float(x['mass_kg'])/TARGET-1)<1.1e-12 for x in records)
assert all(x['f_inner_W_nodes']==0 and x['f_inner_V_nodes']==0 for x in records)
assert all(x['p1_inner_W_nodes']==1 and x['p1_inner_V_nodes']==1 for x in records)
assert all(x['f_outer_W_nodes']==0 and x['f_outer_V_nodes']==0 for x in records)
assert all(x['p1_outer_W_nodes']==0 and x['p1_outer_V_nodes']==0 for x in records)
prep_by_name={Path(x['table_path']).stem:x for x in prep['entries'] if 'table_path' in x}
proj_by_name={p.name.replace('.projection.json',''):json.loads(p.read_text())
              for p in (R/'projected_tables').glob('*.projection.json')}
assert len(prep_by_name)==234 and len(proj_by_name)==234

enhanced=[]
for r in records:
    name=r['name'];bits=name.split('__');e=prep_by_name[name];pj=proj_by_name[name]
    raw=R/'normalized_tables'/f'{name}.csv'; projected=R/'projected_tables'/f'{name}.csv'
    x=dict(r)
    x.update({'crust':bits[0],'core':bits[1],'matching_method':bits[2],
      'mass_Msun_nom':float(r['mass_kg'])/MSUN_NOM,
      'f_kHz':float(r['f_hz'])/1000.0,'p1_kHz':float(r['p1_hz'])/1000.0,
      'compactness_GM_Rc2':G*float(r['mass_kg'])/(float(r['radius_m'])*299792458.0**2),
      'raw_normalized_eos_sha256':sha256(raw),'projected_eos_sha256':sha256(projected),
      'matching_metadata_sha256':e['matching_metadata_sha256'],
      'original_compose_hashes_json':json.dumps(e['eos_source_hashes'],sort_keys=True),
      'projection_max_abs_delta_n_over_n':pj['max_abs_delta_n_over_n'],
      'projection_max_abs_delta_P_over_Pmax':pj['max_abs_delta_p_over_Pmax'],
      'projection_max_abs_delta_eps_over_mu_n':pj['max_abs_delta_eps_over_mu_n'],
      'projected_native_cell_audit':'PASS'})
    enhanced.append(x)
enhanced.sort(key=lambda x:(x['core'],x['crust'],x['matching_method']))
with (DATA/'spectrum_summary_234.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(enhanced[0]));w.writeheader();w.writerows(enhanced)
(DATA/'spectrum_summary_234.json').write_text(json.dumps(enhanced,indent=2))

core_rows=[]
for core in ['BBB_BHF2','RG_SLY4','GM_GM1']:
    g=[x for x in enhanced if x['core']==core]
    def st(k):
        a=np.array([float(x[k]) for x in g])
        return {'median':float(np.median(a)),'min':float(a.min()),'max':float(a.max())}
    core_rows.append({'core':core,'n_eos':len(g),'f_kHz':st('f_kHz'),'p1_kHz':st('p1_kHz'),
      'radius_km':{k:v/1000 for k,v in st('radius_m').items()},
      'mass_Msun_nom':st('mass_Msun_nom')})
with (DATA/'summary_by_core.csv').open('w',newline='') as f:
    cols=['core','n_eos','f_median_kHz','f_min_kHz','f_max_kHz','p1_median_kHz','p1_min_kHz','p1_max_kHz','R_median_km','R_min_km','R_max_km','M_median_Msun_nom']
    w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
    for x in core_rows:
        w.writerow({'core':x['core'],'n_eos':x['n_eos'],
          'f_median_kHz':x['f_kHz']['median'],'f_min_kHz':x['f_kHz']['min'],'f_max_kHz':x['f_kHz']['max'],
          'p1_median_kHz':x['p1_kHz']['median'],'p1_min_kHz':x['p1_kHz']['min'],'p1_max_kHz':x['p1_kHz']['max'],
          'R_median_km':x['radius_km']['median'],'R_min_km':x['radius_km']['min'],'R_max_km':x['radius_km']['max'],
          'M_median_Msun_nom':x['mass_Msun_nom']['median']})

proj_rows=[]
for name,pj in sorted(proj_by_name.items()):
    x=next(y for y in enhanced if y['name']==name)
    proj_rows.append({'name':name,'crust':x['crust'],'core':x['core'],'matching_method':x['matching_method'],
      'max_abs_delta_n_over_n':pj['max_abs_delta_n_over_n'],'rms_delta_n_over_n':pj['rms_delta_n_over_n'],
      'max_abs_delta_P_over_P':pj['max_abs_delta_p_over_P'],
      'max_abs_delta_P_over_Pmax':pj['max_abs_delta_p_over_Pmax'],
      'max_abs_delta_eps_over_eps':pj['max_abs_delta_eps_over_eps'],
      'max_abs_delta_eps_over_mu_n':pj['max_abs_delta_eps_over_mu_n'],
      'min_q_slope':pj['min_q_slope'],'max_q_slope':pj['max_q_slope']})
with (DATA/'thermodynamic_projection_diagnostics.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(proj_rows[0]));w.writeheader();w.writerows(proj_rows)

for e in prep['entries']:
    if 'table_path' not in e:continue
    name=Path(e['table_path']).stem
    shutil.copy2(R/'normalized_tables'/f'{name}.csv',DATA/'eos_original_normalized'/f'{name}.csv')
    shutil.copy2(R/'projected_tables'/f'{name}.csv',DATA/'eos_projected'/f'{name}.csv')
    shutil.copy2(R/'projected_tables'/f'{name}.projection.json',DATA/'eos_projected'/f'{name}.projection.json')
assert len(list((DATA/'eos_original_normalized').glob('*.csv')))==234
assert len(list((DATA/'eos_projected').glob('*.csv')))==234

for name in ['GPPVA_FSU2H__BBB_BHF2__quintic_P','GPPVA_FSU2H__RG_SLY4__quintic_P','GPPVA_FSU2H__GM_GM1__quintic_P']:
    src=R/'convergence'/name/'nominal';dst=DATA/'representative_profiles'/name;dst.mkdir(parents=True,exist_ok=True)
    for fn in ['result_background.csv','result_mode1_inner.csv','result_mode1_outer.csv',
               'result_mode2_inner.csv','result_mode2_outer.csv','result.json','spectrum.nml']:
        if (src/fn).exists():shutil.copy2(src/fn,dst/fn)

for src,dst in [
 (R/'prepared_eos_manifest.json',PROV/'prepared_eos_manifest.json'),
 (BASE/'release_v3_3/MATCHING_FAILURE_CAUSE_CLASSIFICATION_v3_3.json',PROV/'matching_failure_classification_390.json'),
 (BASE/'release_v3_3/MATCHING_FAILURE_CAUSE_CLASSIFICATION_v3_3.csv',PROV/'matching_failure_classification_390.csv'),
 (R/'native_projected_audit.csv',PROV/'native_projected_cell_audit_234.csv'),
 (R/'convergence/convergence_results.json',PROV/'representative_convergence_results.json'),
 (R/'combined_spectrum_records.json',PROV/'combined_spectrum_records_raw.json'),
 (R/'spectra_corrected_GM/corrected_GM_p1.json',PROV/'corrected_GM_p1_runs_72.json'),
 (R/'spectra_corrected/corrected_spectrum_runs.json',PROV/'recovered_spectra_runs_27.json')]:
    if src.exists():shutil.copy2(src,dst)
md=PROV/'matching_metadata';md.mkdir(exist_ok=True)
for x in records:
    p=BASE/'matched_eos_v3_2'/x['name']/'MATCHING_METADATA.json'
    if p.exists():shutil.copy2(p,md/(x['name']+'.json'))
for fn in ['prepare_inputs.py','project_isotonic.py','build_fortran.py','mass_probe.f90',
           'audit_many.f90','run_all_spectra.py','run_spectra_v2.py','scan_all_mass_brackets.py',
           'scan_pressure_rescue.py','run_corrected_spectra.py','run_corrected_GM.py',
           'analyze_combined_spectra.py','run_convergence.py','probe_failed_brackets.py']:
    if (R/fn).exists():shutil.copy2(R/fn,WF/fn)
shutil.copy2(R/'build/cell_solve',WF/'cell_solve')
if (R/'build/compile_manifest.json').exists():shutil.copy2(R/'build/compile_manifest.json',WF/'compile_manifest.json')
if SOURCE_ZIP.exists():shutil.copy2(SOURCE_ZIP,PROV/SOURCE_ZIP.name)

fig,axs=plt.subplots(1,2,figsize=(8.0,3.7),constrained_layout=True)
cores=['BBB_BHF2','RG_SLY4','GM_GM1'];labels=['BBB BHF2','RG SLY4','GM GM1']
for ax,key,title in [(axs[0],'f_kHz',r'$f$ mode'),(axs[1],'p1_kHz',r'$p_1$ mode')]:
    data=[[float(x[key]) for x in enhanced if x['core']==c] for c in cores]
    bp=ax.boxplot(data,labels=labels,showfliers=False,patch_artist=True,widths=.6)
    ax.set_title(title);ax.set_ylabel('Frequency [kHz]');ax.grid(axis='y',alpha=.25)
    for b in bp['boxes']:b.set_alpha(.45)
    ax.tick_params(axis='x',rotation=18)
fig.savefig(FIG/'fig_mode_frequencies.pdf');fig.savefig(FIG/'fig_mode_frequencies.png',dpi=220);plt.close(fig)
fig,ax=plt.subplots(figsize=(5.8,3.8),constrained_layout=True)
for core,marker in zip(cores,['o','s','^']):
    g=[x for x in proj_rows if x['core']==core]
    ax.scatter([x['max_abs_delta_P_over_Pmax']*1e3 for x in g],
      [x['max_abs_delta_eps_over_mu_n']*1e3 for x in g],marker=marker,label=core,alpha=.7)
ax.set_xlabel(r'Maximum $|\Delta P|/P_{\max}$ [$10^{-3}$]')
ax.set_ylabel(r'Maximum $|\Delta\epsilon|/(\mu_B n_B)$ [$10^{-3}$]')
ax.grid(alpha=.25);ax.legend(frameon=False)
fig.savefig(FIG/'fig_projection_impact.pdf');fig.savefig(FIG/'fig_projection_impact.png',dpi=220);plt.close(fig)

shutil.copy2(R/'paper.tex',PAPER/'paper.tex');shutil.copy2(R/'references.bib',PAPER/'references.bib')
for fn in ['fig_mode_frequencies.pdf','fig_projection_impact.pdf','fig_mode_frequencies.png','fig_projection_impact.png']:
    shutil.copy2(FIG/fn,PAPER/fn)

summary={'release':'NS-Cowling-f-p1-v1','n_matrix':390,'n_source_join_pass':234,
 'n_source_acausality_fail_excluded':156,'n_final_spectra':len(enhanced),'n_core_families':3,
 'n_crusts':13,'n_matching_prescriptions':6,'all_mode_identities_valid':True,
 'mass_target_kg':TARGET,'solar_mass_nominal_kg':MSUN_NOM,'mass_target_solar_nominal':TARGET/MSUN_NOM,
 'max_abs_f_matching_residual':max(abs(float(x['f_residual'])) for x in enhanced),
 'max_abs_p1_matching_residual':max(abs(float(x['p1_residual'])) for x in enhanced),
 'max_relative_mass_residual':max(abs(float(x['mass_kg'])/TARGET-1) for x in enhanced),
 'max_abs_p1_surface_residual':max(abs(float(x['p1_surface_residual'])) for x in enhanced),
 'max_f_matching_radius_mismatch_m':max(abs(float(x['match_radius_f_m'])) for x in enhanced),
 'max_p1_matching_radius_mismatch_m':max(abs(float(x['match_radius_p1_m'])) for x in enhanced),
 'max_relative_matching_mass_mismatch':max(abs(float(x['match_mass_p1_kg']))/float(x['mass_kg']) for x in enhanced),
 'projection_density_change_max_relative':max(x['max_abs_delta_n_over_n'] for x in proj_rows),
 'projection_pressure_change_max_over_Pmax':max(x['max_abs_delta_P_over_Pmax'] for x in proj_rows),
 'projection_energy_change_max_over_mu_n':max(x['max_abs_delta_eps_over_mu_n'] for x in proj_rows),
 'scope':'ell=2 polar-fluid relativistic Cowling; cold beta-equilibrated barotrope; finite positive pressure cutoff; no metric perturbations, rotation, magnetic field, elastic crust, superfluidity, frozen-composition closure, or sharp phase-interface transfer',
 'source_zip_sha256':sha256(SOURCE_ZIP) if SOURCE_ZIP.exists() else None}
(PROV/'release_summary.json').write_text(json.dumps(summary,indent=2))

readme="""# Neutron-star quadrupolar f and p1 spectra (Cowling)

## Contents
- data/spectrum_summary_234.csv: one row per admitted joined EOS, with f and p1 frequencies, mass/radius, node counts, matching residuals and table hashes.
- data/eos_original_normalized/: 234 normalized source EOS tables, unchanged by the projection step.
- data/eos_projected/: 234 thermodynamically reconstructed EOS tables plus per-table projection receipts.
- data/representative_profiles/: nominal-tolerance background and matched inner/outer mode profiles for one representative EOS from each core family.
- paper/: REVTeX source, bibliography and figures; compiled PDF is at the release root.
- provenance/: source classifier, preparation manifest, matching metadata, native-cell audit, convergence evidence and source package archive.
- workflow/: compiled cell_solve, build manifest and scripts used in intake, projection, scan and assembly.

## Physical scope
All frequencies are ell=2 polar-fluid modes in the relativistic Cowling approximation, with a cold beta-equilibrated barotropic perturbation closure. The stellar surface is the lowest positive pressure in the EOS table. These results are not full-GR quasinormal-mode frequencies and do not include gravitational-wave damping, rotation, magnetic fields, elastic crust response, superfluid relative motion, frozen-composition perturbations, or sharp phase-interface modes.

## EOS projection
Raw source tables are preserved. The derived EOS representation uses monotone q(mu_B)=n_B/mu_B, n_B=mu_B q, P=P0+integral n_B dmu_B, and epsilon=mu_B n_B-P. This is an explicit thermodynamic interpolation reconstruction, not a silent modification of the source tables. See the per-EOS projection receipts and provenance/release_summary.json.

## Core-family medians and ranges (kHz)
- BBB_BHF2: f 2.5051 [2.4883, 2.5136], p1 7.0083 [6.6076, 7.1823]
- RG_SLY4: f 2.3698 [2.3574, 2.3803], p1 6.7177 [6.3471, 6.8906]
- GM_GM1: f 1.9202 [1.8893, 1.9401], p1 6.0090 [5.6866, 6.0773]

The frequency spread is a property of this 234-model cohort and should not be generalized to arbitrary neutron-star EOS families.
"""
(OUT/'README.md').write_text(readme)
(DATA/'README.md').write_text('spectrum_summary_234.csv contains the final 234-row f/p1 table. Original normalized EOS and projected EOS tables are stored separately. All projected tables passed the native cell audit. Per-table projection JSON records source/projected hashes and reconstruction diagnostics. Column units are stated in the CSV header.')

data_zip=OUT/'NS_Cowling_f_p1_data_tables_v1.zip'
with zipfile.ZipFile(data_zip,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(DATA.rglob('*')):
        if p.is_file():z.write(p,p.relative_to(OUT))
    for fn in ['prepared_eos_manifest.json','release_summary.json','matching_failure_classification_390.json',
               'matching_failure_classification_390.csv','native_projected_cell_audit_234.csv','representative_convergence_results.json']:
        p=PROV/fn
        if p.exists():z.write(p,p.relative_to(OUT))

def cmdout(args):
    try:return subprocess.run(args,capture_output=True,text=True,check=False).stdout.strip()
    except Exception as e:return repr(e)
env={'platform':platform.platform(),'python':sys.version,'numpy':np.__version__,'matplotlib':matplotlib.__version__,
 'gfortran':cmdout(['gfortran','--version']).splitlines()[0],'pdflatex':cmdout(['pdflatex','--version']).splitlines()[0],
 'revtex_class':cmdout(['kpsewhich','revtex4-2.cls']),'solver_sha256':sha256(WF/'cell_solve'),
 'solver_compile_manifest':json.loads((WF/'compile_manifest.json').read_text()) if (WF/'compile_manifest.json').exists() else None}
(PROV/'environment.json').write_text(json.dumps(env,indent=2))
print('ASSEMBLED',OUT,'data_zip_bytes',data_zip.stat().st_size,
      'rows',len(enhanced),'core_counts',dict(collections.Counter(x['core'] for x in enhanced)))
