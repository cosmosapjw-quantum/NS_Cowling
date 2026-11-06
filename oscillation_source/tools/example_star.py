"""Compute one common polytropic star and two identified fluid modes."""
from pathlib import Path
import argparse,json,sys
if __package__ in (None,''):sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.runtime import CaseRunner
from tools.background_io import write_background
ROOT=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--workspace',choices=['learner','solutions'],default='learner');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--dq',type=float,default=.01);a=ap.parse_args()
 if a.output.exists():ap.error('Output directory exists')
 if not 0<a.dq<=.1:ap.error('Require 0 < dq <= 0.1')
 a.output.mkdir(parents=True);result={'status':'incomplete','workspace':a.workspace}
 try:
  runner=CaseRunner(ROOT,a.output/'runs',a.workspace);art=runner.build_fpm()
  c=299792458.;en=1.5e-10;na=1e44
  params={'n_anchor_m3':na,'p_anchor_pa':1e5*((en/c**2)*na/1000)**2/10,
   'gamma':2.,'c_j':en,'n_min_m3':1e30,'n_max_m3':1e47,'pc_pa':6.4e33,'p_cut_pa':6.4e21,'r0_m':1.,
   'rtol':1e-11,'atol_r_m':1e-7,'atol_mass_kg':1e14,'atol_nu':1e-13,'dq_max':a.dq,'emit_profile':True,'max_steps':200000}
  star=runner.run_artifact(art,'tov_baseline',params,timeout=60)
  if star['status']!='ok':raise RuntimeError(star['status']+': '+star.get('detail',''))
  write_background(star,a.output/'background');result.update(inputs=params,mass_kg=star['mass_kg'],radius_m=star['radius_m'],modes={})
  for label,lo,hi in [('f',1800.,2100.),('p1',4100.,4600.)]:
   mode=runner.run_artifact(art,'mode_refine',{'background_path':str((a.output/'background/background.csv').resolve()),'frequency_min_hz':lo,'frequency_max_hz':hi,'ode_rtol':1e-12,'root_xtol':1e-11,'root_rtol':0.,'residual_tol':1e-7,'emit_profile':True,'max_steps':500000},timeout=90)
   if mode['status']!='ok':raise RuntimeError(label+': '+mode['status']+': '+mode.get('detail',''))
   result['modes'][label]={k:mode[k] for k in ['mode_hz','mode_label','mode_identity','mode_residual','mode_bracket_width_hz']}
  result['status']='ok';code=0
 except Exception as e:result['error']=type(e).__name__+': '+str(e);code=1
 (a.output/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,ensure_ascii=False,indent=2));return code
if __name__=='__main__':raise SystemExit(main())
