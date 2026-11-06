#!/usr/bin/env python3
from pathlib import Path
import argparse,csv,json,sys,collections
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'compose_eos_pkg'))
from compose_eos import read_compose,check_consistency
from compose_eos.causality import causality_profile

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,default=ROOT/'eos_data_normalized');ap.add_argument('--matched-root',type=Path,required=True);ap.add_argument('--output-json',type=Path,required=True);ap.add_argument('--output-csv',type=Path)
 a=ap.parse_args(); names=['APR_APR','BBB_BHF2','GM_GM1','RG_KDE0v1','RG_SLY4']
 cores={n:read_compose(str(a.source_root/n)) for n in names}; coreprof={n:causality_profile(d) for n,d in cores.items()}
 crusts={n:read_compose(str(a.source_root/n)) for n in sorted(p.name for p in a.source_root.iterdir() if p.is_dir()) if (a.source_root/n).exists() and n not in names}
 crustprof={n:causality_profile(d) for n,d in crusts.items()}
 rows=[]
 for folder in sorted(p for p in a.matched_root.iterdir() if p.is_dir()):
  try:
   meta=json.loads((folder/'MATCHING_METADATA.json').read_text()); eos=read_compose(str(folder)); rep=check_consistency(eos,verbose=False)
   core=meta['core']; crust=meta['crust']; cp=coreprof[core]; cr=crustprof.get(crust,{'classification':{'robust_acausality':False}})
   n2=float(meta['nb_transition_right_fm3']); cd=cores[core]; mi=int(__import__('numpy').searchsorted(eos.nb,n2)); ci=int(__import__('numpy').searchsorted(cd.nb,n2))
   import numpy as np
   tail_equal=bool(len(eos.nb)-mi==len(cd.nb)-ci and len(eos.nb)>mi and np.allclose(eos.nb[mi:],cd.nb[ci:],rtol=2e-12,atol=1e-14) and np.allclose(eos.p[mi:],cd.p[ci:],rtol=2e-10,atol=1e-12) and np.allclose(eos.eps[mi:],cd.eps[ci:],rtol=2e-10,atol=1e-12))
   if rep['overall']['pass_']: kind='PASS'
   elif rep['C10']['pass_'] is False and cp['classification']['robust_acausality'] and tail_equal: kind='SOURCE_CORE_ACAUSALITY_PRESERVED'
   elif rep['C10']['pass_'] is False and cr['classification']['robust_acausality']: kind='SOURCE_CRUST_ACAUSALITY'
   elif rep['C10']['pass_'] is False: kind='MATCHING_OR_NUMERICAL_CANDIDATE'
   else: kind='NON_C10_GATE_FAILURE'
   rows.append({'directory':folder.name,'core':core,'crust':crust,'method':meta['method'],'status':rep['overall']['status'],'failed_checks':';'.join(rep['overall']['failed_checks']),'root_cause':kind,'source_core_robust_acausal':cp['classification']['robust_acausality'],'source_crust_robust_acausal':cr.get('classification',{}).get('robust_acausality',False),'matched_core_tail_byte_semantic_equal':tail_equal,'source_core_first_acausal_nb':cp['pchip']['first_acausal_nb'],'source_core_max_cs2_pchip':cp['pchip']['max_cs2']})
  except Exception as e: rows.append({'directory':folder.name,'root_cause':'AUDIT_READ_ERROR','error':f'{type(e).__name__}: {e}'})
 counts=dict(collections.Counter(r.get('root_cause') for r in rows)); result={'schema':'compose-eos-failure-cause-classification-v1','n_results':len(rows),'counts':counts,'source_core_profiles':coreprof,'source_crust_profiles':crustprof,'results':rows}
 a.output_json.parent.mkdir(parents=True,exist_ok=True);a.output_json.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 if a.output_csv:
  keys=sorted({k for row in rows for k in row});
  with a.output_csv.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 print(json.dumps({'n_results':len(rows),'counts':counts},indent=2))
 if len(rows)!=390 or 'AUDIT_READ_ERROR' in counts: return 2
 return 0
if __name__=='__main__':raise SystemExit(main())
