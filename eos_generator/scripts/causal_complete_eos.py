#!/usr/bin/env python3
"""Create an explicitly model-dependent CSS causal tail and audit the result."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'compose_eos_pkg'))
from compose_eos import read_compose,write_compose,check_consistency
from compose_eos.causality import causal_completion,causality_profile
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True);p.add_argument('--cs2',type=float,default=1.0);p.add_argument('--n-max',type=float);p.add_argument('--n-min-diagnostic',type=float,default=.08);a=p.parse_args()
source=read_compose(a.input);e,meta=causal_completion(source,cs2=a.cs2,n_max=a.n_max);write_compose(a.out,e)
reloaded=read_compose(a.out);gate=check_consistency(reloaded,verbose=False);profile=causality_profile(reloaded,n_min=a.n_min_diagnostic)
robust=all(profile[k]['max_cs2'] is not None and profile[k]['max_cs2']<1.0 for k in ('pchip','cubic_spline','centered_fd'))
meta['post_write_consistency']={'overall':gate['overall'],'C2':gate['C2'],'C10':gate['C10'],'C11':gate['C11']}
meta['independent_causality_profile']=profile
meta['diagnostic_only']=bool(not gate['overall']['pass_'] or not robust)
meta['release_policy']='Promote only when the full consistency gate passes and all three derivative schemes remain causal; otherwise keep diagnostic_only.'
(Path(a.out)/'eos.causal_completion.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
print(json.dumps({'output':a.out,'overall':gate['overall']['status'],'failed_checks':gate['overall']['failed_checks'],'three_scheme_causal':robust,'diagnostic_only':meta['diagnostic_only']},indent=2))
