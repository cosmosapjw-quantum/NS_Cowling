#!/usr/bin/env python3
"""Run the user's actual R17 DOP853/GL4 engine on the course's frozen native EOS.

The original repository must be supplied explicitly; no production code or saved
frequency is substituted into the course solver. This optional comparison needs
SciPy and the R17 repository. The ordinary course calculation does not.
"""
from __future__ import annotations
import argparse,importlib,json,math,sys,time,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from tools.native_eos_reference import NativeEOS,CONV,G
from tools.atomic_io import write_json,digest

def run_reference(production:Path,output:Path,*,alpha:float=0.,refine:int=4,
                  rtol:float=2e-13,fit:float=.5)->dict:
    production=Path(production).resolve();output=Path(output).resolve()
    if refine not in (1,2,4,8,16) or not 0<rtol<=1e-6 or not .05<fit<.95 or not math.isfinite(alpha):
        raise ValueError('Invalid comparison controls')
    for name in ['r16_model.py','r16_shooting.py','r16_transfer.c','eos_segment.py']:
        if not (production/'tools'/name).is_file():raise FileNotFoundError('Required R17 source: '+name)
    if output.is_relative_to(ROOT) or output.is_relative_to(production):
        raise ValueError('Choose an output outside both source trees')
    output.mkdir(parents=True,exist_ok=False)
    try:
        from scipy.optimize import brentq
        sys.path.insert(0,str(production/'tools'))
        bgmodule=importlib.import_module('r16_model');shoot=importlib.import_module('r16_shooting')
        if not Path(bgmodule.__file__).resolve().is_relative_to(production) or not Path(shoot.__file__).resolve().is_relative_to(production):
            raise ValueError('A different production module was already imported; use a fresh Python process')
        source={str(p.relative_to(production)):digest(p) for p in sorted((production/'tools').glob('*')) if p.suffix in ['.py','.c']}
        tab=ROOT/'data/normalized/compose_69.csv'
        receipt={'production_root':str(production),'source_sha256':source,'table_sha256':digest(tab),
                 'driver_sha256':digest(__file__),'adapter_sha256':digest(ROOT/'tools/native_eos_reference.py'),
                 'alpha':alpha,'refine':refine,'rtol':rtol,'fit':fit,'status':'running'}
        write_json(output/'execution.json',receipt)
        start=time.monotonic();e=NativeEOS(tab,alpha);target=2.644606257493485e30
        records=[];last={}
        def residual(x):
            pc=float(np.exp(x));parts=e.parts(pc,1e20)
            shift,bg=bgmodule.tov(parts,pc*CONV,r0=100.,rtol=rtol)
            mass=bg['mass_msun']*(1.3271244e20/G)
            last.update(parts=parts,shift=shift,bg=bg,pc=pc)
            records.append(dict(pc_pa=pc,mass_kg=mass))
            return mass/target-1
        x=brentq(residual,np.log(1.2e34),np.log(1.8e34),xtol=2e-14,rtol=9e-16);residual(x)
        if abs(records[-1]['mass_kg']/target-1)>1e-12:raise ValueError('Production mass target did not close')
        mesh=shoot.make_mesh(last['parts'],last['shift'],max_dr_cm=1000/refine,max_dx=.04/refine,max_dlogr=.04/refine)
        result={'background':last['bg'],'pc_pa':last['pc'],'alpha':alpha,'refine':refine,'rtol':rtol,'fit':fit,
                'target_mass_kg':target,'target_evaluations':records,'mesh':mesh.settings,'modes':{}}
        for mode,lo,hi in [('f',2000,3000),('p1',6500,8000)]:
            f=brentq(lambda f:shoot.determinant(mesh,f,fit),lo,hi,xtol=1e-10,rtol=1e-14)
            value=shoot.determinant(mesh,f,fit)
            if abs(value)>1e-8:raise ValueError('Production matching determinant did not close')
            shape=shoot.eigenfunction(mesh,f,fit)
            np.savez_compressed(output/(mode+'_shape.npz'),**{k:v for k,v in shape.items() if isinstance(v,np.ndarray)})
            result['modes'][mode]={'frequency_hz':float(f),'determinant':float(value)}
        result.update(seconds=time.monotonic()-start,status='computed_candidates',
                      limitation='Finite frequency windows; labels require physical shape comparison, not frequency order.')
        write_json(output/'result.json',result)
        receipt['status']='complete';receipt['outputs']={p.name:digest(p) for p in output.iterdir() if p.name!='execution.json'}
        write_json(output/'execution.json',receipt)
        return result
    except Exception as exc:
        write_json(output/'FAILURE.json',{'type':type(exc).__name__,'error':str(exc),'traceback':traceback.format_exc()})
        raise

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--production',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--alpha',type=float,default=0.)
    p.add_argument('--refine',type=int,default=4);p.add_argument('--rtol',type=float,default=2e-13)
    p.add_argument('--fit',type=float,default=.5);a=p.parse_args(argv)
    try:
        r=run_reference(a.production,a.output,alpha=a.alpha,refine=a.refine,rtol=a.rtol,fit=a.fit)
        print(json.dumps(r['modes'],indent=2));return 0
    except Exception as e:print(str(e),file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
