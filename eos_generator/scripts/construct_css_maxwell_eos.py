#!/usr/bin/env python3
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'compose_eos_pkg'))
from compose_eos import read_compose
from compose_eos.phase_transitions import build_maxwell_css_eos,write_maxwell_eos
p=argparse.ArgumentParser(description='Sharp Maxwell transition with a phenomenological constant-sound-speed high-density phase.')
p.add_argument('--hadronic',required=True);p.add_argument('--pressure-transition',type=float,required=True);p.add_argument('--delta-epsilon',type=float,required=True);p.add_argument('--cs2-high',type=float,default=1.0);p.add_argument('--n-max',type=float,required=True);p.add_argument('--points',type=int,default=240);p.add_argument('--out',required=True)
a=p.parse_args(); e=build_maxwell_css_eos(read_compose(a.hadronic),a.pressure_transition,a.delta_epsilon,a.cs2_high,n_max_fm3=a.n_max,n_points=a.points);write_maxwell_eos(a.out,e);print(json.dumps(e.transition,indent=2,allow_nan=False))
