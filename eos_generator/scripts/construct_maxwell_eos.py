#!/usr/bin/env python3
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'compose_eos_pkg'))
from compose_eos import read_compose
from compose_eos.phase_transitions import build_maxwell_eos,write_maxwell_eos
p=argparse.ArgumentParser();p.add_argument('--low-phase',required=True);p.add_argument('--high-phase',required=True);p.add_argument('--out',required=True);p.add_argument('--low-name',default='low');p.add_argument('--high-name',default='high');a=p.parse_args()
e=build_maxwell_eos(read_compose(a.low_phase),read_compose(a.high_phase),low_name=a.low_name,high_name=a.high_name);write_maxwell_eos(a.out,e);print(json.dumps(e.metadata['transition'],indent=2))
