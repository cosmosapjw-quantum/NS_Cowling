#!/usr/bin/env python3
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'compose_eos_pkg'))
from compose_eos.phase_transitions import read_maxwell_eos,solve_tov_with_phase_transition
p=argparse.ArgumentParser();p.add_argument('--eos',required=True);p.add_argument('--pc',type=float,required=True);p.add_argument('--out',type=Path);a=p.parse_args();r=solve_tov_with_phase_transition(read_maxwell_eos(a.eos),a.pc);s=json.dumps(r,indent=2,allow_nan=False);print(s)
if a.out:a.out.write_text(s+'\n')
