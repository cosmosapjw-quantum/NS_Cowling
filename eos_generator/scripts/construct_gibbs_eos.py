#!/usr/bin/env python3
"""Input NPZ files must contain mu_b, mu_q and pressure arrays (MeV units)."""
import argparse,json,sys,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'compose_eos_pkg'))
from compose_eos.phase_transitions import ChargeResolvedPhase,gibbs_mixed_phase,write_gibbs_eos
p=argparse.ArgumentParser();p.add_argument('--hadronic-npz',required=True);p.add_argument('--high-density-npz',required=True);p.add_argument('--mu-b-min',type=float,required=True);p.add_argument('--mu-b-max',type=float,required=True);p.add_argument('--points',type=int,default=120);p.add_argument('--out',required=True);a=p.parse_args()
h=ChargeResolvedPhase.from_npz(a.hadronic_npz,'hadronic');q=ChargeResolvedPhase.from_npz(a.high_density_npz,'high-density');r=gibbs_mixed_phase(h,q,np.linspace(a.mu_b_min,a.mu_b_max,a.points));write_gibbs_eos(a.out,r);print(json.dumps({'points':r.eos.n,'chi_min':float(np.min(r.metadata['mixed_fraction_chi'])),'chi_max':float(np.max(r.metadata['mixed_fraction_chi']))},indent=2))
