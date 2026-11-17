"""Generate synthetic core EOS CompOSE tables for the M7 batch demo.

The 390-matrix core models (APR_APR, BBB_BHF2, GM_GM1, RG_SLY4, RG_KDE0v1)
are not present in this sandbox (they would come from the M0 replacement
ZIP). To exercise the M7 batch backend (M5 Maxwell + M5A CSS) on real PASS
paths, we generate synthetic monotone core barotropes from polytropes and
write them as CompOSE tables to inputs/core/<model>/. These are clearly
labelled as SYNTHETIC and are NOT the production core models; they exist
only to validate that the batch orchestrator dispatches and records PASS
receipts when sources are available.
"""
from __future__ import annotations
import sys, os, json
import numpy as np
sys.path.insert(0, '/home/z/my-project/repos/NS_Cowling_eos390/python_frontend')
from compose_eos.eos_data import EOSData
from compose_eos.api import write_compose

CORE_SYNTH = {
    # Each synthetic core: (K, gamma, mn) for a polytrope P = K n^gamma
    # tuned to be a stiff high-density barotrope (representative of nuclear
    # matter above saturation). The exact numbers are placeholders; production
    # uses the real APR/BBB/GM/RG models.
    'APR_APR':    (1500.0, 2.5, 939.0),
    'BBB_BHF2':   (1300.0, 2.4, 939.0),
    'GM_GM1':     (1800.0, 2.6, 939.0),
    'RG_SLY4':    (1400.0, 2.45, 939.0),
    'RG_KDE0v1':  (1600.0, 2.55, 939.0),
}

def make_polytrope_eos(K, gamma, mn=939.0, n_min=0.08, n_max=0.80, n_pts=80):
    nB = np.linspace(n_min, n_max, n_pts)
    P = K * nB ** gamma
    eps = mn * nB + P / (gamma - 1.0)
    mu_B = (eps + P) / nB
    F = eps / nB
    return EOSData(
        nb=nB, mn=mn, mp=mn * 0.9938,
        Q1=P / nB, Q2=np.zeros_like(nB),
        Q3=mu_B / mn - 1.0, Q4=np.zeros_like(nB), Q5=np.zeros_like(nB),
        Q6=F / mn - 1.0, Q7=F / mn - 1.0,
        q_extra=np.zeros_like(nB), Iphase=np.zeros_like(nB, dtype=int),
        Yn=np.zeros_like(nB), Yp=np.zeros_like(nB), Ye=np.zeros_like(nB),
        Ymu=np.zeros_like(nB), Aav=np.zeros_like(nB), Zav=np.zeros_like(nB),
        Yh=np.zeros_like(nB),
        Il=1, iT=0, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_available=False, beta_equilibrated=True,
    )

def main():
    out_root = '/home/z/my-project/download/m7_batch/inputs_core_synth'
    os.makedirs(out_root, exist_ok=True)
    provenance = {}
    for model, (K, gamma, mn) in CORE_SYNTH.items():
        d = os.path.join(out_root, model)
        os.makedirs(d, exist_ok=True)
        eos = make_polytrope_eos(K, gamma, mn)
        write_compose(d, eos)
        provenance[model] = {
            'kind': 'synthetic_polytrope',
            'K': K, 'gamma': gamma, 'mn': mn,
            'n_range': [0.08, 0.80], 'n_pts': 80,
            'note': 'SYNTHETIC core placeholder for M7 batch demo; NOT the production model',
        }
        print(f'  wrote synthetic core {model}: K={K} gamma={gamma} -> {d}')
    with open(os.path.join(out_root, '_synth_provenance.json'), 'w') as f:
        json.dump(provenance, f, indent=2)
    print(f'\nprovenance written: {out_root}/_synth_provenance.json')

if __name__ == '__main__':
    main()
