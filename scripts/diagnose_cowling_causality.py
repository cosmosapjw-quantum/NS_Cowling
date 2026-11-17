"""Reproduce the Fortran ns_table_eos audit on the normalized CSV that the
Cowling solver actually consumes. Locate the rejected (superluminal /
non-monotone g) intervals across all 5 smooth methods.
"""
import sys, os, tempfile
import numpy as np
sys.path.insert(0, "python_frontend")
os.environ["COMPOSE_EOS_FORTRAN_ROOT"] = os.getcwd()
os.environ["COMPOSE_EOS_FORTRAN_EXE"] = os.path.join(os.getcwd(), "build/compose_eos_fortran/compose_eos_fortran")
from compose_eos import read_compose, match_crust_core, write_compose, normalize_compose
from pathlib import Path

CRUST = "inputs/raw_compose/SPG_M1"
CORE = "inputs/raw_compose/APR_APR"
METHODS = ["css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"]

crust = read_compose(CRUST)
core = read_compose(CORE)
print(f"crust n={crust.n}  core n={core.n}")
print(f"crust P [{min(crust.p):.3g},{max(crust.p):.3g}]  core P [{min(core.p):.3g},{max(core.p):.3g}] MeV/fm^3")
print()

for method in METHODS:
    matched, info = match_crust_core(crust, core, method=method, n_bridge=40, drop_composition=True)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td); comp = td/"compose"; comp.mkdir()
        write_compose(comp, matched)
        norm_csv = td/"normalized_eos.csv"
        normalize_compose(comp, norm_csv)
        data = np.loadtxt(norm_csv, delimiter=",", skiprows=1)
        with open(norm_csv) as f:
            cols = [c.strip() for c in f.readline().split(",")]
        jmap = {c:i for i,c in enumerate(cols)}
        P = data[:, jmap["pressure_pa"]]
        eps = data[:, jmap["energy_j_m3"]]
        nb = data[:, jmap["n_b_m3"]]
        n = len(P)
        u = np.log(P); g = eps - P
        h = np.diff(u); d = np.diff(g)/h
        cs2 = 1.0/(d/P[:-1] + 1.0)
        n_super = int(np.sum(cs2 > 1+1e-9))
        n_nonmono = int(np.sum(d <= 0))
        worst = int(np.argmin(d))
        print(f"--- {method} (n={n}) ---")
        print(f"  P>0:{np.all(P>0)} logP inc:{np.all(h>0)} g inc:{np.all(d>0)}")
        print(f"  superluminal:{n_super}/{n-1}  nonmono g:{n_nonmono}")
        print(f"  worst idx={worst} P={P[worst]:.4g}->{P[worst+1]:.4g} dg/du={d[worst]:.4g} cs2={cs2[worst]:.4f}")
        # density gaps to locate bridge
        dnb = np.diff(nb)
        bridge = int(np.argmax(dnb)) + 1
        print(f"  bridge~idx={bridge} (nb {nb[max(0,bridge-1)]:.4g}->{nb[bridge]:.4g})")
        # show intervals around worst
        lo,hi = max(0,worst-2), min(n-1, worst+3)
        for i in range(lo,hi):
            mark = "  <-- bridge" if i==bridge else ""
            dg = d[i] if i<n-1 else float('nan')
            cs = cs2[i] if i<n-1 else float('nan')
            print(f"    i={i:3d} P={P[i]:.4g} eps={eps[i]:.4g} g={g[i]:.4g} dg/du={dg:.4g} cs2={cs:.4f}{mark}")
    print()
