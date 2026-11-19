#!/usr/bin/env python3
"""Recalculate exact polytropic surface extrapolation from frozen BSk22 outer.

E's 217-point BSk22 outer table is *inherited*, not microscopic ETFSI
regenerated. This script independently re-evaluates its anchored
surface E(n),P(n),mu(n) using the declared gamma and first law.
The 10 separate lower bridge functions are recalculated in the
new SPG/PK1 three-layer TOV results, never claimed to be original BSk22.
"""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/stellar_extension_20261117/eos/fine/REPRESENTATION.json"
OUT=ROOT/"data/eos_missing_data_recalc_20261119/bsk22_outer_and_surface"
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def csvout(path,header,rows):
    with path.open("w",newline="") as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)
def main():
    OUT.mkdir(exist_ok=True,parents=True)
    data=json.loads(SOURCE.read_text())
    outer=data["outer"];surface=data["surface"]
    n=np.asarray(outer["n_fm3"],float);P=np.asarray(outer["p_MeV_fm3"],float)
    F=np.asarray(outer["F_MeV"],float)
    assert len(n)==217 and np.all(np.diff(n)>0) and np.all(np.diff(P)>0)
    epsilon=n*F
    mu=F+P/n
    assert all(np.isfinite(epsilon)) and all(epsilon>0)
    csvout(OUT/"BSK22_OUTER_INHERITED_217KNOTS.csv",
       ["n_b_fm-3","P_MeV_fm-3","epsilon_MeV_fm-3","mu_B_Euler_MeV","F_MeV"],
       zip(n,P,epsilon,mu,F))
    gamma=float(surface["gamma"]);n0=float(surface["n0"])
    p0=float(surface["p0"]);eps0=float(surface["epsilon0"])
    rest=eps0/n0-p0/(n0*(gamma-1))
    assert gamma>1 and rest>0 and math.isclose(rest,float(surface["rest_energy_MeV"]),rel_tol=1e-11)
    N=np.r_[0.,n0*np.geomspace(1e-12,1.,255)]
    Ps=p0*(N/n0)**gamma
    Es=rest*N+Ps/(gamma-1)
    Mus=np.zeros_like(N);Fsurf=np.zeros_like(N)
    Mus[1:]=rest+gamma/(gamma-1)*Ps[1:]/N[1:]
    Fsurf[1:]=Es[1:]/N[1:]
    Mus[0]=rest;Fsurf[0]=rest
    assert math.isclose(Ps[-1],p0,rel_tol=1e-13)
    assert math.isclose(Es[-1],eps0,rel_tol=1e-13)
    assert math.isclose(n[-1],outer["n_fm3"][-1])
    assert math.isclose(n[0],n0,rel_tol=2e-8)
    assert abs(epsilon[0]-eps0)/eps0<5e-8
    csvout(OUT/"ANCHORED_SURFACE_GAMMA_RECALCULATED_256KNOTS.csv",
        ["n_b_fm-3","P_MeV_fm-3","epsilon_MeV_fm-3","mu_B_Euler_MeV","F_MeV"],
        zip(N,Ps,Es,Mus,Fsurf))
    receipt=dict(schema="ns_cowling.bsk22_frozen_outer_surface_independent_evaluation.v1",
        source_path=str(SOURCE.relative_to(ROOT)),source_sha256=sha(SOURCE),
        source_outer_rows=217,independent_surface_rows=len(N),
        outer_source_role="E PREVIOUSLY RECONSTRUCTED BSK22; ORIGINAL 217 KNOTS UNCHANGED; NOT REGENERATED ETFSI",
        surface_status="FRESH_ANALYTIC_FIRSTLAW_EVALUATION_FROM_SOURCE_ANCHOR",
        surface_gamma=gamma,surface_rest_mass_MeV=rest,
        exact_anchor_pressure_relative_error=abs(Ps[-1]/p0-1),
        exact_anchor_energy_relative_error=abs(Es[-1]/eps0-1),
        source_neutron_drip_n_fm3=n[-1],surface_bottom_n_fm3=n[0],
        surface_analytic_identity="P=n^2 d(epsilon/n)/dn, epsilon=m*n+P/(gamma-1), mu=(epsilon+P)/n",
        zero_surface=dict(n=0.,P=0.,epsilon=0.,mu=rest),
        surface_and_crust_microscopic_unification=False)
    (OUT/"RECALC_REPORT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("BSK22_OUTER",len(n),"SURFACE",len(N),"ANCHOR relative",
          receipt["exact_anchor_energy_relative_error"],flush=True)
if __name__=="__main__":main()
