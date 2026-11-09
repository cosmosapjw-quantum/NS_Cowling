#!/usr/bin/env python3
"""Strict checks for a pilot or full mass-dependent Cowling artifact."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path

def read(path):
    with path.open(newline="",encoding="utf-8") as f: return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser();p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--expected-models",type=int,required=True);p.add_argument("--required-mass",type=float,default=1.33)
    p.add_argument("--required-low-mass",type=float,default=1.0)
    a=p.parse_args();root=a.output_dir;m=json.loads((root/"manifest.json").read_text())
    models=read(root/"model_branch_summary.csv");spectra=read(root/"spectrum_mass_grid.csv")
    paired=read(root/"paired_mass_sensitivity.csv")
    assert m["selected_model_count"]==a.expected_models and m["completed_model_count"]==a.expected_models
    assert m["failed_model_count"]==0 and not m["model_failures"] and not m["baseline_reproduction_failures"]
    assert len(models)==a.expected_models and len(spectra)==a.expected_models*len(m["mass_grid_msun"])
    assert all(r["baseline_reproduction_status"]=="pass" for r in models)
    for target in (a.required_low_mass,a.required_mass):
        required=[r for r in paired if abs(float(r["mass_msun"])-target)<1e-9]
        assert required and all(r["pair_complete"]=="true" for r in required), \\
            f"required mass {target} lacks a complete six-method pair"
    partial=[r for r in paired if r["pair_complete"]!="true" and int(r["n_valid_methods"])>0]
    assert not partial, f"partially valid six-method cells must be investigated, found {len(partial)}"
    valid=[r for r in spectra if r["valid"]=="true"]
    assert all(abs(float(r["mass_relative_residual"]))<=5e-10 for r in valid)
    coverage={}
    for target in m["mass_grid_msun"]:
        cells=[r for r in paired if abs(float(r["mass_msun"])-float(target))<1e-9]
        coverage[str(target)]={"complete_pairs":sum(r["pair_complete"]=="true" for r in cells),
            "zero_valid_pairs":sum(int(r["n_valid_methods"])==0 for r in cells)}
    print(json.dumps({"status":"PASS","models":len(models),"spectrum_rows":len(spectra),
        "valid_spectra":len(valid),"required_masses":[a.required_low_mass,a.required_mass],
        "mass_coverage":coverage},indent=2))

if __name__=="__main__": main()
