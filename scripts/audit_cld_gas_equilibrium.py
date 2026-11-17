"""Audit committed CLD crust tables: gas mu_n vs tabulated mu_B.

At the variational CLD optimum the cell baryon chemical potential equals the
neutron chemical potential of the uniform gas.  For every emitted table in
CRUST_ROOT the cells are re-minimized (deterministic) at the tabulated
densities and the gas mu_n is compared with mu_B = (Q3 + 1) m_n read from
eos.thermo.  The cluster bulk mu_n is reported separately; it omits the
surface and Coulomb derivative terms.

Usage: python scripts/audit_cld_gas_equilibrium.py CRUST_ROOT OUT_JSON
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.cld_cell import minimize_cld_cell  # noqa: E402
from compose_eos.rmf_solver import load_models  # noqa: E402


def main(crust_root: Path, out_json: Path) -> int:
    models = load_models(ROOT / "data" / "eos390_rmf_parameter_sets.v1.json")
    out = {"schema": "ns-cowling.cld-gas-equilibrium-audit.v1",
           "git_head": subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                      capture_output=True, text=True).stdout.strip(),
           "crust_root": str(crust_root.relative_to(ROOT)) if crust_root.is_relative_to(ROOT) else str(crust_root),
           "models": {}}
    for d in sorted(p for p in crust_root.iterdir() if (p / "eos.thermo").is_file()):
        receipt = json.loads((d / "reconstruction_receipt.json").read_text())
        prov = receipt["provenance"]
        lines = (d / "eos.thermo").read_text().splitlines()
        mn = float(lines[0].split()[0])
        rows = [l.split() for l in lines[1:] if l.strip()]
        nb = [float(x) for x in (d / "eos.nb").read_text().split()[2:]]
        mu_tab = [(float(r[5]) + 1.0) * mn for r in rows]
        model = models[d.name]
        recs = []
        for n, mu in zip(nb, mu_tab):
            c = minimize_cld_cell(model, n, prov["sigma_0"], prov["b"])
            recs.append({"n_B": n, "mu_B_table": mu, "mu_n_gas": c.mu_n_gas,
                         "mu_n_cluster_bulk": c.mu_n,
                         "gas_minus_table_MeV": c.mu_n_gas - mu,
                         "cluster_bulk_minus_table_MeV": c.mu_n - mu})
        gas = [abs(r["gas_minus_table_MeV"]) for r in recs]
        out["models"][d.name] = {
            "eos_thermo_sha256": hashlib.sha256((d / "eos.thermo").read_bytes()).hexdigest(),
            "max_abs_gas_minus_table_MeV": max(gas),
            "max_rel_gas_minus_table": max(g / r["mu_B_table"] for g, r in zip(gas, recs)),
            "max_abs_cluster_bulk_minus_table_MeV": max(abs(r["cluster_bulk_minus_table_MeV"]) for r in recs),
            "rows": recs}
        print(f"{d.name:12s} max|gas-table|={max(gas):.4f} MeV  "
              f"max|cluster_bulk-table|={out['models'][d.name]['max_abs_cluster_bulk_minus_table_MeV']:.3f} MeV")
    out_json.write_text(json.dumps(out, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]).resolve(), Path(sys.argv[2])))
