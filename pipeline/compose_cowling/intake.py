"""Read successful matched tables through the bundled compose_eos v3.3.0 API."""
from __future__ import annotations
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

HEADER = [
    "source_row","i_t","i_nb","i_yq","n_b_fm3","n_b_m3",
    "pressure_mev_fm3","pressure_pa","energy_mev_fm3","energy_j_m3",
    "mu_b_mev","mu_b_j","q1_p_over_nb_mev","q2_entropy_per_baryon",
    "q3_mu_b_over_mn_minus_1","q4_mu_q_over_mn","q5_mu_l_over_mn",
    "q6_f_over_nb_mn_minus_1","q7_e_over_nb_mn_minus_1"
]
MEVFM3_TO_JM3 = 1.602176634e32
MEV_TO_J = 1.602176634e-13

def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare_from_matched(
    matched_root: Path,
    classifier_json: Path,
    output_dir: Path,
    generator_package: Path,
) -> dict[str, Any]:
    generator_package = Path(generator_package).resolve()
    if not generator_package.is_dir():
        raise FileNotFoundError(f"compose_eos package directory not found: {generator_package}")
    sys.path.insert(0, str(generator_package))
    from compose_eos.io import read_compose

    matched_root, classifier_json, output_dir = map(Path, (matched_root, classifier_json, output_dir))
    classification = json.loads(classifier_json.read_text(encoding="utf-8"))
    results = classification.get("results")
    if not isinstance(results, list):
        raise ValueError("failure classifier does not contain a results list")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for result in results:
        dirname = result["directory"]
        directory = matched_root / dirname
        meta_file = directory / "MATCHING_METADATA.json"
        if not meta_file.is_file():
            counts["MISSING_METADATA"] += 1
            continue
        metadata = json.loads(meta_file.read_text(encoding="utf-8"))
        if result.get("status") != "PASS" or metadata.get("status") != "PASS" or metadata.get("diagnostic_only", False):
            counts["NOT_ADMITTED_SOURCE_OR_JOIN"] += 1
            continue
        eos = read_compose(str(directory))
        n, p, eps, mu = (eos.nb, eos.p, eos.eps, eos.mub)
        import numpy as np
        n, p, eps, mu = [np.asarray(x, dtype=float) for x in (n, p, eps, mu)]
        if eos.n < 3 or any(not np.all(np.isfinite(x)) for x in (n,p,eps,mu)):
            raise ValueError(f"{dirname}: nonfinite or undersized EOS")
        identity = np.abs(eps+p-mu*n) / np.maximum(np.abs(eps+p), 1e-300)
        if float(np.max(identity)) > 2e-12:
            raise ValueError(f"{dirname}: source first-law knot identity failed")
        if any(np.any(np.diff(x) <= 0) for x in (n,p,eps,mu)) or np.any(np.r_[n,p,eps,mu] <= 0):
            raise ValueError(f"{dirname}: source barotrope is not positive and strictly increasing")
        hashes = {name:_hash(directory/name) for name in ("eos.nb","eos.thermo","eos.t","eos.yq")}
        target = output_dir / f"{dirname}.csv"
        mn = float(eos.mn)
        with target.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(HEADER)
            for i in range(eos.n):
                vals = [
                    i+1, eos.iT, eos.inb_start+i, eos.iYq,
                    n[i], n[i]*1e45, p[i], p[i]*MEVFM3_TO_JM3,
                    eps[i], eps[i]*MEVFM3_TO_JM3, mu[i], mu[i]*MEV_TO_J,
                    eos.Q1[i], eos.Q2[i], eos.Q3[i], eos.Q4[i], eos.Q5[i], eos.Q6[i], eos.Q7[i]
                ]
                writer.writerow([format(float(v), ".17g") if isinstance(v, (float, int)) else v for v in vals])
        rows.append({
            "directory":dirname,"core":metadata.get("source_core"),"crust":metadata.get("source_crust"),
            "method":metadata.get("method"),"rows":int(eos.n),"table_path":str(target.resolve()),
            "mn_MeV":mn,"mp_MeV":float(eos.mp),"beta_equilibrated":bool(eos.beta_equilibrated),
            "first_law_rel_max":float(np.max(identity)),
            "cs2_sample_min":float(np.min(np.diff(p)/np.diff(eps))),
            "cs2_sample_max":float(np.max(np.diff(p)/np.diff(eps))),
            "source_hashes":hashes,"matching_metadata_sha256":_hash(meta_file),
        })
        counts["PREPARED"] += 1
    manifest = {
        "schema":"compose-cowling-prepared-eos-v1",
        "generator":"compose_eos 3.3.0",
        "matched_root":str(matched_root.resolve()),
        "classifier_sha256":_hash(classifier_json),
        "n_classified":len(results),
        "counts":dict(counts),
        "thermodynamics":"CompOSE source rows: P=n_B Q1, mu_B=m_n(1+Q3), epsilon=n_B m_n(1+Q7); independent knot identity epsilon+P=mu_B n_B",
        "entries":rows,
    }
    (output_dir/"prepared_eos_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    return manifest
