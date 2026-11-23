#!/usr/bin/env python3
"""Evaluate Negele and Vautherin (1972) Table II sets 1/2 independently of the EOS solver."""
import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/eos_nv_dme_20261123/source_contract/nv1972_table1_coefficients.json"
ALTERNATIVES = ROOT / "data/eos_nv_dme_20261123/review/EXACT_DEGREE3_AMBIGUITY.json"
DEFAULT_OUTPUT = ROOT / "data/eos_reconstructed_matching_20261124/source_audit/nv1972_table_ii.json"


def coefficients(source, b04=None, d03_factor=1):
    b = [list(map(Decimal, row)) for row in source["B_rows"]]
    c = [list(map(Decimal, row)) for row in source["C_rows"]]
    upper = source["D_upper_triangle_rows_starting_at_diagonal"]
    d = [[Decimal(0) for _ in upper] for _ in upper]
    for i, row in enumerate(upper):
        for j, entry in enumerate(row, i):
            d[i][j] = d[j][i] = Decimal(entry)
    d[0][3] = d[3][0] = d[0][3] * d03_factor
    if b04 is not None:
        b[0][4] = b04
    return b, c, d


def poly(rows, n, p):
    return sum(value * n**i * p**j
               for i, row in enumerate(rows) for j, value in enumerate(row))


def skyrme_pair(t1, t2):
    return {"t1": float(t1), "t2": float(t2)}


def evaluate(tables):
    b, c, d = tables
    n = p = Decimal(".080")
    rho = n + p
    b1, c1, d1 = (poly(table, n, p) for table in (b, c, d))
    t1 = (8*c1 + 4*d1 + 4*b1/rho)/3
    t2 = (12*b1/rho - 8*c1 - 4*d1)/5
    target1, target2 = Decimal("381.1"), Decimal("15.1")
    set1 = {
        "B": float(b1), "C": float(c1), "D": float(d1),
        "eq_6_12": skyrme_pair(t1, t2),
        "minus_printed_table_ii": skyrme_pair(t1-target1, t2-target2),
        "both_inside_half_last_digit": max(abs(t1-target1), abs(t2-target2)) <= Decimal(".05")
    }

    n, p = Decimal(".093"), Decimal(".064")
    rho, delta = n+p, n-p
    bn, bp = poly(b, n, p), poly(b, p, n)
    alpha, beta = (bn+bp)/2, bn-bp
    target1, target2 = Decimal("391.9"), Decimal("8.3")
    target_alpha = rho*(3*target1+5*target2)/16
    target_beta = delta*(target2-target1)/8
    set2 = {
        "Bn": float(bn), "Bp": float(bp), "alpha": float(alpha), "beta": float(beta),
        "alpha_from_table_ii_and_preceding_identity": float(target_alpha),
        "beta_from_table_ii_and_preceding_identity": float(target_beta),
        "relative_mean_B_difference": float(alpha/target_alpha-1),
        "eq_6_14_as_printed": skyrme_pair(2*alpha/rho-5*beta/delta, 2*alpha/rho+beta/delta),
        "inverted_preceding_tau_identity": skyrme_pair(2*alpha/rho-5*beta/delta, 2*alpha/rho+3*beta/delta)
    }
    return {"set1": set1, "set2": set2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source = json.loads(INPUT.read_text())
    witnesses = json.loads(ALTERNATIVES.read_text())["witnesses"]
    alt = next(w for w in witnesses if w["edited_entry"] == "B04")
    num, den = map(Decimal, alt["values"]["B04"]["rational"].split("/"))

    with localcontext() as ctx:
        ctx.prec = 60
        cases = {
            "printed": evaluate(coefficients(source)),
            "D03x100": evaluate(coefficients(source, d03_factor=100)),
            "B04_alternative": evaluate(coefficients(source, b04=num/den))
        }
    passed = [name for name, values in cases.items() if values["set1"]["both_inside_half_last_digit"]]
    if passed != ["D03x100"]:
        raise RuntimeError(f"unexpected symmetric-source regression: {passed}")
    result = {
        "doi": "10.1103/PhysRevC.5.1472",
        "source_pdf_sha256": source["source"]["pdf_sha256"],
        "coefficient_json_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "formula_locations": ["p. 1487 Table I", "p. 1488 Table II and Eq. (6.12)", "p. 1489 Eq. (6.14) and preceding kinetic-density identity"],
        "units": {"n,p,rho,delta": "fm^-3", "B,alpha,beta": "MeV fm^2", "C,D,t1,t2": "MeV fm^5"},
        "set1_reported": {"n": 0.080, "p": 0.080, "t1": 381.1, "t2": 15.1},
        "set2_reported": {"n": 0.093, "p": 0.064, "t1": 391.9, "t2": 8.3},
        "definitions": {
            "rho": "n+p", "delta": "n-p", "alpha": "(Bn+Bp)/2", "beta": "Bn-Bp",
            "set1_eq_6_12": ["(8C+4D+4B/rho)/3", "(12B/rho-8C-4D)/5"],
            "set2_eq_6_14_print": ["2alpha/rho-5beta/delta", "2alpha/rho+beta/delta"],
            "set2_inverted_preceding_identity": ["2alpha/rho-5beta/delta", "2alpha/rho+3beta/delta"]
        },
        "branches": {"printed": "Literal Table I", "D03x100": "Multiply symmetric D03 and D30 by 100",
                     "B04_alternative": "Replace B04 only with 28658467/200 from the prior degree-three witness"},
        "results": cases,
        "interpretation": [
            "At the quoted symmetric density only D03x100 meets both printed Table II rounding intervals among these three branches.",
            "At the quoted asymmetric densities the printed B mean reproduces the independent kinetic-density combination to about 0.033 percent; its difference does not reproduce the tabulated pair. Changing B04 overshoots the mean strongly.",
            "The printed asymmetric inversion uses +beta/delta in t2. Its preceding kinetic-density identity instead gives +3beta/delta. Neither inversion reproduces both Table II set 2 entries using the printed B.",
            "The symmetric test is evidence for this coefficient among the tested alternatives. A full original polynomial fit, equivalence of Wigner-Seitz energy minimizers and a published author correction are separate claims."
        ]
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
