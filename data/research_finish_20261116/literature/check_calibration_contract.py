#!/usr/bin/env python3
"""Reproduce source-cohort and density-contract checks; no EOS fit or solver import."""
from pathlib import Path
import csv
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    cuter = HERE / "cuter_reference"
    source_hash = digest(cuter / "nuclear_surface_en.c")
    data_hash = digest(cuter / "AME2016.data")
    assert source_hash == "510a03a4643643716e6be67d3b0feac2b1d4c094494f280ab16ec4e857ff014e"
    assert data_hash == "04d55e1d291170f0a7b991cba07e6d4b8d1e764022768cffd801ad4b345dcfd9"
    reference = {}
    for row in (cuter / "AME2016.data").read_text().splitlines():
        z, a, negative_binding = row.split()
        key = (int(a), int(z))
        assert key not in reference
        reference[key] = -float(negative_binding)

    path = ROOT / "data/research_generator_20261116/candidate_surface/candidate_surface_residuals.csv"
    all_rows = list(csv.DictReader(path.open()))
    label = all_rows[0]["model"]
    candidate = {
        (int(row["A"]), int(row["Z"])): float(row["B_AME_MeV"]) / int(row["A"])
        for row in all_rows if row["model"] == label
    }
    assert len(candidate) == len(reference) == 2408
    assert set(candidate) == set(reference)
    cohort_error = max(abs(candidate[key] - reference[key]) for key in reference)
    assert cohort_error < 1e-12

    # This analytic model isolates the meaning of saturation; it is not any SPG EOS.
    nsat, esat, jsym, lsym, ksat, ksym = .16, -16., 32., 60., 240., -100.
    checks = []
    for delta in [0., .2, .4, .6]:
        i2 = delta * delta
        curvature = ksat + ksym * i2
        n0 = nsat * (1. - 3. * lsym * i2 / curvature)
        x0 = (n0 - nsat) / (3. * nsat)
        e = esat + jsym * i2 + lsym * i2 * x0 + .5 * curvature * x0 * x0
        de_dn = (lsym * i2 + curvature * x0) / (3. * nsat)
        pressure = n0 * n0 * de_dn
        # epsilon = n e, hence derivative epsilon = e + n derivative e.
        deps_dn = e + n0 * de_dn
        assert abs(pressure) < 1e-13
        assert abs(deps_dn) > 1.
        checks.append(dict(delta=delta, n0_fm3=n0, P_MeV_fm3=pressure,
                           d_binding_energy_density_dn_MeV=deps_dn))

    # At fixed A,Z and tension depending only on y, Es~n^-2/3,
    # Ec~n^1/3 and Ecurv~n^-1/3. Finite-difference the full mass.
    finite_checks = []
    for a, z in [(16., 8.), (90., 40.), (208., 82.)]:
        y = z / a
        delta = 1. - 2. * y
        curvature = ksat + ksym * delta * delta
        for n in [.08, .12, .16]:
            def terms(density):
                radius = (3. * a / (4. * math.pi * density)) ** (1. / 3.)
                x = (density - nsat) / (3. * nsat)
                e = esat + jsym * delta * delta + lsym * delta * delta * x + .5 * curvature * x*x
                es = 4. * math.pi * radius**2 * 1.05
                ec = 3. * 1.43996448 * z*z / (5. * radius)
                ek = 8. * math.pi * radius * .12
                return a*e, es, ec, ek
            x = (n - nsat) / (3. * nsat)
            de_dn = (lsym * delta*delta + curvature*x) / (3. * nsat)
            bulk, es, ec, ek = terms(n)
            derivative = a * de_dn + (-2.*es + ec - ek)/(3.*n)
            h = 1e-6 * n
            finite = (sum(terms(n+h)) - sum(terms(n-h))) / (2.*h)
            scaled_error = abs(finite-derivative) / max(1.,abs(derivative))
            assert scaled_error < 1e-7, (a,z,n,scaled_error)
            # Same relation rearranged as a pressure residual.
            residual_pressure = n*n*de_dn - n*(2.*es-ec+ek)/(3.*a)
            assert abs(residual_pressure - derivative*n*n/a) < 1e-12
            finite_checks.append(dict(A=a,Z=z,n_fm3=n,scaled_derivative_error=scaled_error))
    result = {
        "status":"PASS_SCOPED",
        "claim":"Cohort identity and finite-nucleus density-stationarity algebra only; not CUTER execution, an SPG fit, or EOS calibration.",
        "source_sha256":source_hash,"data_sha256":data_hash,
        "candidate_AMEmass_rows_sha256":digest(path),
        "cohort_count":len(reference),"cohort_max_binding_per_baryon_error_MeV":cohort_error,
        "bulk_quadratic_examples":checks,
        "finite_nucleus_derivative_checks":finite_checks,
        "max_finite_difference_scaled_error":max(x["scaled_derivative_error"] for x in finite_checks),
    }
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
