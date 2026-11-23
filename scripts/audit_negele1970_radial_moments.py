#!/usr/bin/env python3
"""Evaluate Negele (1970) radial moments and NV72 finite-range B,C,D cross-checks."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn

from audit_nv1972_skyrme_tables import INPUT, coefficients

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/eos_reconstructed_matching_20261124/source_audit"
TABLE = DATA / "negele1970_table_iv.csv"
OUTPUT = DATA / "negele1970_radial_moments.json"
# Negele (1970) Table V(b), columns V0,V1. Eqs. (3.16), (3.24).
FACTORS = {
    "like": {"short": (.671, .840), "direct": (.650, .409), "exchange": (.329, -3.328)},
    "unlike": {"short": (1.329, 1.160), "direct": (1.350, 1.591), "exchange": (1.671, 5.328)}
}
ANCHORS = [(0.080, 0.080), (0.093, 0.064), (0.064, 0.093),
           (0.040, 0.040), (0.030, 0.010), (0.005, 0.005)]


def local_cubic(nodes, values, r):
    """Four-point interpolation on each covered interval; zero outside the grid."""
    start = np.clip(np.searchsorted(nodes, r) - 2, 0, len(nodes) - 4)
    x = nodes[start[:, None] + np.arange(4)]
    y = values[start[:, None] + np.arange(4)]
    value = np.zeros_like(r)
    for i in range(4):
        others = np.arange(4) != i
        li = np.prod((r[:, None]-x[:, others])/(x[:, i, None]-x[:, others]), axis=1)
        value += y[:, i]*li
    return np.where((r >= nodes[0]) & (r <= nodes[-1]), value, 0)


class RadialKernels:
    def __init__(self, quadrature_order, rmax):
        with TABLE.open() as input_file:
            rows = list(csv.DictReader(input_file))
        edges = np.array([0.0] + [float(row["r_fm"]) for row in rows if float(row["r_fm"]) <= rmax])
        x, w = leggauss(quadrature_order)
        self.r = (edges[:-1, None] + np.diff(edges)[:, None]*(x+1)/2).ravel()
        self.w = (np.diff(edges)[:, None]*w/2).ravel()
        self.v = {}
        for channel in FACTORS["like"]:
            for k in (1.0, 1.4):
                col = f"{channel}_{k}"
                pairs = [(float(row["r_fm"]), float(row[col])) for row in rows if row[col]]
                if channel == "short":
                    pairs.insert(0, (0., 0.))
                nodes, values = np.array(pairs).T
                self.v[channel, k] = local_cubic(nodes, values, self.r)

    def potential(self, channel, pair, k):
        power = 3 if channel == "short" else 1
        v1 = (self.v[channel, 1.4]-self.v[channel, 1.0])/(1.4**power-1)
        v0 = self.v[channel, 1.0]-v1
        a, b = FACTORS[pair][channel]
        return a*v0 + b*v1*k**power

    def direct(self, rho, pair):
        k = kf(rho)
        v = self.potential("short", pair, k) + self.potential("direct", pair, k)
        return np.pi*np.dot(self.w*self.r**4, g(k*self.r)*v)

    def exchange(self, a, b, pair):
        x = kf(a)*self.r
        slater = 3*spherical_jn(1, x)/x
        v = self.potential("exchange", pair, kf((a+b)/2))
        return 2*np.pi*np.dot(self.w*self.r**4, slater*g(kf(b)*self.r)*v)

    def B(self, n, p):
        return -2*(n*self.exchange(n, n, "like") + p*self.exchange(p, n, "unlike"))

    def F(self, n, p):
        return n*self.direct(n, "like") + p*self.direct((n+p)/2, "unlike") - self.B(n, p)/4

    def bcd(self, n, p, h):
        perp1 = (self.F(n, p+h)-self.F(n, p-h))/(2*h)
        perp2 = (self.F(p, n+h)-self.F(p, n-h))/(2*h)
        para = (self.F(n+h, p)-self.F(n-h, p))/(2*h)
        c = -self.direct(n, "like")-para
        d = -2*self.direct((n+p)/2, "unlike")-perp1-perp2
        return np.array([self.B(n, p), c, d])


def kf(rho):
    return (3*np.pi**2*rho)**(1/3)


def g(x):
    return 17.5*spherical_jn(3, x)/x**3


def table_values(tables, n, p):
    return [float(np.polynomial.polynomial.polyval2d(n, p, np.asarray(t, dtype=float)))
            for t in tables]


def properties(bcd):
    return dict(zip(("B", "C", "D"), map(float, bcd)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()
    kernels = RadialKernels(12, 7.0)
    refined = RadialKernels(24, 7.0)
    clipped = RadialKernels(12, 6.6)
    references = {(1.0, "short"): 25.091, (1.0, "direct"): -39.730,
                  (1.4, "short"): 26.029, (1.4, "direct"): -33.089}
    checks = []
    for (k, channel), target in references.items():
        value = float(np.dot(kernels.w*kernels.r**2, kernels.v[channel, k]))
        checks.append({"kF": k, "channel": channel, "reference": target, "reconstructed": value,
                       "relative_difference": value/target-1})
    if max(abs(row["relative_difference"]) for row in checks) > .002:
        raise RuntimeError("1970 Table IV transcription or interpolation regression against Table V(a)")

    # The normalization of the exchange weight cancels at a fixed Fermi sphere.
    # Its r dependence is the square of the normalized Fermi-sphere Fourier transform.
    exchange_ratios = []
    for sphere, moment_10, moment_14 in [(1.0, -2.629, -2.560), (1.4, -5.415, -5.311)]:
        x = sphere*kernels.r
        weight = kernels.w*kernels.r**2 * (3*spherical_jn(1, x)/x)**2
        value = float(np.dot(weight, kernels.v["exchange", 1.4])
                      / np.dot(weight, kernels.v["exchange", 1.0]))
        reference = moment_14/moment_10
        exchange_ratios.append({"weight_kF": sphere, "Va_average_moment_k1_4_over_k1_0": reference,
                                "quadrature_ratio": value, "relative_difference": value/reference-1})
    if max(abs(row["relative_difference"]) for row in exchange_ratios) > .002:
        raise RuntimeError("1970 Table IV exchange interaction regression against Table V(a)")

    source = json.loads(INPUT.read_text())
    printed = coefficients(source)
    d03x100 = coefficients(source, d03_factor=100)
    rows = []
    for n, p in ANCHORS:
        raw = kernels.bcd(n, p, 1e-6)
        rows.append({
            "n": n, "p": p,
            "from_1970_average_radial_factors": properties(raw),
            "q24_h_over_two_minus_q12": properties(refined.bcd(n, p, 5e-7)-raw),
            "rmax_6_6_minus_rmax_7_fm": properties(clipped.bcd(n, p, 1e-6)-raw),
            "NV1972_literal": properties(table_values(printed, n, p)),
            "NV1972_D03x100": properties(table_values(d03x100, n, p))
        })

    result = {
        "source": {
            "negele1970_doi": "10.1103/PhysRevC.1.1260",
            "negele1970_original_pdf_sha256": "a94340c7704afd4baa119e2ab19cf91ce4b9fa7fb995a8b96759fe7b2e26132b",
            "locations": ["p.1275 Table IV", "pp.1276–77 Eqs. (3.13–3.24) and Table V", "Appendix A, four-point interpolation"],
            "NV1972_doi": "10.1103/PhysRevC.5.1472",
            "NV1972_locations": ["pp.1476–78 Eqs. (3.4),(3.7),(3.8),(3.11)", "pp.1486,1488 energy adjustment and fitting grid"],
            "NV1973_doi": "10.1016/0375-9474(73)90349-7",
            "NV1973_locations": ["pp.305–306 Table 1 supplies the new bulk polynomial, B/C/D imported from NV72"],
            "transcription": "Table IV averaged potentials in MeV; empty CSV fields indicate no published radial sample in that channel. Table V(b) multiplying factors are in the script.",
            "table_iv_csv_sha256": hashlib.sha256(TABLE.read_bytes()).hexdigest()
        },
        "method": [
            "Equation (3.16) interpolates each averaged radial interaction as v0+v1*kF^a with a=3 short and a=1 otherwise. Equation (3.24) applies the Table V(b) like/unlike factors separately to v0 and v1.",
            "Radial interpolation uses the local four-point Lagrange rule. The short-range value at the unsampled origin is set to zero. Direct and exchange vanish below their tabulated thresholds. All integrations stop at 7 fm.",
            "Each tabulated interval has 12-point Gauss-Legendre quadrature. NV72 uses g(x)=(35/2)j3(x)/x^3 and the Slater function 3j1(x)/x. From the spin averages, VD=(1/4) integral d^3r r² g v_dir, VE=(1/2) integral d^3r r² Slater g v_ex. F uses Eq. (3.8) and C,D use centered differences of F with h=10^-6 fm^-3.",
            "The q24 check doubles quadrature and halves the finite-difference step. The rmax check omits the final radial interval to quantify finite radial coverage."
        ],
        "units": {"n,p": "fm^-3", "kF": "fm^-1", "r": "fm", "source_moment": "MeV fm^3",
                  "B": "MeV fm^2", "C,D": "MeV fm^5"},
        "table_Va_av_source_moments": checks,
        "table_Va_exchange_moment_ratios": exchange_ratios,
        "cross_checks": rows,
        "limitations": [
            "The radial calculation uses the 1970 factorized, spin-averaged forms under Eq. (3.24); the separate unscaled like/unlike radial functions and historical 16 by 16 least-squares fit weights are unavailable here.",
            "The separate starting-energy correction in the 1970 treatment is outside this three-channel check. The 1972 finite-nucleus fit and revised bulk correction are not reconstructed by this radial diagnostic.",
            "Residuals against Table I are significant. At symmetric saturation the reconstructed D is closer to D03x100 than to literal D. This supports a scale check; this quadrature does not establish the exact coefficients or a Wigner-Seitz EOS."
        ]
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
