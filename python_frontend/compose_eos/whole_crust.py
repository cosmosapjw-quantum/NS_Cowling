"""Explicit whole-crust barotropes for the reference-assisted A baseline.

Units are n in fm^-3, P/epsilon in MeV fm^-3, and chemical potentials in MeV.
The surface and interpolation prescriptions here are named approximations,
not a microscopic nuclear calculation. Raw provider files remain immutable.

The public outer table's rounded energy/pressure knots need not admit a
convex energy-preserving interpolant. ``reconstruct_outer_pressure`` therefore
defines a separate pressure-primary representation: keep every n/P knot,
anchor epsilon/n at the first source row, and integrate the first law exactly
on each positive power-law pressure interval. Isotope gaps must be labelled
explicitly by the caller; they are regularizations, not recovered Maxwell
coexistence boundaries. Their interior composition is unspecified.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Sequence

import numpy as np
from numpy.polynomial import Polynomial

from .crust_assembly import ROW_FIELDS, _cold_slice, rebase_header
from .eos_data import EOSData


def _scalar_or_array(value, original):
    value = np.asarray(value)
    return float(value) if np.asarray(original).ndim == 0 else value


def _positive_finite(*values):
    if any(np.any(~np.isfinite(x)) or np.any(np.asarray(x) <= 0) for x in values):
        raise ValueError("finite positive thermodynamic state required")


def _thermo_table(template: EOSData, n, p, epsilon, mu, *, phase_id: int) -> EOSData:
    """Build a thermo-only derived piece; composition masks are all missing."""
    n, p, epsilon, mu = (np.asarray(v, dtype=float).reshape(-1)
                          for v in (n, p, epsilon, mu))
    if len({a.size for a in (n, p, epsilon, mu)}) != 1:
        raise ValueError("thermodynamic columns must have equal lengths")
    _positive_finite(n, p, epsilon, mu)
    if np.any(np.diff(n) <= 0):
        raise ValueError("derived density grid must increase")
    count = n.size
    fields = {name: np.zeros(count, dtype=np.asarray(getattr(template, name)).dtype)
              for name in ROW_FIELDS}
    fields.update(nb=n.copy(), Q1=p/n, Q3=mu/template.mn-1.,
                  Q6=epsilon/(n*template.mn)-1., Q7=epsilon/(n*template.mn)-1.,
                  Iphase=np.full(count, int(phase_id), dtype=np.int64),
                  # A phase-only eos.compo row is PARTIAL, not an observed
                  # zero composition. This is the native reader's canonical
                  # round-trip representation; all species masks stay false.
                  row_status=np.ones(count, dtype=np.int64))
    return replace(template, source_dir=None, q_extra_count=0, nadd_max=0,
                   composition_file_present=True, composition_available=False,
                   # Cold beta-barotrope is a declared effective assumption.
                   beta_equilibrated=True, **fields)


def normalize_public_outer_composition(data: EOSData):
    """Adapt the verified SPG-provider pure-nuclear outer group convention.

    The caller must first verify raw NP=4 with species {10,11,0,1}, one
    nuclear group, raw group value=1, and source identity. In these provider
    outer rows that group value is a baryon fraction, although a standard
    CompOSE number-fraction decoder returns Yh=A*1. Pure nuclei, no free
    nucleons/muons, and Ye=Z/A determine the internal Yh=1 uniquely.

    This deliberately rejects mixed matter and absent composition evidence.
    The 2e-5 absolute charge-closure tolerance is the existing printed-table
    closure convention; the actual source's worst |Ye-Z/A| is 5e-6. Ye and
    every thermodynamic column are preserved without rounding or refitting.
    """
    _cold_slice(data)
    required_masks = ("Yn_present", "Yp_present", "Ye_present", "Ymu_present",
                      "Abar_present", "Zbar_present", "Yh_present")
    if not data.composition_file_present or any(not np.all(getattr(data, name)) for name in required_masks):
        raise ValueError("outer composition adapter requires all seven observed composition fields")
    if np.any(data.Iphase != 0):
        raise ValueError("outer composition adapter is restricted to verified provider phase 0")
    for name in ("Yn", "Yp", "Ymu"):
        if np.any(getattr(data, name) != 0):
            raise ValueError("outer composition adapter requires pure nuclei without free nucleons or muons")
    A, Z = data.Aav, data.Zav
    if (np.any(~np.isfinite(A)) or np.any(~np.isfinite(Z)) or np.any(A <= 0)
            or np.any(Z <= 0) or np.any(Z > A)
            or np.any(abs(A-np.rint(A)) > 1e-8) or np.any(abs(Z-np.rint(Z)) > 1e-8)):
        raise ValueError("outer composition adapter requires an identified integer A,Z isotope")
    if not np.allclose(data.Yh, A, rtol=64*np.finfo(float).eps, atol=0):
        raise ValueError("decoded outer Yh must equal A for the verified raw unit group fraction")
    charge_residual = Z/A-data.Ye
    tolerance = 2e-5
    if np.any(~np.isfinite(charge_residual)) or np.any(abs(charge_residual) > tolerance):
        raise ValueError("outer pure-isotope electron fraction violates printed-table charge closure")
    corrected = replace(data, source_dir=None, Yh=np.ones(data.n),
                        row_status=np.full(data.n, 2, dtype=np.int64),
                        C3_applicable=np.ones(data.n, dtype=bool),
                        C4_applicable=np.ones(data.n, dtype=bool),
                        composition_available=True, composition_file_present=True)
    return corrected, dict(
        representation="verified_provider_outer_unit_baryon_fraction_adapter",
        raw_input_precondition="source hash verified; NP4 IDs10,11,0,1; one group; raw group value1",
        physical_domain="phase0 pure nuclear outer rows only; no inner or mixed-matter conversion",
        decoded_group_semantics="standard number-fraction decoder produced Yh=A",
        corrected_group_semantics="internal baryon fraction Yh=1",
        thermodynamic_fields_changed=[], preserved_electron_fraction=True,
        charge_closure_absolute_tolerance=tolerance,
        max_abs_charge_closure_residual=float(np.max(abs(charge_residual))),
        max_abs_baryon_closure_residual=0.,
        corrections=[dict(source_row=i, A=float(A[i]), Z=float(Z[i]),
                          decoded_Yh=float(data.Yh[i]), corrected_Yh=1.,
                          original_row_status=int(data.row_status[i]), corrected_row_status=2,
                          charge_closure_residual=float(charge_residual[i]))
                     for i in range(data.n)])


def project_euler_mu(data: EOSData) -> tuple[EOSData, dict]:
    """Preserve n, P and epsilon, explicitly replace Q3 by cold Euler mu."""
    _cold_slice(data)
    mu = (data.eps+data.p)/data.nb
    result = replace(data, source_dir=None, Q3=mu/data.mn-1.)
    return result, {
        "representation": "energy_and_pressure_preserving_euler_mu_projection",
        "changed_fields": ["Q3"], "source_energies_changed": False,
        "delta_mu_MeV": (mu-data.mub).tolist(),
        "maximum_absolute_delta_mu_MeV": float(np.max(abs(mu-data.mub))),
    }


@dataclass(frozen=True)
class SurfacePolytrope:
    """P=P0(n/n0)^gamma and epsilon=a*n+P/(gamma-1), including n=0."""
    n0: float
    p0: float
    epsilon0: float
    gamma: float = 5./3.

    def __post_init__(self):
        _positive_finite(self.n0, self.p0, self.epsilon0, self.gamma)
        if self.gamma <= 1 or self.gamma > 2:
            raise ValueError("surface exponent must lie in (1,2]")
        if self.rest_energy_MeV <= 0:
            raise ValueError("surface anchor implies nonpositive limiting chemical potential")
        if self.gamma*self.p0/(self.epsilon0+self.p0) > 1:
            raise ValueError("surface anchor violates causality")

    @classmethod
    def from_anchor(cls, n, p, epsilon, gamma=5./3.):
        return cls(float(n), float(p), float(epsilon), float(gamma))

    @property
    def rest_energy_MeV(self):
        return (self.epsilon0-self.p0/(self.gamma-1))/self.n0

    def evaluate(self, density):
        n = np.asarray(density, dtype=float)
        if np.any(~np.isfinite(n)) or np.any(n < 0) or np.any(n > self.n0):
            raise ValueError("surface density must be in [0,n0]")
        ratio = n/self.n0
        p = self.p0*ratio**self.gamma
        q = (self.p0/self.n0)*ratio**(self.gamma-1)
        mu = self.rest_energy_MeV+self.gamma*q/(self.gamma-1)
        epsilon = self.rest_energy_MeV*n+p/(self.gamma-1)
        cs2 = self.gamma*q/mu
        return {key: _scalar_or_array(value, density) for key, value in
                dict(n=n, p=p, epsilon=epsilon, mu=mu, cs2_over_c2=cs2).items()}

    def enthalpy(self, density):
        """Dimensionless integral h=int_0^P dP/(epsilon+P)=log(mu/a)."""
        n = np.asarray(density, dtype=float)
        self.evaluate(n)  # Validate the declared support, including zero.
        excess = self.gamma/(self.gamma-1)*(self.p0/self.n0)*(n/self.n0)**(self.gamma-1)
        return _scalar_or_array(np.log1p(excess/self.rest_energy_MeV), density)

    def table(self, density, *, template: EOSData, phase_id: int = 900):
        values = self.evaluate(density)
        return _thermo_table(template, values["n"], values["p"], values["epsilon"],
                             values["mu"], phase_id=phase_id)

    def contract(self):
        return dict(representation="anchored_surface_polytrope", n0=self.n0,
                    p0=self.p0, epsilon0=self.epsilon0, gamma=self.gamma,
                    rest_energy_MeV=self.rest_energy_MeV,
                    support_n_fm3=[0., self.n0], exact_surface=dict(n=0., p=0., epsilon=0.),
                    enthalpy="log(mu/rest_energy_MeV)",
                    composition="UNSPECIFIED_EFFECTIVE_BAROTROPE",
                    physical_scope="phenomenological low-density surface extension",
                    source_energy_shift="NONE: extension anchored to unchanged first source state")


class PiecewisePolytrope:
    """Pressure-primary, first-law-integrated log-linear pressure interpolant."""
    def __init__(self, density, pressure, first_energy_per_baryon, *, interval_kinds=None):
        self.nb = np.array(density, dtype=float, copy=True).reshape(-1)
        self.p = np.array(pressure, dtype=float, copy=True).reshape(-1)
        _positive_finite(self.nb, self.p, first_energy_per_baryon)
        if self.nb.size < 2 or self.p.shape != self.nb.shape:
            raise ValueError("at least two corresponding pressure-density knots required")
        if np.any(np.diff(self.nb) <= 0) or np.any(np.diff(self.p) <= 0):
            raise ValueError("pressure and density knots must increase strictly")
        self.log_span = np.log(self.nb[1:]/self.nb[:-1])
        self.gamma = np.log(self.p[1:]/self.p[:-1])/self.log_span
        self.interval_kinds = (list(interval_kinds) if interval_kinds is not None else
                               ["same_phase_pressure_interpolation"]*(self.nb.size-1))
        if len(self.interval_kinds) != self.nb.size-1 or any(not isinstance(k, str) or not k
                                                           for k in self.interval_kinds):
            raise ValueError("one explicit interval kind per source interval required")
        self.F = np.empty_like(self.nb)
        self.F[0] = float(first_energy_per_baryon)
        delta = self._integral(self.p[:-1]/self.nb[:-1], self.gamma-1, self.log_span)
        self.F[1:] = self.F[0]+np.cumsum(delta)
        self.mu = self.F+self.p/self.nb

    @staticmethod
    def _integral(q0, exponent, log_ratio):
        q0, exponent, log_ratio = np.broadcast_arrays(q0, exponent, log_ratio)
        value = exponent*log_ratio
        # exprel is evaluated by a local series near gamma=1 to avoid 0/0.
        factor = np.empty_like(value, dtype=float)
        near = abs(value) < 1e-5
        z = value[near]
        factor[near] = 1+z/2+z*z/6+z**3/24+z**4/120
        factor[~near] = np.expm1(value[~near])/value[~near]
        return q0*log_ratio*factor

    def evaluate(self, density):
        n = np.asarray(density, dtype=float)
        if np.any(~np.isfinite(n)) or np.any(n < self.nb[0]) or np.any(n > self.nb[-1]):
            raise ValueError("density outside pressure-primary interpolation support")
        index = np.clip(np.searchsorted(self.nb, n, side="right")-1, 0, self.nb.size-2)
        x = np.log(n/self.nb[index])
        p = self.p[index]*np.exp(self.gamma[index]*x)
        F = self.F[index]+self._integral(self.p[index]/self.nb[index], self.gamma[index]-1, x)
        mu = F+p/n
        epsilon = n*F
        cs2 = self.gamma[index]*p/(epsilon+p)
        return {key: _scalar_or_array(value, density) for key, value in
                dict(n=n, p=p, epsilon=epsilon, mu=mu, cs2_over_c2=cs2).items()}

    def validate(self):
        # On each power-law interval c_s^2 has a derivative of constant sign,
        # including gamma=1, so both one-sided endpoints give exact extrema.
        left = self.gamma*self.p[:-1]/(self.nb[:-1]*self.F[:-1]+self.p[:-1])
        right = self.gamma*self.p[1:]/(self.nb[1:]*self.F[1:]+self.p[1:])
        lo, hi = np.minimum(left, right), np.maximum(left, right)
        bad = (~np.isfinite(lo)) | (~np.isfinite(hi)) | (lo <= 0) | (hi > 1)
        return dict(status="PASS" if not np.any(bad) else "FAIL", interval_count=self.nb.size-1,
                    first_law="ANALYTIC_IDENTITY", extrema_method="both_one_sided_endpoints_exact_for_power_law",
                    min_cs2_over_c2=float(np.min(lo)), max_cs2_over_c2=float(np.max(hi)),
                    failed_intervals=np.flatnonzero(bad).tolist(),
                    intervals=[dict(index=i, kind=k, gamma=float(self.gamma[i]),
                                    min_cs2_over_c2=float(lo[i]), max_cs2_over_c2=float(hi[i]))
                               for i, k in enumerate(self.interval_kinds)])

    def contract(self):
        return dict(representation="pressure_primary_piecewise_power_law",
                    n_fm3=self.nb.tolist(), p_MeV_fm3=self.p.tolist(),
                    F_MeV=self.F.tolist(), gamma=self.gamma.tolist(),
                    interval_kinds=self.interval_kinds,
                    energy_anchor=dict(row=0, n=float(self.nb[0]), F_MeV=float(self.F[0])),
                    interpolation="P=P_i*(n/n_i)^gamma_i; F'=P/n^2",
                    isotope_gap_semantics="explicit effective regularization, composition undefined in gap")


def reconstruct_outer_pressure(data: EOSData, *, interval_kinds=None):
    """Return projected table, exact evaluator and a full physical correction ledger."""
    _cold_slice(data)
    if np.any(data.Q2 != 0) or np.any(data.Q5 != 0):
        raise ValueError("cold one-charge outer reconstruction requires zero entropy and Q5")
    model = PiecewisePolytrope(data.nb, data.p, data.eps[0]/data.nb[0],
                               interval_kinds=interval_kinds)
    validation = model.validate()
    if validation["status"] != "PASS":
        raise ValueError("pressure-primary outer reconstruction is not causal and stable")
    q_energy = model.F/data.mn-1.
    result = replace(data, source_dir=None, Q3=model.mu/data.mn-1.,
                     Q6=q_energy.copy(), Q7=q_energy.copy())
    delta_F = model.F-data.eps/data.nb
    report = dict(representation="public_derived_pressure_primary_outer", changed_fields=["Q3", "Q6", "Q7"],
                  unchanged_fields=["nb", "Q1", "Q2", "Q4", "Q5", "composition", "source_phase"],
                  physical_anchor="first raw source epsilon/n; no fitted or partner-dependent energy shift",
                  max_abs_delta_F_MeV=float(np.max(abs(delta_F))),
                  validation=validation,
                  corrections=[dict(source_row=i, n_fm3=float(data.nb[i]),
                                    delta_F_MeV=float(delta_F[i]),
                                    delta_epsilon_MeV_fm3=float(data.nb[i]*delta_F[i]),
                                    source_mu_MeV=float(data.mub[i]),
                                    reconstructed_mu_MeV=float(model.mu[i]),
                                    delta_mu_MeV=float(model.mu[i]-data.mub[i]))
                               for i in range(data.n)])
    return result, model, report


def _minimum_on_unit_interval(polynomial: Polynomial):
    points = [0., 1.]
    for value in polynomial.deriv().roots():
        if abs(value.imag) < 1e-12 and 0 < value.real < 1:
            points.append(float(value.real))
    values = np.asarray([polynomial(t) for t in points], dtype=float)
    j = int(np.argmin(values))
    return float(values[j]), points[j]


class HermiteInterval:
    """A named cubic P(mu) interval preserving both endpoint n, P and epsilon."""
    def __init__(self, n_left, p_left, epsilon_left, n_right, p_right, epsilon_right):
        _positive_finite(n_left, p_left, epsilon_left, n_right, p_right, epsilon_right)
        self.n_left, self.n_right = float(n_left), float(n_right)
        self.p_left, self.p_right = float(p_left), float(p_right)
        self.epsilon_left, self.epsilon_right = float(epsilon_left), float(epsilon_right)
        self.mu_left = (epsilon_left+p_left)/n_left
        self.mu_right = (epsilon_right+p_right)/n_right
        self.width = self.mu_right-self.mu_left
        if n_right <= n_left or p_right <= p_left or self.width <= 0:
            raise ValueError("bridge endpoints must increase in n, P and cold Euler mu")
        d0, d1, dp = self.width*n_left, self.width*n_right, p_right-p_left
        self.pressure = Polynomial([p_left, d0, 3*dp-2*d0-d1, -2*dp+d0+d1])
        self.number = self.pressure.deriv()/self.width
        self.susceptibility = self.number.deriv()/self.width
        self.chemical = Polynomial([self.mu_left, self.width])
        self.causal = self.chemical*self.susceptibility-self.number

    @classmethod
    def from_endpoints(cls, n_left, p_left, epsilon_left, n_right, p_right, epsilon_right):
        return cls(n_left, p_left, epsilon_left, n_right, p_right, epsilon_right)

    def validate(self):
        nn, xn = _minimum_on_unit_interval(self.number)
        ss, xs = _minimum_on_unit_interval(self.susceptibility)
        cc, xc = _minimum_on_unit_interval(self.causal)
        ee, xe = _minimum_on_unit_interval(self.chemical*self.number-self.pressure)
        # Bound roundoff in causal polynomial evaluation only. Strict positivity
        # of density and susceptibility is never replaced by a tolerance.
        # Scale the terms BEFORE their cancellation. A luminal EOS has an
        # exactly zero causal polynomial; scaling by the already-subtracted
        # residual would turn harmless cancellation roundoff into failure.
        before_subtraction = self.chemical*self.susceptibility
        cscale = max(1e-300, float(sum(abs(before_subtraction.coef))+sum(abs(self.number.coef))))
        bound = 64*np.finfo(float).eps*cscale
        passed = nn > 0 and ss > 0 and ee > 0 and cc >= -bound
        cs_min = cs_max = None
        if nn > 0 and ss > 0:
            denominator = self.chemical*self.susceptibility
            derivative_numerator = self.number.deriv()*denominator-self.number*denominator.deriv()
            points = [0., 1.]
            for value in derivative_numerator.roots():
                if abs(value.imag) < 1e-12 and 0 < value.real < 1:
                    points.append(float(value.real))
            speeds = [float(self.number(t)/denominator(t)) for t in points]
            cs_min, cs_max = min(speeds), max(speeds)
        return dict(status="PASS" if passed else "FAIL", first_law="ANALYTIC_IDENTITY",
                    min_n_fm3=nn, min_dn_dmu=ss, min_causal_margin=cc, min_epsilon=ee,
                    min_cs2_over_c2=cs_min, max_cs2_over_c2=cs_max,
                    causal_roundoff_bound=bound,
                    extrema_locations=dict(n=xn, dn_dmu=xs, causal=xc, epsilon=xe),
                    mean_density=(self.p_right-self.p_left)/self.width,
                    support_n_fm3=[self.n_left, self.n_right],
                    support_mu_MeV=[self.mu_left, self.mu_right])

    def evaluate_mu(self, chemical_potential):
        mu = np.asarray(chemical_potential, dtype=float)
        if np.any(~np.isfinite(mu)) or np.any(mu < self.mu_left) or np.any(mu > self.mu_right):
            raise ValueError("chemical potential outside bridge support")
        t = (mu-self.mu_left)/self.width
        p, n = self.pressure(t), self.number(t)
        epsilon = mu*n-p
        cs2 = n/(mu*self.susceptibility(t))
        return {key: _scalar_or_array(value, chemical_potential) for key, value in
                dict(n=n, p=p, epsilon=epsilon, mu=mu, cs2_over_c2=cs2).items()}

    def table(self, count: int, *, template: EOSData, phase_id: int = 901):
        if isinstance(count, bool) or int(count) != count or count < 3:
            raise ValueError("at least three bridge samples required")
        if self.validate()["status"] != "PASS":
            raise ValueError("bridge violates continuous stability or causality")
        mu = np.linspace(self.mu_left, self.mu_right, int(count))
        values = self.evaluate_mu(mu)
        # Analytic endpoint identities may lose a last bit in polynomial
        # evaluation. Preserve the exact physical input states at both ends.
        values["n"][[0, -1]] = [self.n_left, self.n_right]
        values["p"][[0, -1]] = [self.p_left, self.p_right]
        values["epsilon"][[0, -1]] = [self.epsilon_left, self.epsilon_right]
        return _thermo_table(template, values["n"], values["p"], values["epsilon"],
                             values["mu"], phase_id=phase_id)

    def contract(self):
        return dict(representation="cubic_pressure_in_euler_chemical_potential",
                    mu_left_MeV=self.mu_left, mu_width_MeV=self.width,
                    pressure_coefficients_in_t=self.pressure.coef.tolist(),
                    endpoint_states=[dict(n=self.n_left, p=self.p_left, epsilon=self.epsilon_left),
                                     dict(n=self.n_right, p=self.p_right, epsilon=self.epsilon_right)],
                    t_definition="(mu-mu_left)/mu_width", composition="UNSPECIFIED_IN_BRIDGE",
                    validation=self.validate())


def certify_piecewise_hermite(data: EOSData, *, sharp_interfaces: Sequence[int] = ()):
    """Check smooth intervals and caller-declared coexistence doublets.

    A phase label alone does not authorize a Maxwell boundary. The caller
    supplies the left-row indices from the construction receipt explicitly.
    """
    _cold_slice(data)
    sharp = set(sharp_interfaces)
    if any(isinstance(i, bool) or int(i) != i or i < 0 or i >= data.n-1 for i in sharp):
        raise ValueError("sharp interface index outside table")
    results = []
    for i in range(data.n-1):
        if i in sharp:
            n0, n1 = data.nb[i:i+2]
            p0, p1 = data.p[i:i+2]
            e0, e1 = data.eps[i:i+2]
            mu0, mu1 = (data.eps[i:i+2]+data.p[i:i+2])/data.nb[i:i+2]
            pscale = max(abs(p0), abs(p1), 1e-300)
            mscale = max(abs(mu0), abs(mu1), 1e-300)
            latent = abs((e1-e0)-mu0*(n1-n0))
            passed = (n1 > n0 and abs(p1-p0) <= 2e-12*pscale
                      and abs(mu1-mu0) <= 2e-12*mscale
                      and latent <= 2e-12*max(abs(e0), abs(e1), 1e-300))
            results.append(dict(index=i, kind="declared_sharp_interface",
                                status="PASS" if passed else "FAIL",
                                pressure_difference=float(p1-p0), mu_difference=float(mu1-mu0),
                                latent_identity_absolute_residual=float(latent),
                                interpolation="FORBIDDEN_ACROSS_DENSITY_JUMP"))
            continue
        try:
            interval = HermiteInterval(data.nb[i], data.p[i], data.eps[i],
                                       data.nb[i+1], data.p[i+1], data.eps[i+1])
            result = interval.validate()
            result["contract"] = interval.contract()
        except ValueError as exc:
            result = dict(status="FAIL", reason=str(exc))
        results.append(dict(index=i, **result))
    failures = [r["index"] for r in results if r["status"] != "PASS"]
    return dict(status="PASS" if not failures else "FAIL", interval_count=len(results),
                failed_intervals=failures, intervals=results,
                source_mu_euler_residual_max_MeV=float(np.max(abs(data.mub-(data.eps+data.p)/data.nb))))


@dataclass
class EOSPiece:
    data: EOSData
    label: str
    origin: str
    source_rows: Sequence[int] | None = None


def assemble_shared_pieces(pieces: Sequence[EOSPiece]):
    """Join shared physical endpoint rows, preserving source rows and masks.

    All gaps must already be represented by an explicit piece. The output
    records each piece's inclusive row range for its analytic evaluator; a
    shared endpoint may have the neighboring piece's row provenance.
    """
    pieces = list(pieces)
    if len(pieces) < 2 or len({p.label for p in pieces}) != len(pieces):
        raise ValueError("at least two uniquely labelled pieces required")
    first = pieces[0].data
    tables = []
    for piece in pieces:
        _cold_slice(piece.data)
        if piece.data.Il != first.Il or piece.data.temperature_mev != first.temperature_mev:
            raise ValueError("piece slice conventions differ")
        if piece.source_rows is not None and len(piece.source_rows) != piece.data.n:
            raise ValueError("source row map has wrong length")
        tables.append(rebase_header(piece.data, mn=first.mn, mp=first.mp))
    rows, provenance, spans, phases = [], [], [], {}
    boundary_reports = []
    priority = lambda origin: 2 if origin in ("public_derived", "generated", "native_generated") else 1
    for piece_index, (piece, table) in enumerate(zip(pieces, tables)):
        start = len(rows)
        if piece_index:
            previous = tables[piece_index-1]
            residuals = {name: float(abs(getattr(previous, name)[-1]-getattr(table, name)[0]))
                         for name in ("nb", "p", "eps", "mub")}
            for name in ("nb", "p", "eps", "mub"):
                a, b = getattr(previous, name)[-1], getattr(table, name)[0]
                if abs(a-b) > 128*np.finfo(float).eps*max(abs(a), abs(b), 1e-300):
                    raise ValueError("shared endpoint mismatch in "+name)
            boundary_reports.append(dict(left=pieces[piece_index-1].label, right=piece.label,
                                         law="SHARED_STATE_P_EPS_N_MU_CONTINUOUS", absolute_residuals=residuals))
            start -= 1
            if priority(piece.origin) > priority(pieces[piece_index-1].origin):
                rows.pop(); provenance.pop()
                first_local = 0
            else:
                first_local = 1
        else:
            first_local = 0
        for j in range(first_local, table.n):
            phase_key = (piece.label, int(table.Iphase[j]))
            if phase_key not in phases:
                phases[phase_key] = len(phases)+1
            row = {key: getattr(table, key)[j] for key in ROW_FIELDS}
            row["Iphase"] = phases[phase_key]
            rows.append(row)
            provenance.append(dict(origin=piece.origin, piece=piece.label,
                                   source_row=None if piece.source_rows is None else int(piece.source_rows[j]),
                                   piece_row=j))
        spans.append(dict(label=piece.label, origin=piece.origin, first_row=start,
                          last_row=len(rows)-1, endpoint_policy="shared physical state; source-priority row retained"))
    fields = {key: np.asarray([row[key] for row in rows], dtype=np.asarray(getattr(first, key)).dtype)
              for key in ROW_FIELDS}
    combined = replace(tables[0], source_dir=None,
                       composition_file_present=any(t.composition_file_present for t in tables),
                       composition_available=any(t.composition_available for t in tables),
                       beta_equilibrated=all(t.beta_equilibrated for t in tables), **fields)
    _cold_slice(combined)
    return combined, dict(status="ASSEMBLED_SHARED_ENDPOINT_PIECES", row_count=combined.n,
                          family="reference_assisted_hybrid", provenance=provenance, pieces=spans,
                          boundaries=boundary_reports,
                          phase_map=[dict(piece=k[0], source_phase=k[1], output_phase=v)
                                     for k, v in phases.items()],
                          composition_interpolation="undefined inside derived pieces and isotope-gap regularizations",
                          consumer_rule="use the supplied piece contracts; no opaque spline over all output rows")


class WholeCrustBarotrope:
    """Continuous pressure consumer for the explicit A piece contracts.

    The finite-positive-density export is a sample of this representation.
    This evaluator reaches the analytic vacuum surface P=0 and uses exactly
    the recorded region laws, including caller-declared sharp interfaces.
    It supplies equilibrium barotropic quantities, not composition or a
    frozen-composition perturbation closure.
    """
    _state_keys = ("n", "p", "epsilon", "mu", "cs2_over_c2", "enthalpy")

    @classmethod
    def from_contract(cls, record: dict):
        self = cls()
        surface = record["surface"]
        self.surface = SurfacePolytrope(surface["n0"], surface["p0"], surface["epsilon0"], surface["gamma"])
        outer = record["outer"]
        self.outer = PiecewisePolytrope(outer["n_fm3"], outer["p_MeV_fm3"], outer["F_MeV"][0],
                                        interval_kinds=outer["interval_kinds"])
        if not np.allclose(self.outer.F, outer["F_MeV"], rtol=64*np.finfo(float).eps, atol=0):
            raise ValueError("outer contract energy knots differ from its first-law law")
        if self.outer.validate()["status"] != "PASS":
            raise ValueError("outer pressure law is not causal and stable")
        left, right = record["lower_bridge"]["endpoint_states"]
        self.bridge = HermiteInterval(left["n"], left["p"], left["epsilon"],
                                       right["n"], right["p"], right["epsilon"])
        if self.bridge.validate()["status"] != "PASS":
            raise ValueError("lower bridge is not causal and stable")
        upper = record["upper"]
        self.upper_n = np.asarray(upper["n"], dtype=float)
        self.upper_p = np.asarray(upper["p"], dtype=float)
        self.upper_epsilon = np.asarray(upper["epsilon"], dtype=float)
        if (self.upper_n.ndim != 1 or self.upper_n.size < 2
                or self.upper_p.shape != self.upper_n.shape or self.upper_epsilon.shape != self.upper_n.shape):
            raise ValueError("upper contract requires corresponding one-dimensional state arrays")
        _positive_finite(self.upper_n, self.upper_p, self.upper_epsilon)
        if np.any(np.diff(self.upper_n) <= 0):
            raise ValueError("upper densities must increase")
        self.upper_mu = (self.upper_epsilon+self.upper_p)/self.upper_n
        self.sharp_interfaces = set(upper.get("sharp_interfaces", []))
        if any(isinstance(i, bool) or int(i) != i or i < 0 or i >= self.upper_n.size-1
               for i in self.sharp_interfaces):
            raise ValueError("upper sharp interface index outside table")
        if any(i+1 in self.sharp_interfaces for i in self.sharp_interfaces):
            raise ValueError("adjacent sharp interfaces require a separate multiphase contract")
        self.upper_intervals = []
        self.sharp_pressure_bands = {}
        self._upper_search_pressure = self.upper_p.copy()
        for i in range(self.upper_n.size-1):
            if i in self.sharp_interfaces:
                p0, p1 = self.upper_p[i:i+2]
                mu0, mu1 = self.upper_mu[i:i+2]
                if (abs(p1-p0) > 2e-12*max(abs(p0), abs(p1))
                        or abs(mu1-mu0) > 2e-12*max(abs(mu0), abs(mu1))):
                    raise ValueError("declared upper interface lacks common pressure and chemical potential")
                self.upper_intervals.append(None)
                # Admission accepts a finite equal-P encoding discrepancy.
                # Its entire observed signed band belongs to the one-sided
                # interface selection, never to a smooth root solver.
                pad = 8*np.finfo(float).eps*max(abs(p0), abs(p1))
                self.sharp_pressure_bands[i] = (min(p0, p1)-pad, max(p0, p1)+pad)
                # This is a search-coordinate convention at an admitted
                # coexistence doublet, not a change of either source row.
                self._upper_search_pressure[i+1] = p0
                continue
            interval = HermiteInterval(self.upper_n[i], self.upper_p[i], self.upper_epsilon[i],
                                       self.upper_n[i+1], self.upper_p[i+1], self.upper_epsilon[i+1])
            if interval.validate()["status"] != "PASS":
                raise ValueError(f"upper smooth interval {i} is not causal and stable")
            self.upper_intervals.append(interval)
        if np.any(np.diff(self._upper_search_pressure) < 0):
            raise ValueError("upper pressure order differs from its phase-boundary contract")
        self._check_same_state(
            (self.surface.n0, self.surface.p0, self.surface.epsilon0),
            (self.outer.nb[0], self.outer.p[0], self.outer.nb[0]*self.outer.F[0]), "surface/outer")
        self._check_same_state(
            (self.outer.nb[-1], self.outer.p[-1], self.outer.nb[-1]*self.outer.F[-1]),
            (self.bridge.n_left, self.bridge.p_left, self.bridge.epsilon_left), "outer/lower bridge")
        self._check_same_state(
            (self.bridge.n_right, self.bridge.p_right, self.bridge.epsilon_right),
            (self.upper_n[0], self.upper_p[0], self.upper_epsilon[0]), "lower bridge/upper")
        self.pressure_domain = (0., float(self.upper_p[-1]))
        self.density_domain = (0., float(self.upper_n[-1]))
        return self

    @staticmethod
    def _check_same_state(left, right, label):
        for a, b in zip(left, right):
            if abs(a-b) > 128*np.finfo(float).eps*max(abs(a), abs(b), 1e-300):
                raise ValueError(label+" contracts do not share the same physical state")

    @staticmethod
    def _interval_at_pressure(interval: HermiteInterval, pressure: float):
        from scipy.optimize import brentq
        if pressure == interval.p_left:
            t = 0.
        elif pressure == interval.p_right:
            t = 1.
        else:
            def residual(t):
                if t == 0.:
                    return interval.p_left-pressure
                if t == 1.:
                    return interval.p_right-pressure
                return float(interval.pressure(t))-pressure
            t = brentq(residual, 0., 1., xtol=5e-15, rtol=4*np.finfo(float).eps)
        mu = interval.mu_left+t*interval.width
        n = float(interval.number(t))
        epsilon = mu*n-pressure
        if t == 0.:
            n, epsilon = interval.n_left, interval.epsilon_left
        elif t == 1.:
            n, epsilon = interval.n_right, interval.epsilon_right
        return dict(n=n, p=pressure, epsilon=epsilon, mu=float(mu),
                    cs2_over_c2=float(interval.number(t)/(mu*interval.susceptibility(t))))

    def _outer_at_pressure(self, pressure):
        index = int(np.clip(np.searchsorted(self.outer.p, pressure, side="right")-1,
                            0, self.outer.nb.size-2))
        if pressure == self.outer.p[index]:
            n, F = self.outer.nb[index], self.outer.F[index]
        elif pressure == self.outer.p[index+1]:
            n, F = self.outer.nb[index+1], self.outer.F[index+1]
        else:
            x = np.log(pressure/self.outer.p[index])/self.outer.gamma[index]
            n = self.outer.nb[index]*np.exp(x)
            F = self.outer.F[index]+self.outer._integral(
                self.outer.p[index]/self.outer.nb[index], self.outer.gamma[index]-1, x)
        epsilon = n*F
        mu = F+pressure/n
        return dict(n=float(n), p=pressure, epsilon=float(epsilon), mu=float(mu),
                    cs2_over_c2=float(self.outer.gamma[index]*pressure/(epsilon+pressure)))

    def _sharp_state(self, index, pressure, side):
        row = index if side == "low" else index+1
        adjacent = index-1 if side == "low" else index+1
        if not 0 <= adjacent < len(self.upper_intervals) or self.upper_intervals[adjacent] is None:
            raise ValueError("one-sided sound speed lacks an adjacent smooth interval")
        interval = self.upper_intervals[adjacent]
        endpoint_fraction = 1. if side == "low" else 0.
        mu = (self.upper_epsilon[row]+pressure)/self.upper_n[row]
        return dict(n=float(self.upper_n[row]), p=pressure, epsilon=float(self.upper_epsilon[row]),
                    mu=float(mu), cs2_over_c2=float(interval.number(endpoint_fraction)/(
                        mu*interval.susceptibility(endpoint_fraction))))

    def _state_one(self, pressure: float, transition_side: str):
        if pressure <= self.surface.p0:
            n = self.surface.n0*(pressure/self.surface.p0)**(1/self.surface.gamma)
            # Retain the exact anchor at equality; exponent evaluation can
            # otherwise wander outside the closed analytic support by 1 ULP.
            if pressure == self.surface.p0:
                n = self.surface.n0
            result = self.surface.evaluate(n)
            result["p"] = pressure
            result["enthalpy"] = self.surface.enthalpy(n)
            return result
        if pressure <= self.outer.p[-1]:
            result = self._outer_at_pressure(pressure)
        elif pressure <= self.bridge.p_right:
            result = self._interval_at_pressure(self.bridge, pressure)
        else:
            for index in self.sharp_interfaces:
                lower, upper = self.sharp_pressure_bands[index]
                if lower <= pressure <= upper:
                    result = self._sharp_state(index, pressure, transition_side)
                    break
            else:
                index = int(np.clip(np.searchsorted(self._upper_search_pressure, pressure, side="right")-1,
                                    0, self.upper_n.size-2))
                interval = self.upper_intervals[index]
                if interval is None:
                    raise ValueError("pressure selects a density jump without a declared side")
                result = self._interval_at_pressure(interval, pressure)
        result["enthalpy"] = float(np.log(result["mu"]/self.surface.rest_energy_MeV))
        return result

    def state_at_pressure(self, pressure, *, transition_side="high"):
        """Evaluate scalar or array P, choosing a one-sided state at sharp P_t.

        Returned fields are n, p, epsilon, mu, c_s^2/c^2 and dimensionless
        enthalpy. P=0 is the analytic vacuum state. Values beyond the actual
        highest core pressure raise; there is no high-density extrapolation.
        """
        if transition_side not in ("low", "high"):
            raise ValueError("transition_side must be 'low' or 'high'")
        values = np.asarray(pressure, dtype=float)
        if (np.any(~np.isfinite(values)) or np.any(values < 0)
                or np.any(values > self.pressure_domain[1])):
            raise ValueError("pressure outside whole-crust contract support")
        states = [self._state_one(float(p), transition_side) for p in values.reshape(-1)]
        return {key: _scalar_or_array(np.asarray([s[key] for s in states]).reshape(values.shape), pressure)
                for key in self._state_keys}
