# Correction: CompOSE normalization must preserve gravitating energy

The pre-audit implementation at `d90d3dca033f759c79ed034b3b8ff939acb43145`
added `dm = mn_core - mn_crust` to every retained crust free energy per baryon
and chemical potential, and added `n_B * dm` to its energy density. The code,
regression tests and machine-readable mass contract called this a
representation-only operation. That interpretation was incorrect.

For CompOSE reference energies `m_old` and `m_ref` (in MeV, i.e. rest mass
multiplied by `c^2`), normalization alone means

```
1 + Q6_ref = (m_old / m_ref) * (1 + Q6_old)
1 + Q7_ref = (m_old / m_ref) * (1 + Q7_old)
1 + Q3_ref = (m_old / m_ref) * (1 + Q3_old)
Q4_ref = (m_old / m_ref) * Q4_old
Q5_ref = (m_old / m_ref) * Q5_old
```

The density, pressure, free energy, total energy density and physical chemical
potentials are unchanged. This agrees with the normalization formula already
recorded in `docs/MAXWELL_MASS_AUDIT_20261114.md`, section 8.1. The source tables
are not changed by this correction.

The old operation preserves the derivative of free energy per baryon but does
not preserve gravitating energy or sound speed. Writing `s = (c_s/c)^2`,

```
s_old = (dP/dn_B) / (d epsilon/dn_B)
s_shifted = (dP/dn_B) / (d epsilon/dn_B + dm)
```

Thus a reversible additive energy transformation is not proof of physical
identity. It also changes coexistence chemical potentials, the chosen matching
window, and potentially the existence of a positive density jump.

## Implementation and report semantics

Native `find_window` and `match_crust_core` now decode each input with its own
reference mass. They retain its physical energy through endpoint selection,
bridge construction and output assembly. The reference Python implementation
uses the same rule. Output Q columns are encoded with the core reference.

The legacy `dm_mass_gauge_meV` and frontend `dm_crust` report fields now mean
the **applied physical energy-per-baryon shift**, which is zero. The frontend
separately reports `normalization_reference_delta_mev = mn_core - mn_crust`;
this difference is not applied to the physical EOS. No sharp-Maxwell mass gate
is relaxed and no source-dependent constituent-mass correction is introduced.

The existing output Euler re-derivation `mu_E = F + P/n_B` is separate from mass
normalization. The reported retained-source Euler residual continues to disclose
that change; source Q3 is not silently declared numerically identical to it.

Native matching now also rejects nonfinite temperature/charge-grid metadata
and, for `Il=1`, nonzero Q5 beyond the existing `1e-8` beta-equilibrium bound,
even when an input caller sets `beta_equilibrated = true`.

## Regression evidence and historical scope

The newly added Python physical-energy regression failed on the old code:
all 300 retained GPPVA_FSU2H rows changed by `+0.56533 MeV` per baryon. It passes
after the correction for CSS and all four cubic/quintic methods. Native
re-encoding invariance regressions independently compare physically identical
inputs stored with reference masses differing by `0.01 MeV`; all five failed
against the old executable. Final native verification is recorded by the
parent audit run, not inferred from these local Python results.

The old SPG_M1+RG_SLY4 positive CSS test depended on the energy shift. Without
it, the Python-selected candidate has `n_start = 0.06878987 fm^-3` below
`n_crust = 0.074131024 fm^-3`, so the positive-jump gate correctly rejects it.
That rejection remains an explicit regression. Positive interface tests use
the existing GPPVA_FSU2H+GM_GM1 benchmark. This is not a no-go theorem for other
SPG/SLy4 windows, and the older survivor census cannot establish the corrected
implementation's survivor set. Historical receipts are retained, while the
new census must use the corrected code and record the resulting changes.
