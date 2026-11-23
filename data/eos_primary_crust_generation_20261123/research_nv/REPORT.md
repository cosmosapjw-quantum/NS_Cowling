# Negele--Vautherin primary-source recovery and native generation

Date: 2026-11-23. Scope: the previously pressure-missing NV source. The native
generator is `NV1974_author_fit`, not a new Hartree--Fock solution and not a
claim of byte/numerical identity with the unavailable 1973 original pressure.

## Source recovery and authority

- Original NV paper: J. W. Negele and D. Vautherin, *Nuclear Physics A* **207**
  (1973), 298--320, DOI [10.1016/0375-9474(73)90349-7](https://doi.org/10.1016/0375-9474(73)90349-7).
  DOI verified against the publisher, INSPIRE record 69604, and OSTI record
  4480535. Public full text was not recovered during this task.
- Recovered primary author publication: J. W. Negele, *The Equation of State of
  Matter at Sub-Nuclear Density*, IAU Symposium **53** (1974), 1--25,
  DOI [10.1017/S0074180900099885](https://doi.org/10.1017/S0074180900099885).
  This supplies a complete author analytic EOS expression and coefficients.
- [Publisher PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C5DF32C365CA76E65A4530174B7A4249/S0074180900099885a.pdf/div-class-title-the-equation-of-state-of-matter-at-sub-nuclear-density-div.pdf)
  downloaded to private execution scratch:
  `/tmp/ns_eos_primary_crust_20261123/sources/nv/Negele_1974_IAU53.pdf`.
  SHA256 `328d47923c38d0531ce409c83f827114bdc751f207b8b893553537a45226fa71`.
  The full PDF is not included in the repository.
- Table II, p.20: 11 cell configurations and energies; Table III, p.23:
  ground-state and uniform-gas fit coefficients; p.24: energy fit and
  thermodynamic derivatives. These three pages were rasterized and visually
  checked; the exponent and signs were not taken from OCR alone.
- Section 4, pp.14--17: cold locally neutral matter, electrons included; the
  tabulated neutron/proton chemical potentials are last occupied HF orbital
  eigenvalues. The beta-stability expression on p.20 includes rest masses
  separately, confirming that the tabulated chemical potentials omit them.
- The author explicitly warns on p.24 that second derivatives of an energy fit
  to 11 points amplify fit noise. No coefficient refitting was performed here.

The immutable pre-existing input remains
`data/eos_missing_data_recalc_20261119/external_original/1973-Negele-Vautherin.dat`,
SHA256 `0a00852b60a02f4b6dfdc8f6305c343ca188c41f9cf0f2eb31983294520efe1e`.
Its density column is in cm^-3, its pressure is missing, and it has 11 rows and
9 columns. A new factual transcription lives in `nv1974_primary_source.dat`.

## Formula and computational contract

Retain physical rest energy explicitly: let m_n c^2 be the neutron rest energy,
n the baryon number density, and x = ln[n/(10^-4 fm^-3)]. Define

```text
Q(x) = c1 + c2*x + c3*x^2 + c4*x^3 + c5*x^4 + c6*x^5 + c7*x^6
e(n) = E_T - m_n*c^2 = c0 + exp[Q(x)] MeV
```

The primary coefficients are:

| i | Ground-state c_i |
|---|---:|
| 0 | -4.0 |
| 1 | 2.8822899e-1 |
| 2 | 5.9150523e-1 |
| 3 | 9.0185940e-2 |
| 4 | -1.1025614e-1 |
| 5 | 2.9377479e-2 |
| 6 | -3.2618465e-3 |
| 7 | 1.3543555e-4 |

c0 is in MeV; the polynomial is dimensionless in the source's MeV energy
convention. Since d/dn = n^-1 d/dx, let F = exp(Q) MeV,
A = F Q', and B = F[(Q')^2+Q'']. Then the directly derived formulas are

```text
epsilon = n (m_n*c^2 + c0 + F)
P       = n A
mu_b    = m_n*c^2 + c0 + F + A
dP/dn   = A + B
c_s^2/c^2 = (A+B)/mu_b
Gamma   = (A+B)/A .
```

Thus d epsilon/dn = mu_b and epsilon+P=n mu_b hold analytically. There is no
pressure finite differencing in the production path.

Production evaluates Q, Q', Q'' together by Horner recurrences. Cost is O(N)
for N generated points and O(N) storage for diagnostic output. Per-point work
contains six polynomial recurrence steps plus logarithm/exponential evaluation;
no iterative HF solver, quadrature, Python numerical kernel, or thread pool is
required. The parent campaign parallelizes independent source/core/method jobs.

The generator is deliberately restricted to the 11-point source-table interval
2.79e-4 <= n <= 7.89e-2 fm^-3. It does not extend to densities suggested only by
the author's plotted curves. The numerical neutron rest energy
939.56542052 MeV and proton header value 938.27208816 MeV are the existing
project convention, explicitly recorded; the paper writes the neutron mass
symbolically. The fitted energy already includes electrons, so no electron
term is added again. No interpolated composition is fabricated.

## Implementation and evidence

Implemented files:

- `fortran_eos/src/compose_eos_nv_primary.f90`
- `python_frontend/tests/test_nv_primary_native.py`
- `data/eos_primary_crust_generation_20261123/research_nv/nv1974_primary_source.dat`

Public entry point:

```fortran
call generate_nv_primary_native(source_path,n_points,eos,output_dir,status)
```

The output directory must already exist; success is `status='ok'`.
All failure paths return a diagnostic status instead of terminating the parent
program. Allowed grid sizes are 3..100000. Existing evidence files are not
overwritten. The strict source parser checks exact row/column counts, trailing
data, finite numbers, primary coefficient/table identity, density convention,
and the declared mass convention. A test caught a real exponent-overflow
failure under IEEE traps; lexical exponent validation now rejects it before
numeric conversion.

Outputs: `nv_primary_states.csv`, `nv_primary_source_comparison.csv`,
`NV_PRIMARY_GENERATION.json`. Parent integration owns CompOSE output,
generation dispatch, matching campaigns, source catalog and publication.

Verification command:

```text
python3 -m pytest -q --basetemp=/tmp/ns_eos_primary_crust_20261123/nv_build/pytest_final python_frontend/tests/test_nv_primary_native.py
```

Actual result: **30 passed in 2.11 s**, exit code 0. The native test compilation
uses gfortran, F2018, O2, array/runtime checking, and invalid/zero/overflow traps.
Python only builds, launches, hashes and reads results.

An independent Fortran real128 five-point finite-difference oracle
differentiates the primary energy expression at 31 interior densities. It is
separate from the production analytic derivative/Horner implementation.
A 4097-point domain scan checks pressure, compressibility, causal sound speed
and the Euler identity.

| Check | Observed maximum/range |
|---|---:|
| Relative P difference against real128 energy derivative | 1.3544720900e-14 |
| Absolute mu difference | 1.1368683772e-13 MeV |
| Relative sound-speed-squared difference | 1.7874590696e-14 |
| Euler arithmetic residual | 2.2737367544e-13 MeV |
| c_s^2/c^2 over 4097 points | 8.6083008019e-4 .. 9.1537622492e-3 |
| Polynomial absolute-sum/absolute-value cancellation ratio | <= 60.7366359955 |

These checks establish implementation fidelity to the published analytic
representation. They do not establish that the representation reproduces all
original HF pointwise quantities.

## Source fidelity and a physically relevant limitation

The 1974 table's row at n=0.00577 fm^-3 prints e=1.966 MeV; the pre-existing
1973/NUDA file has 1.996 MeV. Both are retained in the comparison diagnostic.
The original file was not modified.

At the 11 source densities, the generated analytic curve has:

| Diagnostic | Maximum |
|---|---:|
| Absolute energy offset difference from 1974 Table II | 0.1154114073 MeV |
| Absolute difference from historical 1973/NUDA energy | 0.1154114073 MeV |
| Analytic mu offset minus Table II HF neutron eigenvalue | 1.3189488036 MeV |
| Relative difference between P_fit and n(mu_n_table-e_table), normalized by P_fit | 0.2720378940 |

The Euler reconstruction from rounded Table II columns is retained as a
diagnostic, with a rounding-only envelope. It is not asserted to be the
author's pressure or an exact fit constraint: the tabulated quantities come
from finite-cell HF minimization, whereas the published pressure differentiates
the author's smooth energy fit.

At the lowest sampled density the fitted mu_b-m_n*c^2 is
**-0.2168231026 MeV** despite a tabulated neutron eigenvalue of +0.2 MeV.
Seven of 257 grid points lie below the free-neutron threshold under the fitted
chemical potential. This is a limitation of the published representation near
drip, not floating-point error, and it has not been hidden by altering its
coefficients or cutting the table. Mechanical stability and first-law closure
do not by themselves validate the inner-crust phase interpretation there.
The parent must keep this distinction in its matching/coverage conclusions.

## Claim status and next action

- **literature-supported:** recovered primary author expression, coefficients,
  thermodynamic convention, and the Table II discrepancy.
- **derived:** closed analytic P, mu, compressibility and sound-speed formulas.
- **implementation-verified / numerically checked:** generated Fortran
  barotrope and 30 tests, with the quantitative scope above.
- **not claimed:** fresh NV Hartree--Fock minimization; exact identity with
  original 1973 pressure; full outer+inner crust coverage; all physical
  composition conditions; successful core matching before the parent's run.
- **source access unresolved:** the original 1973 NPA full text. User provision
  of DOI **10.1016/0375-9474(73)90349-7** would permit direct comparison of the
  1973 and 1974 coefficients/tables. It is not required to run this separately
  labeled, already recovered author-fit generator.

The parent can now use the new source identity in the remaining matching
campaign while retaining the original pressure-missing record and the above
scientific distinction.
