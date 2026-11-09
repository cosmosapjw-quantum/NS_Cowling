# Thermodynamic Surface Continuation Implementation Plan

> **For agentic workers:** Use the host's available task-by-task implementation workflow. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Construct thermodynamically consistent positive-state low-density continuations for the six-method pilot and quantify their effect on fixed-central-pressure TOV mass and radius.

**Architecture:** A pure-Python continuation component reconstructs the cutoff tangent from the first native cubic-Hermite P(mu) cell, defines two C1 curvature families, and obtains baryon density by integrating the cold first law. It prepends positive tail knots to the frozen 19-column normalized table without changing any source knot. A dedicated native-Fortran TOV runner compares the original finite-cutoff background with a tangent tail at two lower cutoffs and two distinct curvature-deformation families.

**Tech Stack:** Python 3.12, NumPy, pytest, GNU Fortran, frozen v5 `cell_solve task=star`, GitHub Actions.

## Global Constraints

- Do not modify the frozen v5 sources, projected source EOS knots, release summary, old spectrum table, or PR #3.
- Preserve all original projected rows and their ordering byte-for-byte at the CSV cell-value level; only prepend positive continuation rows.
- Keep all native table primary fields strictly positive and strictly increasing: `n_b`, `P`, `epsilon`, and `mu_b`. The native reader rejects literal zeros and extrapolation.
- Enforce the cold first law `d epsilon = (epsilon + P) d n_b / n_b`, `mu_b=(epsilon+P)/n_b`, and hence `dP=n_b dmu_b`.
- Use the exact right-sided derivative of the native cubic-Hermite `P(mu)` cell at the cutoff. Do not substitute a P(epsilon) spline derivative.
- The baseline and shape-deformed continuations must have identical `P`, `epsilon`, `n_b`, `mu_b`, and `dP/depsilon` at the source cutoff.
- Use the same continuation family and parameter values for all six matching methods; first five projected knots were verified identical across them.
- TOV comparisons use each method's frozen reference central pressure, `r0=10 m`, `r_limit=1e8 m`, `rtol=1e-9`, and `step_cap=0.05`. First reproduce the existing finite-cutoff baseline before interpreting extended runs.
- This phase validates equilibrium barotropic thermodynamics and TOV structure only. It does not classify radial stability, calculate radial `omega_0^2`, or admit new Cowling f/p1 spectra.

---

### Task 1: Thermodynamic continuation model and focused tests

**Files:**
- Create: `python_frontend/compose_cowling/surface_continuation.py`
- Test: `python_frontend/tests/test_surface_continuation.py`

**Interfaces:**
- `ThermodynamicSurfaceContinuation.from_projected_rows(rows, *, alpha=0.0, family="power_law")`
- `state_at_energy(epsilon_mev_fm3) -> dict[str, float]`
- `state_many(epsilon_mev_fm3) -> dict[str, numpy.ndarray]`
- `build_extended_projected_rows(rows, *, alpha, family, tail_points, min_energy_ratio)`
- `write_extended_projected_csv(source_path, output_path, *, alpha, family, tail_points, min_energy_ratio) -> receipt`

- [x] **Step 1: Add the focused failing tests**

Test cutoff C1 matching; first-law derivative closure; `mu*n=epsilon+P`; exponential and polynomial curvature families preserving the same cutoff tangent; positivity, monotonicity and subluminal `cs2`; low-pressure table ordering; unchanged original source rows; native 19-column header; and SHA receipts.

- [x] **Step 2: Verify the relevant failure**

Run: `PYTHONPATH=. pytest -q python_frontend/tests/test_surface_continuation.py`
Expected: collection fails because `surface_continuation.py` is not yet implemented.

- [x] **Step 3: Implement the minimum thermodynamic model**

For `x=epsilon/epsilon_cut`, use `P=P_cut*x**gamma_s*F(x)`, where `gamma_s=(epsilon_cut/P_cut)*cs2_cut`. Support three forms: `F=1` (local tangent power law), `F=exp(alpha*(1-x)**2)`, and `F=1+alpha*(1-x)**2`. Keep the power exponent `gamma_s=(epsilon_cut/P_cut)*cs2_cut` distinct from the equilibrium adiabatic index `Gamma_eq=((epsilon_cut+P_cut)/P_cut)*cs2_cut`. Each `F` satisfies `F(1)=1` and `F'(1)=0`. Integrate `d ln(n_b)/d ln(epsilon)=epsilon/(epsilon+P)` with fixed Gauss-Legendre quadrature and define `mu_b=(epsilon+P)/n_b`. Set the auxiliary normalized columns from the cutoff state, recomputing `P/n_b`, `mu_b/m_n-1`, and cold `epsilon/(n_b*m_n)-1` while leaving composition indicators fixed and clearly documented as a cutoff-frozen auxiliary representation.

- [x] **Step 4: Verify the focused pass**

Run: `PYTHONPATH=. pytest -q python_frontend/tests/test_surface_continuation.py`
Expected: all cutoff, first-law, C1-family, causality, monotonicity and writer/provenance tests pass.

- [x] **Step 5: Run affected integration checks**

Run: `python -m compileall -q python_frontend/compose_cowling/surface_continuation.py`
Expected: no syntax errors.

- [ ] **Step 6: Commit the independently passing component**

Commit the module and its focused tests together.

### Task 2: Six-method native TOV continuation experiment

**Files:**
- Create: `scripts/run_surface_continuation.py`
- Create: `.github/workflows/surface-continuation.yml`
- Consume: `oscillation_release_v1/data/spectrum_summary_234.csv`
- Consume: `oscillation_release_v1/data/eos_projected/GPPVA_FSU2__BBB_BHF2__*.csv`
- Consume: `python_frontend/compose_cowling/solver.py::solve_v5_tov_star`

**Interfaces:**
- CLI `--output-dir`, `--tail-points`, `--coarse-energy-ratio`, `--fine-energy-ratio`.
- Output `surface_continuation_summary.json`, `surface_continuation_summary.csv`, and generated extended normalized EOS tables.

- [ ] **Step 1: Add a failing CLI/unit integration test for input cohort selection and case design**

Require exactly six methods for `GPPVA_FSU2 × BBB_BHF2`; central pressures and reference mass/radius come from the frozen summary. Cases are original finite-cutoff baseline; tangent power-law at `x_min=1e-8` and `1e-12`; exponential quadratic deformation at `alpha=±0.5`; and polynomial quadratic deformation at `alpha=±0.5`, all shape-deformed cases at `x_min=1e-12`.

- [x] **Step 2: Implement the runner**

Reproduce each finite-cutoff baseline with v5 `cell_solve task=star`. Reject baseline if relative mass or radius differs from the frozen summary by more than `2e-5`. Generate each continuation table, run at the same reference `Pc`, record `M`, `R`, surface pressure, residual, backend build key, per-case relative changes, continuation-parameter differences, and `x_min` refinement differences. Write partial JSON after each attempt so failures preserve casewise evidence.

- [x] **Step 3: Add a dedicated Actions workflow**

Install GNU Fortran, NumPy and pytest; run the focused tests; execute the six-model TOV continuation experiment; upload all summary/table outputs even on failure. Do not alter the existing mass-dependent workflow or its 384-point pilot.

- [ ] **Step 4: Verify the exact-head workflow and artifact**

Expected: six baseline replays and all requested continuation integrations either produce an explicit casewise result or a preserved failure; no implicit extrapolation and no modification to source knots.

### Task 3: Research DAG, result audit and coding-agent handoff

**Files:**
- Modify: `docs/MASS_SEQUENCE_BRANCH_STABILITY_AUDIT_20261109.md`
- Modify: `docs/MASS_DEPENDENT_COWLING_DAG_20261109.md`
- Modify: `docs/MASS_DEPENDENT_COWLING_RESEARCH_PLAN_20261109.md`
- Modify: `docs/LOCAL_AGENT_RESUME_MASS_DEPENDENT_COWLING_20261109.md`

- [ ] **Step 1: Record the exact source SHA, cutoff reconstruction and analytic closure**
- [ ] **Step 2: Record all six-model TOV mass/radius shifts and cutoff/shape sensitivity with artifact digest**
- [ ] **Step 3: Keep the selected surface model unresolved unless the cutoff and family sensitivities support a predeclared tolerance**
- [ ] **Step 4: Keep radial `omega_0^2`, `Gamma_1`, and f/p1 work blocked pending the next explicit physics decision**
- [ ] **Step 5: Push a docs-only closeout commit and mirror the checkpoint to Dropbox Paper**

---

## Explicit unresolved decisions

- No literature-sourced outer-crust EOS is present in the projected table. The local tangent law and the two C1 curvature families are controlled mathematical continuations, not physical crust reconstructions.
- The acceptance tolerance for the `x_min=1e-8` versus `1e-12` surface-limit shift and for between-family TOV response is not predeclared. The first run must report these sensitivities without declaring a physical model selection from arbitrary tolerance choices.
- Frozen-composition `Gamma_1` remains unavailable from the projected table; an eventual `Gamma_1=Gamma_eq` radial calculation may only be described as barotropic-equilibrium radial stability.


## Execution status at push preparation (2026-11-09)

- Component and runner files are implemented locally; the focused test command completed with `8 passed`. `compileall` completed without syntax errors.
- The observed red failure was a missing continuation-module import before the implementation was added; the same focused suite is green after implementation.
- The six-method native TOV experiment is configured for 42 runs: six finite-cutoff baseline reproductions and 36 positive-tail continuations. It has not yet executed in this workspace; the dedicated exact-head GitHub Actions run is the integration evidence gate.
- `gamma_s=(epsilon_cut/P_cut) c_s^2` is the power-law exponent. `Gamma_eq=((epsilon_cut+P_cut)/P_cut)c_s^2=gamma_s+c_s^2` is the equilibrium adiabatic index. They are distinct and both are recorded.
- The tail terminates at a positive `epsilon/epsilon_cut` ratio (1e-8 or 1e-12), because the frozen native table reader rejects literal zero primary fields. The exact P=0 limit is asymptotic; the lower-cutoff comparison tests approach to that limit.
- No continuation has yet been selected as a physical outer-crust EOS. No radial `omega_0^2`, radial stability, or new Cowling f/p1 claim is part of this implementation.
