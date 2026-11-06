# Adversarial audit and bounded repair report

## Audit contract

The audit sequence was executed as an evidence workflow: self-discover, step-back, CRAG/source retrieval, metacognitive self-ask, CoVe/claim verification, adversarial self-ask, chained derivation and pre-deployment review. These labels identify review functions, not private model reasoning. The review was bounded to five repair cycles. No native parallel subagent facility was invoked, and no simulated delegation is claimed.

## Principal correction

The initial paper wording could be read as saying that weighted isotonic regression materially repaired the source q knots. That is false for this cohort. The source q=n_B/mu_B sequences are already nondecreasing for all 234 tables, and the stored PAV receipts show zero pooled knots. The actual modification is the between-knot closure: q is declared linear, pressure is reintegrated exactly, and epsilon is recalculated. This changes P and epsilon knot values while leaving n_B unchanged to roundoff. The paper, manual, README and code API now call this a q-linear solver-adapter EOS, not a source EOS repair.

## Cycle 1 — Source discovery and provenance

Finding: the earlier preparation script imported compose_eos through a hard-coded absolute path and read matched table trees outside the delivered oscillation release. This was not portable and obscured which EOS generator had authority.

Repair: the integrated release vendors the complete compose_eos 3.3.0 release and the exact seminar source package. The new intake command has explicit matched-root, classifier, generator-package and output-directory arguments. It stores hashes for each source EOS file and matching metadata. The previous normalized and adapter tables are retained without overwrite.

Verification: source package, generator package, prior 390-way matrix provenance and all 234 accepted mode tables are present in the integrated release. Existing spectra were reused, not recalculated.

## Cycle 2 — Physical and thermodynamic step-back

Finding: first-law consistency at EOS knots and cellwise interpolation admissibility had been conflated. The Fortran solver uses cubic Hermite P(mu_B), not an abstract continuous thermodynamic function returned by the EOS generator.

Repair: the documentation now separates the source master function, the raw table knots, the native Hermite interpolator and the q-linear adapter. The q-linear closure is derived explicitly. It is not described as a physical phase transition or as a microscopic EOS reconstruction.

Verification: for every source table, epsilon+P=mu_B n_B holds at the knots within the existing 2e-12 intake gate. For the adapter, epsilon+P=mu_B n_B is rechecked after construction. The cell formula c_s^2=q/(q+mu_B q') is checked for positive q and nonnegative q' and is exactly representable by the native cubic Hermite pressure polynomial.

## Cycle 3 — Implementation, API and edge cases

Finding: the original helper functions were not a stable public pipeline API and assumed one fixed file layout. The original PAV operation was not distinguished from the subsequent thermodynamic re-integration.

Repair: the integrated pipeline provides importable adapter functions, a CLI, a generator-backed intake step, receipts, explicit PAV policy, finite/positive/strict-monotonicity gates and refusal on invalid source identity or noncausal output. No source table is overwritten.

Verification: synthetic EOS tests cover exact q-linear integration, first-law closure, causal sound speed, PAV no-op and rejection of an invalid first-law input. An actual cohort table is rebuilt and compared against its saved adapter.

## Cycle 4 — Numerical and result consistency

Finding: a successful mode root is not enough; it must also have the expected eigenfunction node identity and pass stellar matching and mass gates. A source representation change must not silently replace old mode inputs.

Repair: the release keeps source CSVs, adapter CSVs and receipts separately. The verifier reads the existing 234-row summary and explicitly checks unique identifiers, 78 EOSs per core family, valid status, native cell audit, positive frequencies and (0,0)/(1,1) node counts.

Verification: existing 234-mode calculations were not rerun. The saved results pass the recorded target-mass, matching, surface and node diagnostics. The adapter comparison records the original/adapter differences. The q PAV stage pools zero knots in this dataset; pressure/energy re-integration is the nontrivial change.

## Cycle 5 — Literature, manuscript, build and release

Finding: the original manuscript cited Finn 1988 around the Cowling derivation without clearly identifying that the code uses the McDermott-Van Horn-Scholl prescription in which all metric perturbations are neglected. A related Finn formulation retains part of the metric perturbation structure.

Repair: McDermott et al. 1983 is added as the primary source for the equations implemented here. Finn 1988 is retained as a related but distinct formulation. Yoshida and Kojima 1997 support the model-accuracy caveat. Ranea-Sandoval et al. 2018 and Tonetto and Lugones 2020 support the discussion of sharp-interface hybrid-star modes and conversion-rate-dependent junction conditions. Montefusco et al. 2025 supports the frozen-composition versus beta-equilibrated distinction.

Verification: the updated REVTeX manuscript and technical documents are built where the local TeX toolchain is available. BibTeX records are checked with the academic reference validator. The integrated package is tested, hashed and archived. The saved f/p1 spectrum is read-only during this cycle.

## Claim gates and residual limitations

Established by derivation: under the cold one-component barotropic first law, q>0 and q'>=0 imply positive compressibility and 0<c_s^2<=1 for the q-linear representation.

Numerically checked: source-knot first-law residuals, q monotonicity, PAV pooled-knot counts, adapter first-law residuals, pressure/energy changes, analytic adapter sound-speed bounds and cubic-Hermite reproduction are checked by the builder or its tests.

Previously computed and reused: the 234 f/p1 spectra, node identities, TOV masses, matching diagnostics and representative convergence table.

Not claimed: full-GR frequencies, damping times, frozen-composition response, exact vacuum-surface frequencies, sharp-interface discontinuity modes, a global error bound from three convergence examples, or exact physical equivalence between source EOS and adapter EOS. The adapter is a declared numerical thermodynamic realization and its pressure/energy changes are model systematics.
