# Native finite-table retained-source audit

Source: `fortran_eos/tests/source_prefix_audit.f90`.
Tests: `python_frontend/tests/test_source_prefix_audit_native.py`.

Compile without the production dispatcher:

```sh
gfortran -O2 -ffree-line-length-none -fcheck=all \
  -ffpe-trap=invalid,zero,overflow \
  fortran_eos/tests/source_prefix_audit.f90 -o source_prefix_audit
./source_prefix_audit input.tsv certificate.json
```

The first nonblank, noncomment record is `mn <reference rest energy in MeV>`.
Remaining records have exactly four whitespace-separated fields:
`nb Q1 Q7 Iphase`. Blank lines and lines beginning with `#` are allowed.
The frontend copies these fields from the actual retained CompOSE source rows;
it does not perform the physical conversion. The standalone Fortran program
calculates

\[
P_i=n_i Q1_i,\qquad
\varepsilon_i=n_i m_n c^2(1+Q7_i).
\]

All physical arithmetic and interval comparison use `real128`. The source
uncertainty factor is exactly `b=64*epsilon(binary64)`; the halfwidths are
`b*abs(n)`, `b*abs(Q1)`, `b*abs(mn)`, and `b*max(1,abs(Q7))`. The Q7 absolute
floor covers cancellation in the subtract-one CompOSE normalization. This is
an explicitly declared roundoff/serialization allowance, not a new matching
tolerance, arbitrary physical repair, or uncertainty of the original nuclear
model. Decimal input tokens must be finite and within the binary64 range.

For each interval, the stable cold first law implies the necessary bounds

\[
P_1(1/n_1-1/n_2)\le
\varepsilon_2/n_2-\varepsilon_1/n_1
\le P_2(1/n_1-1/n_2).
\]

The auditor retains correlations in the original fields rather than assigning
independent errors to every converted quantity:

- `Delta(epsilon/n)=mn*(Q7_2-Q7_1)`, with the same uncertain reference mass;
- the lower bound is `Q1_1*(1-n_1/n_2)`;
- the upper bound is `Q1_2*(n_2/n_1-1)`.

Density uncertainty intervals must remain strictly ordered. Corner extrema of
these expressions enclose the source boxes; an additional `128*epsilon(real128)`
evaluation allowance is reported. A violation is certified only when the
entire energy-difference interval lies outside the necessary pressure-integral
interval. Pressure decrease is separately refuted only when `P2_max<P1_min`.
Nominal violations that overlap the declared roundoff budget are counted
separately and do not become a false positive rejection.

Different phase IDs with overlapping P and baryon-chemical-potential intervals
are identified as coexistence-compatible. Exact one-sided common-P, common-mu
endpoints are accepted. This check does not infer an unrepresented mixed phase.

## JSON and exit contract

The stable fields are:

- `status`: `FINITE_TABLE_NECESSARY_BOUNDS_SATISFIED`,
  `REFUTED_BY_FINITE_TABLE`, or `INVALID_INPUT`;
- `rows`, `pressure_decrease_count`, `firstlaw_violation_count`;
- `first_bad_interval_one_based` (zero when none);
- `nominal_pressure_decrease_count`, `nominal_firstlaw_violation_count`;
- `coexistence_compatible_interval_count`;
- `intervals`, containing original row values, converted epsilon/P, source
  halfwidths, nominal and conservative bounds, and the individual decisions.

Parsed inputs return exit0, including scientifically refuted input. Malformed
input returns exit2 with an `INVALID_INPUT` JSON certificate. Missing or
unwritable output is an execution failure, never a scientific certificate.

Passing means only that these finite-table necessary conditions did not refute
the supplied rows. It does not establish a continuous first law, stability
between knots, original-paper fidelity, or a global energy minimum. The
generator's analytic/branch derivative diagnostics remain separate evidence.
The parent must use exactly the retained source prefix when qualifying a
matched result; a failure in a row that was discarded by the match is not a
certificate about that match.

## Verified examples

The standalone native suite has **20 passing tests**, including an exactly
integrable quadratic energy, explicit first-law and pressure counterexamples,
coexistence with and without binary64 serialization roundoff, mass correlation,
invalid/null/overflow/nonfinite input, unordered density, malformed records,
and hidden trailing data beyond the bounded line buffer.

An integration smoke copied n,Q1,Q7 and phase IDs from real serialized BBP
CompOSE artifacts produced by the final native generator tests:

| Full source | Rows | Pressure-decrease refutations | First-law refutations | First bad interval |
|---|---:|---:|---:|---:|
| Original published closure |257|0|221|1|
| Raw full-derivative variational path |257|1|1|4|
| Variational Maxwell completion |253|0|0|0|

The Maxwell source has one coexistence-compatible interval. These are full-source
smoke results; the parent separately executes the same auditor on the actual
retained prefixes in candidate matched artifacts. No source name is whitelisted.
