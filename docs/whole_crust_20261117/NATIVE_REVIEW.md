# Whole-crust native kernel: independent review and controls

Reviewed source: `fortran_eos/app/whole_crust_kernel.f90`.
Scientific scope: `docs/whole_crust_20261117/PLAN_KO.md`.
Controls: `python_frontend/tests/test_whole_crust_native.py`.

## Formula review

The pressure-primary mode holds the source density and pressure knots fixed,
sets only the first energy-per-baryon anchor, and integrates the cold first law.
For one interval with (P(n)=P_L(n/n_L)^\Gamma\),

\[
F_R-F_L=\int_{n_L}^{n_R}\frac{P(n)}{n^2}\,dn
=\frac{P_L}{n_L}\log(r)\,
\frac{\exp[(\Gamma-1)\log r]-1}{(\Gamma-1)\log r},
\quad r=n_R/n_L.
\]

The exact \(\Gamma=1\) limit is \((P_L/n_L)\log r\).
The source energy values beyond the first row are diagnostic inputs and are not
additional anchors. The returned chemical potential and sound speed obey
\(\mu=F+P/n\), \(c_s^2/c^2=\Gamma P/(\epsilon+P)\).
For a positive power-law pressure with positive initial F, the energy and
enthalpy remain positive. Sound speed is monotone inside an individual interval,
so checking its two one-sided endpoint values covers that interval.

The surface mode uses

\[
P=K n^\gamma,\quad \epsilon=a n+\frac{P}{\gamma-1},\quad
\mu=a+\frac{\gamma P}{(\gamma-1)n},\qquad
a=\frac{\epsilon_a-P_a/(\gamma-1)}{n_a}>0.
\]

It preserves the supplied anchor, has \(\mu\to a\) and \(P,\epsilon\to0\) as
\(n\to0\), and is causal for \(1<\gamma\le2\). It is an explicit surface
prescription, with physical-composition inference outside the kernel's scope.

For the Hermite bridge, the code's normalized polynomial coefficients give
\(n=P'\), \(\epsilon=\mu n-P\). The extrema checks of n, n', and
\(\mu n'-n\) have the correct algebraic coefficients. They certify the whole
bridge between the endpoints; output sampling alone is not used as the proof.
The code returns the original physical endpoint n, P and epsilon exactly rather
than replacing them with a rounded polynomial evaluation.

## Concrete numerical issue identified

An exactly luminal quadratic, \(P=0.005\mu^2\), has
\(n=0.01\mu\), \(\epsilon=P\), and \(\mu n'-n=0\). Decimal endpoints
\((n,P,\epsilon)=(0.1,0.5,0.5)\) and \((0.12,0.72,0.72)\) produce a small
negative causal-margin value when the Hermite coefficients are formed in binary
floating point. The original strict comparison `cqmin < 0` can therefore reject
the analytically causal boundary. A dedicated native regression below records
the outcome on the actual Fortran compiler.

Any numerical bound must be based on **pre-cancellation terms**, such as the
separate magnitudes of \(\mu n'\) and n. A bound proportional only to the already
cancelled polynomial coefficients can be effectively zero at this boundary and
would not address the error. Preserve the signed residual and bound; pressure,
energy or the reported sound speed must not be silently clipped to manufacture
agreement. Production outer/lower pieces are well inside the causal bound, so
this control exercises a numerical edge separately from their physical result.

## Actual control coverage

The test file contains 29 native cases:

- Four pressure integrals, including gamma exactly 1 and 1±1e-12, compared with
  independent SciPy quadrature in log density.
- Exact preservation of n/P knots and invariance to all unanchored source-energy
  entries.
- Three analytic surface polytropes, including independent differentiation of
  epsilon rather than only checking Euler closure.
- An exact, strictly subluminal quadratic P(mu) and the luminal roundoff control.
- Ten rejected density/range/nonconvexity/acausality cases.
- Six NaN/Infinity cases across all three modes.
- Three existing-output preservation cases.

Local syntax compilation passed. The local pytest collection produced **29
skips**, because that environment has neither gfortran nor a configured native
binary. This local result is not evidence of native execution success. The root
production execution supplies `WHOLE_CRUST_KERNEL_EXE` or compiles the standalone
source and preserves the actual native red/green log separately.

## Bounded reuse of existing upper EOS

The original SPG_M1–SLy4 smooth table and its generation/matching receipts are
unchanged. The independent source audit checked that its all-knot P(mu_E)
Hermite representation is mechanically stable and causal on every smooth
interval. Its upper domain remains n=0.5 fm^-3 in A. This supports reuse for A's
whole-crust assembly without rerunning the unrelated 1,301-case historical
suite. The broader source support and stellar-consumer suitability belong to
the next B/D tasks. M2/M4's representation issue and the original public outer
failures remain in the source audit; they are not reclassified by M1's success.
