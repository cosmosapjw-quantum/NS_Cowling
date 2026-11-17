# Native SLSQP adapter

The unmodified fixed-form `slsqp_optmz.f` is from SciPy v1.15.3,
commit `e29dcb65a2040f04819b426a04b60d44a8f69c04`.
`PROVENANCE.json` records immutable upstream URLs and SHA-256 identities.
The Algorithm 733 BSD-release permission remains in its original header;
SciPy's BSD license is included as `SCIPY_LICENSE.txt`.

The free-form adapter is
`fortran_eos/src/compose_eos_optimize_native.f90`.
It accepts a native analytic objective/gradient callback, finite box bounds,
and any number of linear equality rows `A*x=b` (including zero rows).
It does not invoke Python or SciPy at runtime.

Compile the fixed-form source with `-std=legacy -ffixed-line-length-none`;
compile the adapter after `compose_eos_types.f90`. The source includes its
own required level-1 BLAS routines. No external library is needed by SLSQP.
Use one object copy when linking with the physics modules.

```fortran
call minimize_slsqp_native(x,lower,upper,aeq,beq,evaluate,ftol,maxiter,result)
```

The callback signature is:

```fortran
subroutine evaluate(x,f,g,ok)
  real(wp), intent(in) :: x(:)
  real(wp), intent(out) :: f,g(:)
  logical, intent(out) :: ok
end subroutine
```

Callbacks must supply both objective and gradient when `ok=.true.`.
A failed callback terminates immediately (mode -101) at the retained trial
point; nonfinite objective/gradient gives mode -102. No finite penalty is
substituted into SLSQP. In these cases `objective_valid=.false.`; the finite
sentinel in `objective` is unavailable data, not a computed energy.
Modes from original SLSQP are preserved, including 9 for iteration limit.
Mode -100 denotes invalid wrapper input, -103 an exhausted reverse-
communication guard. The `maxiter-1` initialization follows the legacy
SciPy adapter.

`success` means SLSQP mode 0 only. Callers must independently verify physical
feasibility, projected KKT conditions, radius constraints and local curvature.
`max_constraint_error` and `max_bound_violation` are diagnostics, not new
scientific tolerances. Initial points are projected into the finite box;
roundoff-sized trial excursions are clipped only for objective evaluation.
Repeated identical callback points are cached within one invocation.

Native tests include an analytically solved equality-constrained quadratic,
an active-bound minimum, Rosenbrock, iteration-limit and callback failures,
an infeasible equality and invalid bounds, and a 129-variable
box/equality-constrained analytic projection. RED/GREEN logs and commands are
under `data/research_native_20261116/optimizer/`.
