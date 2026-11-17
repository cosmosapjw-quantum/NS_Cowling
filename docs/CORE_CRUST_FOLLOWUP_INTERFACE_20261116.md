# Source-interface intake and bounded follow-up, 2026-11-16

Base `4df4d36a18fb7e3ba7b1cfe2d7bbcc659a56ec07`; isolated detached worktree. Project status remains **PARTIALLY_COMPLETE**. This node supplies and checks only the spherical interface functions; it does not wire a new CLD solver or claim source calibration. The v6 proxy and its receipts are untouched.

## Source and calibration boundary

[Scurto et al. 2024, Eqs. 24–25](https://arxiv.org/html/2402.15548v1) specify spherical Phi and sigma, and state that sigma parameters are fitted to nuclear masses for each bulk model. Their Eq. 31 states 2,408 masses and the fit weighting. A copy of the source HTML was obtained with SHA-256 `7f1d4227db6df5f462d13db1f514303eab79f5724ea04bc4c1aa7b4e026bc69b`. [AME2016](https://www-nds.iaea.org/amdc/ame2016/mass16.txt) raw file was acquired, SHA-256 `2167f57a2a98331e4649b2dd2b658a9006ed4fba1975729ebfe52a42b4b9218a` (418937 bytes, 3475 lines); local copy `build/source_intake/ame2016/mass16.txt` is not committed. Its identity does not establish which 2,408 rows were fitted, the precise mass-fit implementation, or fitted model coefficients. Source bulk-model identity must also be resolved: the local SPG surrogates and the source's GDFM model are not automatically interchangeable. No mass fit or calibrated EOS was produced.

## Implementation and validation

`python_frontend/compose_eos/scurto_interface.py` implements Eqs. 24–25 and the analytic composition derivative independently of the proxy. Inputs are explicit; no fitted values are supplied. The module is not connected to CLD reconstruction or output emission. Geometry is factored near f=1; the surface denominator is rearranged and scaled to reduce cancellation and overflow. The implementation fails closed outside conservative numerical bounds 1e-8 <= y <= 1-1e-8, sigma_0 <= 1e6 and b <= 1e100; these are computational support limits, not physical constraints. These checks are float64 checks, not a proof over every representable value.

Independent review found the first implementation returned zero for `f=nextafter(1,0)` and NaN for the derivative at `(y,sigma_0,b)=(0.1,1.3,1e300)`; it also identified cancellation near `b=-16`. Those failures were reproduced in `logs/ns_scurto_review_red.log`, and decimal-oracle and endpoint tests passed after repair. A second review found silent underflow/overflow for further extreme finite inputs; four refusal tests failed before the numerical support bounds were imposed and passed afterward. The original missing-module test failure is also retained. Logs are in `data/receipts/core_crust_followup_20261116/logs/`.

Commands from repository root, base before interface edits:

```sh
bash fortran_eos/build.sh build/compose_eos_fortran
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 -m pytest -q python_frontend/tests pipeline/tests
PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg python3 -m pytest -q eos_generator/tests
```

Build succeeded; frontend/pipeline 460 passed, 9 skipped; generator 142 passed with an existing Matplotlib Axes3D warning. After the interface change: `PYTHONPATH=python_frontend python3 -m pytest -q python_frontend/tests/test_scurto_interface_functions.py` gave 17 passed. Python 3.12.3, NumPy 2.4.2, SciPy 1.17.0, gfortran 13.3.0. Log SHA-256 values are in `SHA256SUMS.txt` adjacent to the logs. The full suites are baseline tests, not a post-change combined-suite result. No census or stellar calculation was rerun.

## Remaining obligations

HOLD source-calibrated CLD pending exact model bulk inputs, selected mass dataset and fit protocol/coefficients; wire and validate the complete energy and all variational derivatives only once identified. Multistart and density convergence, source-backed outer crust and physical surface, separately established mode identity and radial stability, DD2/RG conventions, composition-aware Maxwell and EOS390 DAG campaign remain open. Do not promote proxy outputs or infer spectral completeness or global minimum from this node.

Independent read-only review returned `PASS_SCOPED` for the equation-only helpers after these repairs, reporting a 630-case Decimal comparison with maximum relative error 1.02e-15; the detailed sweep output was not retained, so the reproducible committed evidence is the targeted test log rather than that sweep. The review does not admit a calibrated CLD or stellar result.

## Post-change verification and test repair

The interface was published as 594552d7a597dd150be14b9a860a8f5eaa571d7f.
On that clean commit, a fresh native build and full local suites completed:
**477 frontend/pipeline passed, 9 skipped (38.71 s); 142 generator passed (14.53 s)**,
with the same pre-existing Matplotlib warning. Thus the earlier paragraph describes
the original baseline runs, not the latest available verification.

[CI #614](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/43786808560)
also succeeded on merge 3bf8b045271e752eb4cb705a0111a0ca97665a16
(head 594552d, base b12ba451), with 477/9 and 142 test results.
That workflow ran the existing proxy census; no source-calibrated census or new stellar result follows.

Independent review found a test weakness: the near-f=1 pytest.approx assertion's
default absolute tolerance admitted a value 100 times its expected value. The assertion
now sets abs=0.0; all 17 targeted tests pass after this repair. No helper formula changed.

A new 400-digit Decimal check retained its script and all **1,722 scalar outputs**:
zero failures, maximum relative error for normal outputs **1.0180684467550535e-15**,
maximum subnormal/zero-reference error **1 ULP**. It includes the support boundaries
and extremely small sigma_0. These finite samples do not prove all floating-point cases.
Tiny exact outputs may round to zero. This is separate evidence from the earlier
unretained 630-case sweep.

See data/receipts/core_crust_followup_20261116/verification/RECEIPT.json,
the adjacent raw logs, compressed Decimal results, replay script and checksum manifest.
The nine original follow-up logs verified against their hashes. The v6 proxy and prior
audit receipts remain byte-identical to base 4df4d36. AME2016 and source-HTML
identities remain the supplied intake provenance; they were not reacquired in this check.
Equation-only verdict **PASS_SCOPED**; scientific project **PARTIALLY_COMPLETE**.
