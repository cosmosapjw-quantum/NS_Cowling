# COWLING-EOS-SOURCE-CALIBRATION-04

Date: 2026-11-19 KST. Branch
`research/cowling-eos-source-calibration-20261119`.
Authoritative parent:
`0d802d09fff57df1cc5557916be824cbc457d2b7`, the
`COWLING-MAXWELL-NATIVE-03` scientific and coding record.

This report records an executed, source-locked and **bounded scientific
calibration** of original cold GPPVA neutron-star EOS data, phase-local
first-law closures and Fortran DVODE sharp-interface Cowling f/p1 spectra.
Do not infer publication acceptance, microscopic phase minimization,
full metric-QNM spectroscopy or global stellar radial stability.

## Direct results

1. **Source identity recovered:** All three original CompOSE unified
   FSU2, TM1e and DDME2 source ZIPs were independently downloaded from
   the official Zenodo records and are **byte-for-byte identical** to the
   previously downloaded CompOSE ZIP files. The separately published
   CompOSE `eos.zip_checksum.txt` SHA-256 nevertheless differs for
   all three, so preserve that provider-index inconsistency rather
   than overwriting evidence. The Zenodo records' MD5 and file sizes
   match the actual original archives:
   - [FSU2 record 14808903](https://zenodo.org/records/14808903),
     MD5 `0ec951fe76c0b637e6cbeeb209237fdd`;
   - [TM1e record 14808958](https://zenodo.org/records/14808958),
     MD5 `e73b477f731f1be05c3ae0e1cb2ad663`;
   - [DDME2 record 14508967](https://zenodo.org/records/14508967),
     MD5 `3ee738c7dea32222002b34b0e4a5df16`.
   This **settles the archive byte identity**, not whether the
   supplied nuclear EOS satisfies every thermodynamic or physical
   assumption used by our solver.

2. **Same named RMF is not automatically source identical:** The
   original provider's homogeneous phase-6 rows were compared
   at EXACTLY their tabulated baryon densities with the local
   functional solvers, without any regridding. The maximum
   `abs(P_RMF/P_provider-1)` in `0.08 <= n_B <= 0.45 fm^-3`:
   - FSU2: `0.0422125627754` (4.2213%);
   - TM1e: `0.1551819453627` (15.5182%);
   - DDME2: `0.0000515592193` (0.00516%).
   Pressure mismatches materially affect phase coexistence. DDME2
   uses original native Fortran density-dependent RMF evaluator;
   the other two use the original stored RMF parameter solver.
   **Classification:** FSU2/TM1e are *same named parameter family,
   not source-level numerically equivalent*; DDME2 shows high
   numerical agreement but is not byte-identical either.
   See `source_audit/RESULT.json` for every evaluated density,
   `P,epsilon,mu_B` and beta/Euler residuals.

3. **True phase-local first-law repair:** The cold thermodynamic
   identity is
   ```text
     F(n)=epsilon(n)/n,   dF/dn = P(n)/n^2,
     mu_B=F+P/n,         c_s^2/c^2 = dP/d(epsilon).
   ```
   Perform `P(n)`-primary reconstruction separately on the
   original `Iphase=0` outer, `Iphase=7` inner and
   `Iphase=6` homogeneous core. Preserve original source `n,P`
   exactly at knots and anchor `F` to the ORIGINAL first point
   of each segment. No fitted `mu_B` offset. **Do not integrate
   across a provider's phase boundary**: the previous cross-seam
   pressure-primary interpolation produced a different effective
   physical barotrope.

4. **Sharp Maxwell vs CSS remain separate physics.** With
   phase-local repair of the *official unified* FSU2, TM1e and
   DDME2 source, native CSS finite matching passes on all three;
   original separated phase data do not admit a Maxwell crossing
   on the verified overlapping support. The CSS transition data:
   | Official model | Pt [MeV/fm³] | n_low [fm^-3] | n_high [fm^-3] |
   |:--|--:|--:|--:|
   | FSU2 | 0.2546799995 | 0.054000001 | 0.064722370 |
   | TM1e | 0.3455899893 | 0.089000002 | 0.124896752 |
   | DDME2 | 0.4067599943 | 0.071999997 | 0.090984385 |

   This is **a declared phenomenological CSS jump**, not proof
   that the original provider full unified source contains a
   thermodynamic first-order Maxwell equilibrium. The
   earlier *source-original crust + newly evaluated FSU2
   homogeneous core* model **does** have a positive-density-jump
   finite Maxwell crossing, but section 2 forbids identifying it
   with the original provider's homogeneous phase-6 physics.

5. **At equal 1.4 M_sun**, the corrected `Iphase=7` and
   `Iphase=6` separate first-law models have measured
   differences from the earlier cross-seam first-law closure:
   | Model | R1.4 [km] phase-local | Delta R [km] | Delta slow p1 [Hz] | Four strict nodes? |
   |:--|--:|--:|--:|:--|
   | FSU2 | 13.8491568 | -0.0214378 | +47.3837 | yes |
   | TM1e | 13.1267068 | -0.0374269 | conditional slow W1/V2 | no |
   | DDME2 | 13.2308341 | -0.0311919 | +220.5755 | yes |

   These changes are MODEL-DEPENDENT PHYSICAL EFFECTS from the
   phase-local `epsilon,mu` construction, not precision errors
   of the Fortran integrator. The two successful strict-nodal
   cases also have native TOV parity with Python in the input
   campaign. TM1e has a real-looking and numerically stable
   but **strict-node-inconclusive** overtone on part of the
   central-density range.

6. **Error decomposition identifies source coefficient sampling.**
   At `n_c=0.30 fm^-3`, source-original FSU2 crust + generated
   FSU2-core Maxwell example, slow `p1`:
   | Method | Frequency [Hz] | Absolute error vs independent Python DOP853 |
   |:--|--:|--:|
   | Independent Python DOP853 Cowling | 4897.917834808 | reference |
   | Python TOV + Fortran Cowling, 4097 coefficients | 4898.004398533 | 0.08656373 Hz |
   | Python TOV + Fortran Cowling, 32769 coefficients | 4897.919675627 | 0.00184082 Hz |

   At fixed 32769 Cowling coefficient samples, switching among
   4097 to 32769 native-TOV background samples moves `p1`
   by typically `1e-4 - 3.4e-4 Hz` relative to the Python
   background reference. This tests **source coefficient
   interpolation** independently from the ODE integrator and
   achieves an observed ~47.0-fold reduction in absolute error
   compared to 4097 coefficients. It does NOT establish
   a universal convergence order or error estimate.

7. **Unnecessary breakpoint forcing rejected.** Added an
   optional `NativeSharpCowlingDVODE(..., knot_waypoints=True)`
   numerical method that forces DVODE output at every original
   one-sided `h` breakpoint, never at a density discontinuity
   and without modifying EOS values. Four source models × three
   resolutions × four modes = **48** measured frequency pairs.
   All had the same eigenfrequencies (maximum observed absolute
   change `0.0 Hz` at the stored precision), and forced
   waypoints increased median warm run time by ~8.4%.
   Therefore default remains **OFF**, and the method is not
   asserted to improve the solution.

8. **Physical mode-family continuity, with an explicit
   counterexample to a rigid node label.** Three FSU2
   models × both conversions × five central densities
   `(0.28,0.30,0.33,0.36,0.39) fm^-3` return a total of
   **48** adjacent-configuration f/p1 mode associations
   (24 star-to-star comparisons, two modes each).
   The phase-wise sign-invariant `(W,rV)` shape
   overlap minimum is **0.9909994**; these are not
   merely frequency-ordered assignments.

   For the separate **TM1e phase-local CSS** source, the
   slow `W1` eigenfrequency varies smoothly as
   `n_c=0.34,0.35,0.3668085,0.38 fm^-3`:
   ```text
      6553.1706  ->  6623.3697  ->  6692.2004  ->  6719.2903 Hz
      (W,V): 1,1       1,2           1,2           1,1
   ```
   Each node count is reproducible at 160,450,900,1800
   output points. At `n_c=0.3668085`, the interface
   radius is `11.81918083 km`, with the extra two
   `V` zero crossings at `11.68210 km` inside the
   core and `11.82232 km` on the exterior crust side.
   The overlap across the W1/V1 -> W1/V2 change at
   `0.34 -> 0.35` is **0.9990642**; for the return
   `0.3668 -> 0.38`, **0.9974330**.

   Default mode labeling remains **strict** and reports
   `UNRESOLVED` at both V-node changes. A new
   *opt-in* `v_node_policy="shape_continuity"` allows
   only `|Delta V_nodes|=1`, unchanged `W_nodes`
   and **overlap >=0.995**, with explicit
   `CONTINUED_WITH_V_NODE_CHANGE` status. This
   is a *conditional numerical continuation*,
   **not** a proof of general Sturm-Liouville mode identity
   or a claim about physical damping/elasticity.

## Source authority, repair limitations and claim gate

- The three official GPPVA archive bytes are independently
  authenticated against their published Zenodo records.
  Their separate CompOSE `eos.zip_checksum.txt` values
  remain inconsistent with the actual ZIP bytes; a provider
  metadata repair is still required to close that
  administrative mismatch.
- BBP still uses the secondary, approximately 3-significant-figure
  Diener 2008 Table 10.4. Original BBP 1971 physical theory:
  Baym, Bethe & Pethick, *Nucl. Phys. A* **175**, 225 (1971),
  DOI `10.1016/0375-9474(71)90281-8`. Do not assign a
  fabricated original BBP digital table to this stage. The
  03 BBP+BPS physics and accepted CSS models remain frozen;
  phase-local first-law does not silently make them microscopic.
- Official GPPVA DDME2 full EOS contains a high-density
  nonmonotonic tail (last two original density rows);
  the separate first-413-row numerical projection is
  source-hashed and was never merged back into the original.
  In addition, auditing the whole retained DDME2
  homogeneous phase over its full density support finds a
  supra-causal interval. The bounded `n_c ~0.3-0.35 fm^-3`
  stellar runs do not extrapolate to that high-density region.
- A same-named local RMF core does not itself prove source
  equivalence for FSU2/TM1e; the measured exact-density
  differences in section 2 restrict claims accordingly.
- Both the sharp slow and rapid fluid conversion conditions
  are those documented in 03; no full GR, superfluid,
  elasticity, finite conversion-time or microscopic
  phase competition problem is solved here.

**Status:** `COWLING-EOS-SOURCE-CALIBRATION-04:
SCOPED_NUMERICAL_AND_SOURCE_IDENTITY_VERIFIED` for the
documented 1.4-Msun bounded stars, original source bytes,
phase-local source reconstruction, node and coefficient
convergence experiments. `TM1e strict slow p1 mode
identity: UNRESOLVED`; `TM1e W1 branch numerical
continuity: CONDITIONAL`; `source-original FSU2/TM1e
homogeneous core functional identity: NOT ESTABLISHED`;
`original BBP digital numerical source: UNRESOLVED`.
No global physical EOS proof or publication claim.

## Reproduce and inspect

Clone the exact branch, then run on the authorized
machine with GNU Fortran 13 and Python 3.12:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend
# Source archives and secondary BBP PDF are not redistributed on Git:
python3 scripts/fetch_cowling_m03_sources_20261119.py
# Independently fetch official Zenodo archives by record and verify bytes,
# per-record MD5 and physical file size, without overwriting source bytes.
python3 scripts/verify_official_zenodo_archive_identity_20261119.py
python3 scripts/audit_native_rmf_source_calibration_20261119.py
python3 scripts/audit_phase_local_thermodynamic_repairs_20261119.py
python3 scripts/compare_phase_local_fixed_mass_20261119.py
python3 scripts/decompose_sharp_tov_cowling_error_20261119.py
python3 scripts/benchmark_sharp_waypoint_calibration_20261119.py
python3 scripts/continue_sharp_mode_families_20261119.py
python3 scripts/trace_tm1e_mode_identity_20261119.py
python3 scripts/localize_tm1e_v_nodes_20261119.py
python3 scripts/summarize_cowling_calibration_20261119.py
# Focused stage-04 and relevant native TOV/Cowling regression (see test log):
pytest -q python_frontend/tests/test_cowling_mode_continuation_20261119.py \
  python_frontend/tests/test_cowling_sharp_knot_waypoints_20261119.py \
  python_frontend/tests/test_cowling_source_calibration_20261119.py \
  python_frontend/tests/test_cowling_calibration_evidence_20261119.py
```

### Fresh implementation verification and exact source binding

The final scoped acceptance was actually executed with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend \
  python3 scripts/verify_calibration_acceptance_20261119.py
```

Result: **122 tests passed, 3 skipped, exit code 0**. GNU Fortran
13 compiled the actual TOV and Cowling native source modules under
`-O0 -Wall -Wextra -fcheck=all -ffree-line-length-none` (exit 0;
unused DVODE dummy-argument warnings only). `compileall` and `git diff
--check` passed. This command binds the observed test outputs to a SHA-256
snapshot of all 16 scoped test files and the six principal implementation
modules, before and after the run. Source changes after the receipt invalidate
its equivalence, even if the recorded log remains positive.

Evidence: `data/cowling_eos_source_calibration_20261119/verification/
VERIFICATION_RECEIPT.json` and `verification/TEST_RUN.log`, with all
original source bytes stored locally, not republished on GitHub. The
machine-readable overall outcome also verifies matching source hashes and
reports `SCOPED_CALIBRATION_SOURCE_AND_IMPLEMENTATION_VERIFIED`.

The Zenodo independent binary readback was executed as
part of this stage and the resulting hash receipt is
`data/cowling_eos_source_calibration_20261119/zenodo/ZENODO_INDEPENDENT_ARCHIVE_READBACK.json`.
Run `pytest` with the stage-04 tests and the complete
03/native regression suite before reproducing any
implementation-verified claim.

Single summary/SSOT:
`data/cowling_eos_source_calibration_20261119/CALIBRATION_OUTCOME.json`.
Detailed evidence, per-density source rows and case timings
are in the section-specific `RESULT.json` records.
Source raw ZIP/PDF bytes are locally preserved but
not automatically re-published in Git.

**Next bounded science task if selected:** identify the
FSU2/TM1e provider-core RMF parameter and equilibrium
convention discrepancy at exact densities, and derive a
source-equivalent homogeneous core, *before* reinterpreting
same-functional Maxwell existence as source-original
phase coexistence. The primary numerical bottleneck
for f/p1 accuracy is now known: Cowling source
coefficient representation rather than native TOV.
