# Finite cold-EOS causal endpoint audit

The native read-only checker takes exported one-dimensional CSV files with the header `nb_fm3,epsilon_MeV_fm3,pressure_MeV_fm3,Iphase,region`. Number density is in fm⁻³, energy density (including rest energy) and pressure are in MeV fm⁻³. The assumed EOS is cold, has one independent baryon density and obeys the first law along the selected composition path. At fixed charge fraction the relevant derivative is the compression enthalpy along that fixed-composition family.

Let h=(ε+P)/n>0. The cold first law gives dε=h dn, dP=n dh, and the sound-speed ratio is dP/dε=d ln h/d ln n. Hence, on a thermodynamically stable causal branch, h and n/h are nondecreasing. At a first-order jump h and P can remain constant while n increases. For h₂≥h₁ define A=(h₂−h₁)(h₂+h₁)/2. Integrating n dh then requires

    n₁ A/h₁ ≤ P₂−P₁ ≤ n₂ A/h₂.

The inequalities are necessary endpoint conditions. They allow hidden jumps and do not certify the chosen interpolant, microscopic beta equilibrium, extensions of the source to other densities, or stellar perturbation physics. This is the causal pressure-area constraint applied to adjacent knots. For the thermodynamic construction of crust-core matching see Fortin et al., Phys. Rev. C 94, 035804 (2016), §III, DOI 10.1103/PhysRevC.94.035804. For the general causal endpoint pressure bounds see Komoltsev and Kurkela, Phys. Rev. Lett. 128, 202701 (2022), DOI 10.1103/PhysRevLett.128.202701.

## Native entry and reproduction

Compile with a compiler providing real128 (the recorded run used GNU Fortran). The shell entry compiles in an isolated temporary directory and accepts several input tables in one invocation. A positive violation count is a physical result reported per table. A malformed input produces an unsuccessful process exit. `--intervals` requires a new output path and writes both bounding states, pressures, allowed area, guaranteed excess, region and failure category. `--relative-uncertainty r` treats each of n, ε and P at every endpoint as independently varying by ±r. The minimum numerical allowance is 128 binary64 epsilons.

~~~sh
bash fortran_eos/run-causal-audit.sh --self-test
# Set PHYSICAL_CSV_DIR to the extracted PUBLIC_MATCHED_ROWS.zip directory.
bash fortran_eos/run-causal-audit.sh --intervals failures.csv "$PHYSICAL_CSV_DIR"/*/*.csv > summary.csv
bash fortran_eos/run-causal-audit.sh --relative-uncertainty 1e-5 "$PHYSICAL_CSV_DIR"/*/*.csv > sensitivity_1e-5.csv
~~~

The native audit streams tables, retaining only the preceding state. A pressure-area violation means every point consistent with the supplied endpoint boxes violates the required interval. At positive uncertainty the checker takes conservative outer bounds for h and for each allowed area, and ignores shared-variable correlations. Failure is conclusive at that error budget. A non-exclusion at enlarged uncertainty does not assert that one set of model parameters or a common curve can realize all of the individually allowed endpoint pairs.

## Results

The full 787-table public archive from the preceding commit provided 567,013 physical rows. The nominal run rejected 177 tables in 79,465 intervals. All failures are inside the inherited core segments. The 215,870 crust pairs, 62,305 bridge-interior pairs and 1,574 seam pairs had zero violations. The 79,465 repeated intervals reduce to 2,361 distinct (core identity, left density, right density) keys. All 566,226 pairs pass the ordinary adjacent sound-speed secant test.

`PUBLIC_ENDPOINT_AUDIT_787.csv` joins each native outcome to the original census key and source identity. `CORE_UNIQUE_INTERVALS.csv` reports one representative, the one with greatest relative excess, for each violating key. It gives the endpoints, their enthalpies and the bound from the native checker. `FULL_CRUST_110.csv` gives the first bracket of offending pressures among the 110 public full-crust matches. A zero left-row index denotes no violation. `SUMMARY.json` gives uncertainty sensitivities and the original ZIP identity. At the default error allowance 85 of 110 public full-crust candidates pass this additional endpoint check. In the separate protected set, all 35 tables pass, including all five full-crust cases. Thus 90 of 115 original full-crust candidates remain unexcluded across the two sets. No private source states are included here.

`NV_NATIVE_CHECKS.csv` checks new Fortran NV output in the unchanged benchmark and at 1, 2 and 4 OpenMP threads. The printed branch has 49 rows, the adjusted branch 49 and the uniform beta branch 257. All 352 adjacent intervals of each triplet satisfy the endpoint constraints. The physical fixed-composition inputs retain their beta-equilibrium qualification from the generator.

The relative source-column perturbations in the sensitivity survey are hypothetical budgets. The nominal violations measure inconsistency of the literal published numeric representation at the consumed precision; they do not identify the historical mechanism or a unique physical correction. For the stellar interpretation, a candidate whose first failing bracket lies above the chosen central pressure does not visit that bracket in that stellar model.
