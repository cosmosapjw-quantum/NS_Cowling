# Independent composition-aware core generation

The new `scripts/generate_core_eos.py` entry point generates the five original
core models on a caller-defined density grid. It reads no CompOSE target
thermodynamics. The BBB route reads the rounded SNM/PNM energy samples published
in the original paper; it reconstructs that energy representation rather than
re-solving the Brueckner many-body kernel.

## Definitions and execution

All states are cold, locally charge neutral and beta equilibrated, except that
APR explicitly permits a globally neutral two-charge bulk Gibbs coexistence.
The species are npemu for APR, BBB2, KDE0v1 and SLy4; the declared GM1 default is
npe. Energies are MeV, densities fm^-3, pressure/energy density MeV fm^-3.
For equilibrated matter, epsilon + P = n_B mu_B and d epsilon/dn_B = mu_B.
Consequently c_s^2/c^2 = (dP/dn_B)/mu_B = n_B(dmu_B/dn_B)/mu_B.

The native export stores the actual composition, presence masks and APR phase
identities. Q3=mu_B/m_n-1, Q4=-mu_e/m_n, and Q6=Q7=epsilon/(n_B m_n)-1.
The header masses follow each declared model, not an arbitrary common file value.
Tests and the independent review check these conventions and the first law.

```sh
bash fortran_eos/build.sh build/compose_eos_fortran
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. \
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
python3 scripts/generate_core_eos.py build/my_core_generation
```

Use a fresh output directory. `--models`, `--points`, `--n-min` and
`--n-max` select the actual calculation. An APR subwindow does not acquire
coexistence boundaries outside that window. A narrow high-density window
exposed an audit-sampling bug; its red/green evidence and CLI regression remain.
The derivative audit records its actual finite-difference steps after clipping
them to the requested window.

## Executed evidence

The committed campaign has 807 generated rows: APR163 and four other models161.
All five tables passed native write/read, species-presence and phase checks.
Charge/baryon closure is below7.51e-14. The largest independently differenced
first-law discrepancy is5.184e-8MeV. Seven focused tests pass. A separate
72-state off-equilibrium derivative review is in the theory_review evidence.

| Model | Executed n_B interval, fm^-3 | Representation and interpretation |
|---|---|---|
|APR_APR|0.085–1.34|Published21-parameter functional and explicit bulk Gibbs prescription; original-author code cross-check is retained from the prior loop.|
|BBB_BHF2|0.08–1.5|Published Paris+Urbana SNM/PNM energy samples, parabolic asymmetry, PCHIP. Above0.8fm^-3 the paper's entries already involve extrapolation.|
|GM_GM1|0.08–1.0|Published rounded RMF couplings and source-compatible npe species.|
|RG_KDE0v1|0.08–1.8|Full asymmetric published Skyrme functional with rearrangement.|
|RG_SLY4|0.08–1.4|Full asymmetric published Skyrme functional with rearrangement.|

These executed windows are not assertions of microscopic validity over every
density, missing species or phase. No source-table pressure mismatch rejects
an independently checked implementation. Source/parameter/phase differences
remain recorded in the earlier comparison evidence.

## Density-resolved causal support

The table's adjacent-row secants are sampling diagnostics. In particular, the
first acausal secant interval need not bracket the pointwise crossing.
A separate direct derivative of the generated chemical potential locates:

| Model | c_s^2/c^2=1 density, fm^-3 | Root spread over three finite-difference steps |
|---|---:|---:|
|APR_APR|0.8689531516|1.49e-12|
|RG_KDE0v1|1.4109658933|2.67e-12|
|RG_SLY4|1.2650056715|7.96e-12|

The derivative is fourth-order, with relative steps2e-4,1e-4,5e-5. Values
immediately below/above each root and the full calculation are retained in
`data/research_completion_20261116/core/CAUSAL_BOUNDARIES.json` and its adjacent
replay script. The small step spread measures numerical root convergence, not
physical uncertainty in the EOS. BBB2 and GM1 showed no acausal sampled interval
inside the executed windows; that observation is not a universal proof.

High-density acausal extensions remain available for diagnostics. Physical use
must respect the causal interval and the declared microscopic assumptions.
Neither a high-density defect nor a CompOSE discrepancy deletes the lower-density
generation route. No stellar result or core–crust matching was calculated here.

The exact code/input identities, native checks, raw states and tables are in
`data/research_completion_20261116/core/generated/`; the build and regression
logs are in its sibling `logs/` directory.
