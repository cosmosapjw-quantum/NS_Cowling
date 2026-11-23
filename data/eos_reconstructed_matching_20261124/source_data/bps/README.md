# Reconstructed BPS outer crust for matching

`BPS_MTW1964_MS1966` solves the cold, charge-neutral ion/electron system on the 130 nuclear candidates in BPS (1971), Table 1, p. 303. The generated EOS covers total mass-energy density 10⁴ g cm⁻³ through neutron drip. The historical mass selection and the 1966 formula convention specify this reconstruction. Exact archival equivalence would require the 1965 report UCRL-11980, Appendix D, and the authors' experimental mass selection. [BPS article](https://doi.org/10.1086/151216)

## Nuclear input

`BPS_Table1_MTW1964.csv` gives A, Z and a neutral-atom mass in unified atomic mass units. A zero instructs the Fortran code to calculate the mass from Myers and Swiatecki (1966). The 23 nonzero entries come from the Mattauch–Thiele–Wapstra 1964 evaluation as reprinted on pp. 51–57 of the January 1970 revised *Radiological Health Handbook*. The intersection with the 130 candidates uses the evaluation's tabulated standard errors below 1,070 micro-u. This threshold approximates the historical choice of measured masses and does not recover a per-nuclide measured flag. The reprint PDF has SHA-256 `18d6b629775e0150f9dbe82026247dedcb76b4eb876e64b8a843ead156a32a23`. [MTW article](https://doi.org/10.1016/0029-5582(65)90114-8) · [Federal reprint](https://www.osti.gov/servlets/purl/4708654)

The remaining 107 mass excesses are calculated from the MS66 shell and deformation model, using the light magic number 14, the model constants on p. 34 and the lower of the stable spherical/prolate branches of Eq. (10). The choice of 14 was checked against the BPS Table 2 binding energies in the previous reconstruction. The nuclear masses entering the EOS are found as A u + atomic mass excess − Z mₑ, with rest energies in MeV and u c² = 931.49410242 MeV. The reconstruction retains atomic electronic binding. The source PDF used for the equations has SHA-256 `23d5031429ac8d635fca19db993972d1527b3be03a05f6647042faea01078561`. [Myers–Swiatecki 1966](https://doi.org/10.1016/S0029-5582(66)80001-9)

The input contract can be checked against the previous reconstruction bundle. Keep columns A and Z, in order, from `sources/ms66/ms66_bps_table1_light14.csv` (SHA-256 `67c62a5cdc6f66044c481fb65484b2f0b539c520afff94bfe9492f43fb4623f0`). Join `sources/historical64/mtw1964_z20to50.csv` (SHA-256 `7bba0854571d0a13a23e44ec23b9465b5600677f3301244e0d94540c19eac650`) on A and Z after selecting use_in_hybrid=1. Preserve six printed decimals of the selected atomic masses and write zero for other candidates. The resulting input SHA-256 is `452352a6f2438d94ff4f9c04562f01b505aea44cd94d294d214201b73bba1f73`.

## Generation and verification

From the repository root, run the following with a Fortran compiler, BLAS and LAPACK. The first argument after the model is a new output directory.

```sh
PYTHONPATH=python_frontend python -m compose_eos.cli generate BPS_MTW1964_MS1966 /tmp/bps-reconstructed-new
python data/eos_reconstructed_matching_20261124/check_bps_reconstruction.py --native /tmp/bps-reconstructed-new/generated --out /tmp/bps-comparison-new
```

The generator minimizes the Gibbs energy per baryon at fixed total pressure over all 130 candidates, including the zero-temperature relativistic electron gas and bcc Coulomb lattice. It obtains neutron drip from g = mₙ c². The default 601 pressure samples plus the two one-sided states at each of eleven crossings give 623 samples in 12 phases. Native composition labels are 201–212. A match from these outer-crust rows to a core may span the entire unsupported inner crust. The separately recorded matching diagnostics determine the validity of each construction.

`Table5_BPS_outer.csv` is the 43-row outer-crust subset of the printed BPS Table 5, converted to fm⁻³ and MeV fm⁻³. `compose68/eos.nb` and `compose68/eos.thermo` are the separately downloaded CompOSE #68 files with SHA-256 `35216571b367befd05388f804016b0950da9455cc83d96f4723cc3b0316cdc33` and `45011d7e94607825b6f00f9e70ba48470aa21652a138aafaaa6b31bc0f2f291a`. #68 attributes the outer crust to BPS. The comparison includes #68 rows 1–49 through the largest pressure in the printed outer-crust table, within its rounded precision. Six rows lie below the regenerated surface. The file contains no nuclear composition. Many of the tabulated reference pressures are shared by the printed table and the CompOSE distribution. [CompOSE #68](https://compose.obspm.fr/eos/68) · [Contributor sheet](https://compose.obspm.fr/download/1D/NS/Classical/APR/eos.pdf)

The check evaluates P and ε at each reference baryon density. Within a single layer it uses log-Hermite interpolation with d log P/d log n = Γ and dε/dn = g. In the coexistence gap it uses the constant transition pressure and the linear mixture energy. The checks record out-of-support rows and perform no extrapolation. For an independent comparison to the previous Python/SciPy reconstruction, add `--prior-model /path/to/BPS_reconstructed/models/historical_1964_1966`. The SHA-256 of the exact reference files is recorded with the resulting mass, state and crossing residuals.
