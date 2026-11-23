# Fixed-cell orbital HF reconstruction

The eleven densities and cell baryon counts follow Table 3 of Negele and Vautherin (1973), DOI 10.1016/0375-9474(73)90349-7. The existing orbital solver uses the density-matrix-expansion functional of Negele and Vautherin (1972), DOI 10.1103/PhysRevC.5.1472. Mode 0 retains the printed polynomial coefficients. Mode 1 multiplies the printed D03 coefficient by 100, a conditional transcription hypothesis. The derivation, independent rounding counterexamples and source transcriptions are in the existing `eos_nv_dme_20261123` evidence.

The file `native_triplets_both_modes.csv` contains the actual 66 runs at the eleven reference densities with fixed A and charges Z₀−2, Z₀, Z₀+2 for both coefficient choices. The triplet table `d03_triplet_native` selects the lowest energy per baryon of the three converged mode-1 orbitals at each density on a 96-point radial grid. The archived source values include the occupied neutron and proton kinetic frontiers and the total electron chemical potential. In each Wigner–Seitz cell, the pressure comes from a fixed-particle-number dilation. The binary results include the independently evaluated central dilation error. Energies are measured relative to the neutron rest mass, including the proton–neutron mass difference and electron energy.

In the native table Q3 encodes the occupied neutron chemical potential, Q4 the occupied proton–neutron difference and Q5 the total weak-reaction affinity with the electron contribution. The NML records `beta_equilibrated=F`. At an integer shell the addition and removal energies can differ. The full interval from these orbitals is in the cell CSV. Five of the eleven selected triplets bracket zero. The 1D Euler enthalpy differs from the occupied neutron level by −0.100 through +1.708 MeV. The pressure–energy and convexity secants pass across the ten adjacent density pairs. Cell A, global charge equilibrium and continuum convergence remain separate questions.

The wider fixed-A search tests 109 distinct clustered mode-1 charge points across the eleven densities. The energies of three selected triplets decrease at other charges. On this mesh all eleven extended winners lie below both immediate integer neighbors. Six bracket zero in their frontier interval. Replaying each triplet winner from a diffuse seed at mesh 96 gives energies within 2.1×10⁻¹² MeV per baryon. Going from mesh 96 to 160 at fixed charge shifts the energy by 0.0067–0.1251 MeV per baryon and the pressure by up to 5.94%. The original historical EOS and the equilibrium composition are unresolved.

The clustered literal run fails its density-fit-domain line search in all 33 primary charge cells. Seed and mesh changes give eight converged literal runs across four distinct sites. At mesh 96, the uniform seed converges at sites 1, 2, 3 and 10. The weak-affinity intervals of the first three are positive; the two adjacent energy secants decrease and lie below both endpoint Euler enthalpies. The remaining converged site has six missing sites between it and site 3. The isolated cells are preserved in `literal_partial_cells.csv`. The failed intervals are in `literal_partial_secants.csv`.

## Reproduction

Build the repository's native programs with `fortran_eos/build.sh BUILD`. The campaign runner accepts plans in the `plans` directory. Use the original input under `data/eos_original_papers_20261123/research_nv/nv1973_original_source.dat`. For each plan in the mapping below, place its copy under the same name in a fresh campaign root and run the existing `scripts/run_nv_hf_source_campaign_20261123.py` with `--plan ROOT/PLAN_FILE --executable BUILD/nv_hf_probe --source INPUT --out ROOT/CAMPAIGN --workers 4 --timeout 360`.

| Plan file | Campaign subdirectory |
|---|---|
| PLAN.json | site96 |
| LITERAL_ALT_PLAN.json | literal_alt |
| CHARGE_EXT_PLAN.json | charge_ext |
| FOLLOWUP_PLAN.json | followup |

Run `python scripts/collect_nv_orbital_diagnostics.py --campaign-root ROOT --output NEW_OUTPUT --importer BUILD/compose_eos_import_cold_table`. The collector preserves the source grid, projects through the native CompOSE writer and reader, and compares Q5, μn, ε and P against the orbital results. Each numerical output digest is recorded in `all_native_cells.csv`. Binary and source digests are in `orbital_summary.json`. Failed orbital results have empty energy and frontier columns. The separately reported raw iteration peaks describe trial densities.
