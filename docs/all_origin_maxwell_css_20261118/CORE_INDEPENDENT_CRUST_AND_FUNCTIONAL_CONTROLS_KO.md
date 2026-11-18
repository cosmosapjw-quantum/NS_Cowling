# Crust independence and same-functional controls

## Distinct concepts

**Crust-only dataset**: the published EOS table is truncated at the crust–core transition or otherwise distributed as a crust model, not necessarily developed without core functional input.

**Independent crust/core functional**: the chosen crust and selected matched core use different nuclear energy-density functionals, parameter sets and high-density calculations. This does **not** establish that the crust by itself is independent of all core functionals, because same-functional unified versions can be published separately.

**Strict whole neutron-star EOS**: the crust source has a verified outer-to-inner EOS through the actual new transition, without interpolation across rejected phase intervals, and the resulting cold barotrope is TOV admissible. A finite doublet alone does not meet this definition.

**Conditional reconstructed whole star**: a locally reconstructed native inner crust is joined to a separately reconstructed BSk22 outer crust and an explicitly causal phenomenological low-density connector, with a bounded TOV calculation; this is not a single microscopically unified EOS.

## Prioritized crust candidates actually inventoried

| Corpus | Crust models | Original thermodynamic data | Source outer-to-inner status | Matching status |
|---|---:|---|---|---|
| Negele–Vautherin 1973 | 1 | number density, energies, composition; no published P column in downloaded nuda raw | inner only | BLOCKED missing pressure, not NO_CROSSING |
| BSk 2018 | 4 | raw numeric multi-column tables (BSk22/24/25/26), metadata and normalization not converted to verified cold barotrope | supplied inner crust, not source-validated full surface | BLOCKED not ingested as original CompOSE |
| Gogny MVCD 2020 | 3 | raw density, P, E/A−mn for D1S/D1M/D1MS | inner only, no matched outer | BLOCKED source-complete whole EOS |
| GMRS 2022 | 15 | raw density, P, E/A−mn for varied Skyrme/χEFT-calibrated meta-models | low-density sampling exists, but no verified surface thermodynamic closure | BLOCKED source/normalization/surface contract |
| Baym–Pethick–Sutherland / BBP | 2 | named historical candidates, no sealed standalone table integrated | BPS outer-only; BBP inner-only | BLOCKED incomplete candidate |
| Earlier user ZIP MYN13 / OMEG1 | 1 crust + 1 core | original four-file ZIP with known SHA-256, not a sealed remote CompOSE input | no approved Maxwell transition in prior ZIP audit | BLOCKED remote source identity |
| Original CompOSE SPG(M1–M5) | 5 | native original provider tables including original BSk22 outer | whole low-density outer plus inner exists in source, but screened inner phase components can be incomplete | zero dual-match strict whole |
| Original CompOSE GPPVA(DD2, DDME2, TW, FSU2, FSU2R, FSU2H, NL3wrL55, TM1e) | 8 | native provider cold tables with BSk22 outer; pasta geometry retained | originally has inner+outer, but phase-by-phase physical admission not global | zero dual-match strict whole |
| CompOSE energy-primary corrected | 13 | original source energy/density retained, selected P and chemical columns reconstructed from same version | original outer retained; adjusted inner thermodynamic derivative | zero dual finite matches |
| Code-generated native crust (SPG + GPPVA) | 13 | generated finite branch, explicit source-limit metadata | original native output **lacks** an outer crust | 6 dual finite; 5 conditionally completed whole |
| New original XMLSLZ(NL3, PK1, TM1, DDME-X, PKDD) | 5 | provider immutable original cold rows + geometries | source includes low-density phases, gaps may remain | 1 dual finite TM1+GM1; no full admission |

nucleardatapy raw evidence: <https://pypi.org/project/nucleardatapy/1.0.2/>. The downloaded exact wheel is not silently republished as a source file; its SHA-256 and per-model source-row hashes are captured in `SOURCE_CATALOG.json`. The source library identifies 23 crust data inputs. Representative source paper: *The nucleardatapy toolkit for simple access to experimental nuclear data, astrophysical observations, and theoretical predictions* (EPJA, 2026), DOI 10.1140/epja/s10050-025-01760-w.

## The seven verified finite pairings and the SAME-FUNCTIONAL controls

| Matched crust + core (not the same functional) | Same-functional core reference exists? | Verified control source | Interpretation |
|---|---|---|---|
| native GPPVA(FSU2) + native GM1 | **Yes**: GPPVA(FSU2) unified | <https://compose.obspm.fr/eos/215> | Different RMF parameterizations; separate corresponding unified FSU2 control |
| native SPG(M1) + native PK1 | **Yes**: SPG(M1) unified | <https://compose.obspm.fr/table_selection> | Different GDFM/RMF couplings; SPG native surface coefficients inferred from public composition |
| native SPG(M2) + native PK1 | **Yes**: SPG(M2) unified | <https://compose.obspm.fr/eos/312> | Same caveat; external unified control exists |
| native SPG(M3) + native PK1 | **Yes**: SPG(M3) unified | <https://compose.obspm.fr/table_selection> | Same caveat |
| native SPG(M4) + native PK1 | **Yes**: SPG(M4) unified | <https://compose.obspm.fr/table_selection> | Same caveat |
| native SPG(M5) + native PK1 | **Yes**: SPG(M5) unified | <https://compose.obspm.fr/table_selection> | Same caveat |
| original XMLSLZ(TM1) crust + original GM1 core | **Yes**: same-provider XMLSLZ(TM1) uniform phase 0 | original `XMLSLZ_TM1` full provider table | Different RMF parameterizations; raw source crust contains an unadmitted component gap |

Note: having a same-functional *reference* does not mean that the actually selected core was same-functional. All seven matched pairs use different nuclear functionals. The named SPG(Mi) unified model exists on the CompOSE table listing; the locally reconstructed SPG(Mi) surface parameter variant still has unresolved identity with the original authors' calibration. The latter cannot be inferred merely from a shared model name.

## Scientific boundary of the screening

For a single cold locally neutral baryonic barotrope, a Maxwell candidate must have equal `P_-` and `P_+`, equal `mu_B`, a positive density jump, latent energy `mu_t*(n_+-n_-)`, and correct pressure preference. A CSS transition instead imposes a different phase prescription. Passing both constructions as *separate alternatives* does not establish they are both microscopic phase equilibria.

**Only Maxwell and CSS are matching methods in the present experiment.** Cubic/quintic Hermite are two alternative EOS input interpolants, tested against both matching methods (2×2). Earlier cubic_F/P and quintic_F/P native smooth-matching outputs are historical context only and are **not** included as an additional matching method or counted toward any of the dual-success outcomes. Positive local Maxwell crossing scans are supplementary evidence, not a proof that all phases have been examined or that no other root exists. The source-unmodified TOV checks for 5 conditional full stars apply to those original-source Maxwell/CSS EOSs; they must not be attributed to newly interpolated 2×2 variants without separate stellar runs.

The strict user requirement remains unresolved for original-source same-functional-independent pairs: **zero unconditional original inner+outer complete NS EOSs have passed BOTH CSS and Maxwell on the sealed source support**. Conditional code-completed SPG(M1–M5)+PK1 TOV curves do pass, but their outer crust is an explicitly introduced reconstructed component.
