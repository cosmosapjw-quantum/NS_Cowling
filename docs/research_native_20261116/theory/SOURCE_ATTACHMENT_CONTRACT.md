# Original crust attachment contract

Base: 4de49430560930be84d9a4147d71c3c9841da2db. Read-only inspection of all thirteen pinned provider PDFs, source arrays and _CRUST_ONLY_SOURCE_SPECS on 2026-11-16. Exact model PDF hashes and measurements are in data/research_native_20261116/theory/SOURCE_ATTACHMENT_CONTRACT.json.

## Decision

The source-prescription construction and independently thermodynamic phase matching must be separate operations. **A common-tangent root is not a necessary admission gate for recreating a provider-prescribed density attachment.** The supplied primary PDFs specify BSk22 below the nominal nB=0.002 fm^-3 boundary and model-specific CLD/TF above. Requiring a new Gibbs root changes that construction. Conversely, the source PDFs do not establish exact pressure/chemical-potential continuity, interpolation, or branch-specific energy shifts at the join; these remain diagnostics and source uncertainties.

Thus a native fixed-density attachment may proceed under an explicit identity such as source_prescription_attachment, with the nominal boundary, left/right thermodynamic states, independent physical residuals, and the exact lower representation recorded. It must not claim Maxwell/Gibbs equilibrium merely because it follows the documented attachment.

## Per-model contract

Every listed pinned eos.pdf names BSk22 below 0.002 fm^-3. Each source table's last common outer row is nB=0.000269021 fm^-3 at row217. The next source row, row218, is:

| Model | Native inner functional | First inner source nB | Terminal source nB |
|---|---|---:|---:|
| SPG_M1 | GDFM + spherical CLD | .002511886432 | .07413102413 |
| SPG_M2 | GDFM + spherical CLD | .002290867653 | .08912509381 |
| SPG_M3 | GDFM + spherical CLD | .002630267992 | .07943282347 |
| SPG_M4 | GDFM + spherical CLD | .002630267992 | .09772372210 |
| SPG_M5 | GDFM + spherical CLD | .002511886432 | .07079457844 |
| GPPVA_DD2 | DD2 TF | .0020000001 | .067000002 |
| GPPVA_DDME2 | DDME2 TF | .0020000001 | .071999997 |
| GPPVA_FSU2 | FSU2 TF | .003 | .054000001 |
| GPPVA_FSU2H | FSU2H TF | .0020000001 | .086999997 |
| GPPVA_FSU2R | FSU2R TF | .0020000001 | .082999997 |
| GPPVA_NL3wrL55 | NL3wrL55 TF | .0020000001 | .082000002 |
| GPPVA_TM1e | TM1e TF | .0020000001 | .089000002 |
| GPPVA_TW | TW TF | .003 | .074000001 |

Densities are in fm^-3. Terminal source density is stored table support, not a newly derived physical transition.

For TW, the pinned PDF explicitly names BSk22 but the previously inspected official live page212 names BPS. The paper/PDF-selected BSk22 path can be implemented with that declared choice while preserving the conflict. It is not necessary to pretend the disagreement has disappeared, and a source-text conflict alone does not prove no physically valid TW EOS can be generated.

Provider links:
- SPG M1–M5: https://compose.obspm.fr/download/1D/Crust/M1/eos.pdf (replace M1 with M2–M5).
- GPPVA: https://compose.obspm.fr/download/1D/Crust/DD2/eos.pdf and corresponding model path.
- TW page: https://compose.obspm.fr/eos/212 ; FSU2R page: https://compose.obspm.fr/eos/203 .
- Method: Grill et al., PRC90,045803(2014), https://arxiv.org/abs/1404.2753 ; Scurto et al., PRD109,103015(2024), https://arxiv.org/html/2402.15548v1 .
The stored provider PDFs, not a guessed modern URL's availability, are the exact evidence objects for this inspection.

## What the source actually samples

The source table contains **no stored rows** in the open interval between .000269021 and the first model-specific density above. The word “below .002” in the PDF therefore cannot be read as proof that a resolved BSk22 lower-inner sequence is included in the archive.

This source-grid observation does not, by itself, establish an EOS defect: tabulated EOSs may delegate interpolation to users. However, no original interpolation prescription was recovered here. It must not be silently promoted to independently computed support.

The two endpoint pressures are not equal. For DD2 they decode from P=nB*Q1 to .0004981910045986 and .001960000017999996 MeV fm^-3. Across the thirteen models the second pressure is .0018700000935–.0033399999 MeV fm^-3, with the same lower endpoint pressure. Hence these particular adjacent rows are not the two endpoints of a zero-pressure-jump Maxwell coexistence segment. A density jump alone must not be labeled as a physical transition.

## Minimal truthful native construction

1. **Regenerate source sample locations:** compute native outer states at the source outer densities and native model-specific TF/CLD at the source inner densities; export all successful samples and preserve the interval-support metadata. This can finish full source-grid coverage without satisfying a newly imposed common-tangent condition. Original-author equivalence still depends on calibration and approximation provenance.
2. **Generate continuous density coverage with declared approximation:** extend the published BSk22 analytic lower representation to .002, then switch to the native model-specific inner functional at the nominal provider boundary. Record which BSk22 thermodynamic/chemical representation is used, left/right P, epsilon, mu, charge and baryon constraints, and any jump. If no source stitching details are known, label it a documented-boundary implementation choice; do not claim the exact author join was recovered.
3. **Generate an independently thermodynamic alternative:** solve a same-functional Gibbs/Maxwell construction where applicable and name it separately. This is optional research, not a replacement gate for steps1–2.

For step2, a piecewise table with a pressure or chemical jump can be exported as source-prescription data with an explicitly unsupported cross-boundary interpolation interval. Treating it as a smooth single barotropic physical EOS still needs a chosen physical or interpolation model; raw exportability does not establish that stronger claim.

SPG sigma0,b remain missing original inputs. A candidate refit gives a reproducible native CLD branch but must preserve that identity even if it covers every density. The lower attachment issue and the SPG coefficient issue are independent.

## Consequence for project policy

Remove only the unintended requirement that all original provider attachments possess a newly found Gibbs root. Preserve physical residual reporting, source identity, candidate labels, unresolved interpolation support and the exact TW conflict. Do not add branch-dependent offsets or species mass changes merely to pass a physical check on a different source construction.
