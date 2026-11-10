# C7X Independent Source Audit

Exploratory 16-source audit, not official C7 admission.

46/46 raw files matched pinned Git blob SHA-1. Strict Q3/Euler gate: 9 pass, 7 fail. Q5 gate: RG_SLY4 and RG_KDE0v1 fail. Robust acausal core tails: APR_APR and RG_KDE0v1. All 55 pairs have intersecting tabulated mu_B ranges, but all require rest-mass rebase; no root search, TOV or Cowling calculation was run.

The full 60,377-byte bundle is in Dropbox at /NS_Cowling/artifacts/sharp_maxwell/C7X_INDEPENDENT_SOURCE_AUDIT_20261110T123004Z.zip. SHA-256: 42c3f149eca8f5bf82abba524d4ffbaf4bd4fe0ac841475ac81ad087f61b49ed.

Native project build and native parser tests are blocked because no usable local checkout could be materialized. Generic Fortran smoke test passed; synthetic tests passed 8/8. PR #4 and the canonical 65-pair manifest are unchanged.