# Final release validation

Status: PASS for the declared numerical and file-integrity checks.

- 234 unique EOS rows; 78 per high-density core family.
- All f-mode eigenfunctions have zero internal nodes in both W and V.
- All p1-mode eigenfunctions have one internal node in both W and V.
- All 234 thermodynamically projected tables pass the native EOS cell audit.
- 234 original normalized EOS tables and 234 projected EOS tables are included separately.
- The data-only ZIP passes unzip -t.
- REVTeX 4.2 paper compiles to a 5-page PDF; PDF text extraction succeeds; rendered pages were visually inspected.
- Maximum relative target-mass residual: 9.914e-13.
- Maximum absolute normalized matching residual: f 4.560e-08, p1 3.089e-10.

Nonfatal typesetting warnings remain for two displayed equations (approximately 1.2 pt and 7.2 pt overfull hboxes). The rendered pages show no clipped or overlapping content. The author line is intentionally a placeholder and must be replaced before submission.

Physical scope: ell=2 polar-fluid modes in the relativistic Cowling approximation; cold beta-equilibrated barotropic perturbation closure; finite positive pressure cutoff. These are not full-GR mode frequencies, do not provide gravitational-wave damping times, and omit frozen-composition, rotation, magnetic, elastic-crust, superfluid and sharp-interface mode physics.
