# Explicit thermo-only Maxwell input transport

The post-normalization actual-source campaign reached a second interface defect: the pipeline said it explicitly dropped composition for both methods, but its Maxwell branch forwarded composition-bearing inputs directly into the strict native kernel. This produced COMPOSITION_UNSUPPORTED before any coexistence search. It was an orchestration defect, not a failure of the thermodynamic phase construction.

The public construct_maxwell_transition API now accepts the explicit boolean drop_composition=False. The default retains the native rejection of composition-bearing input. With True, a private frontend transport helper makes independent copies and removes only abundance arrays, presence/applicability masks, row status and composition-availability/file flags. Zero placeholders carry false presence masks; they are not asserted zero fractions.

The helper preserves nb, all Q1–Q7, Iphase, extra quantities, axes, rest-energy references, temperature, beta-equilibrium metadata and source_dir. The existing support check runs on the original inputs before any view is created. Nadd remains subject to the native unsupported-input gate. Thus the option neither erases source phase boundaries nor silently broadens thermodynamic support.

The already-declared thermo-only matching pipeline uses the same helper explicitly for Maxwell. It preserves its original work/crust_internal and work/core_internal files and writes separate crust_thermo_internal and core_thermo_internal views. composition_transport.json records source composition presence and view hashes; the pair result records whether composition was dropped, including on later native failures. The output remains a native phase-only barotrope with two phase IDs and no fabricated abundances. Local charge neutrality is still an explicit caller/pipeline assumption; the transport does not establish it from nuclear composition.

No production physical formula or Fortran acceptance criterion changed in this repair. Python performs copying, serialization and provenance only.

## Evidence

The new public/API/pipeline test file is python_frontend/tests/test_maxwell_explicit_drop_20261117.py. Before the change it produced 7 failed, 1 passed: the opt-in/helper did not exist and the supposedly thermo-only pipeline reached the native composition rejection. Afterward all 8 pass. Controls cover:

- unchanged default rejection and boolean-only opt-in;
- deep-copy preservation of every thermodynamic column, phase labels and source identity;
- input abundance/mask immutability;
- native Maxwell agreement with an equivalent composition-free analytic input;
- phase-only CompOSE write/readback without abundance presence claims;
- continued Nadd rejection;
- actual pipeline execution preserving original files and admitting the analytic pair.

One existing mocked pipeline-failure test now emits the internal CSV/metadata files that the real native read operation creates, allowing the new explicit transport stage to run before its injected native failure. Its category/stage assertion is unchanged; no physical validation is mocked by the new end-to-end test.

The native executable is unchanged from the reference-normalization repair. Red/green and affected-suite logs plus identities are under data/research_matching_20261117/native_matching/maxwell_composition/.
