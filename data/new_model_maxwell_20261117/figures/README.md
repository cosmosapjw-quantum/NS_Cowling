# Nine finite-domain Maxwell candidates

Descriptive plot of the nine successful candidate reports in ../campaign/CAMPAIGN.json.
No EOS evaluation, fitting or thermodynamic reconstruction is performed by the plotting script.
The CSV preserves exact parsed native-report float values and their JSON paths.
The figure rounds density/pressure labels to five decimal places.

- new_model_maxwell_successes.png: 3432 by 1936 pixels, 220 dpi.
- new_model_maxwell_successes.svg: vector figure with text retained.
- new_model_maxwell_successes.pdf: one-page vector figure.
- new_model_maxwell_successes.csv: nine extracted native-report records.
- FIGURE_RECEIPT.json: source/script/output hashes, versions, units and scopes.
- PLOT_RUN.json and PLOT.*.log: actual process command and exit status.
- VISUAL_QA.json and PDFINFO.log: render/visual review record.

From the repository root, reproduce with:

    python scripts/plot_new_model_maxwell_20261117.py --campaign data/new_model_maxwell_20261117/campaign/CAMPAIGN.json --output-dir data/new_model_maxwell_20261117/figures

Dependencies: Python, numpy, matplotlib. PDF rendering for visual inspection used Poppler.
Successful outputs count source pairings; the campaign contains four new physical RMF models
and five explicitly named native parameter variants. The plot retains admitted finite-component,
thermodynamics-only scope and makes no global-ground-state claim.
