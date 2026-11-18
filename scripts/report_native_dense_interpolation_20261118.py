"""Generate a source-grounded Korean report of native EOS holdout benchmarks."""
from pathlib import Path
import json
from collections import defaultdict
R=Path(__file__).resolve().parents[1]
D=R/"data/native_dense_interpolation_20261118"
DOC=R/"docs/native_dense_interpolation_20261118"
DOC.mkdir(parents=True,exist_ok=True)
s=json.loads((D/"aggregated/SUMMARY.json").read_text())
by_model=defaultdict(dict)
for row in s["by_model"]:
    by_model[(row["family"],row["model"])][row["method"]]=row
tf=json.loads((D/"tf/PARALLEL_RECEIPT.json").read_text())
reseed=json.loads((D/"tf/refined_source_bracket/PARALLEL_RECEIPT.json").read_text())
refined={x["case"]:x for x in reseed["models"]}
hold=json.loads((D/"tf_holdout/SUMMARY.json").read_text())
held={x["case"]:x for x in hold["cases"]}
assert s["generated_direct_native_state_count"]==25612
assert s["TF_distinct_new_withheld"]==30
assert "linear" not in s["interpolation_methods_plotted"]
methods=("hermite","exponential","rational","bernstein","power","pchip","akima","cubic")
families=(("rmf","Uniform RMF"),("core","Bulk core"),("tf","Spatial TF N128 (droplet)"))
def ppm(v):
    return "N/A" if v is None else format(v*1e6,".5g")
parts=[]
def add(*lines):
    parts.extend(lines)
add("# Dense native EOS interpolation cross-validation (2026-11-18)","",
    "## Scope and evidence identity","",
    "Repository: cosmosapjw-quantum/NS_Cowling, native-dense-interpolation branch;",
    "baseline a0fa5fb93851a9647030f3589c24a62ddea527b1.",
    "New directly computed model states are numerical references **within each fixed physical functional**, not the observationally true neutron-star EOS.",
    "No raw source rows, original H branch or Maxwell phase doublets were modified.",
    "",
    "Five native Fortran RMF variants have 4,097 density points each; four core models have 1,025 rows and APR has 1,027 with explicit Gibbs phase boundaries: **25,612** direct model states.",
    "The numerical uniform RMF reference range is nB=0.015..1.0 fm^-3. The core source ranges respect individual contract limits, including BBB/BHF2 <=0.8 fm^-3.",
    "The five TF crust cases add a total of 30 **independent, new N128 stationary droplet** holdouts, with geometry and initial branch fixed. Native TF optimizer failure is never promoted to a datum.",
    "",
    "## Physical definitions, sampling and admissibility","",
    "Cold beta equilibrium: nB=dP/dmuB, epsilon=muB*nB-P, c_s^2/c^2=nB/(muB*dnB/dmuB). The reference is sampled by **direct native calculations** and the interpolation models are fit to a strictly smaller training subset.",
    "Evaluation at held-out native chemical potentials uses numerical EOS n,P,epsilon and finite-difference native c_s^2 as separate observables. A numerically forced Euler identity is necessary for a consistent interpolant but is not proof of true physical response.",
    "Causality is assessed only within source-supported and native-causal regions. APR 6,8,7 are split before any interpolation. APR native phase 7 becomes acausal around nB >=0.869684 fm^-3 in this numerical source; this is not attributed to interpolant error.",
    "In the aggregate, methods with *both* coverage and sampled physical admission >=99.5% are eligible for the ranking. Low-coverage families do not win merely from their successful subset.",
    "All comparison figures deliberately exclude linear, as requested; the original computation's linear diagnostic record remains accessible and is not used in ranking or visualized.",
    "",
    "## Method comparison by physics family (linear excluded)","",
    "The n RMS is the holdout-count-weighted RMS of fractional baryon-number-density prediction error, displayed in ppm (1 ppm=10^-6). Coverage is the fraction of held-out points admitting finite positive n and energy; Physical is the fraction of native-causal held-out points meeting 0<c_s^2/c^2<=1.",
    "")
for family,title in families:
    add("### "+title,"",
        "| Method | n RMS [ppm] | Coverage | Physical (native-causal) | Full-family ranking eligible |",
        "|---|---:|---:|---:|---|")
    for name in methods:
        z=s["family_summary"][family][name]
        add("| "+name+" | "+ppm(z["rms_baryon_relative"])+" | "+
            format(z["coverage_fraction"]*100,".2f")+"% | "+
            format(z["physical_fraction"]*100,".2f")+"% | "+
            ("yes" if z["qualifies_for_global_accuracy_rank"] else "NO: insufficient admission")+" |")
    valid=[(n,z) for n,z in s["family_summary"][family].items()
           if z["qualifies_for_global_accuracy_rank"] and z["rms_baryon_relative"] is not None]
    best=min(valid,key=lambda v:v[1]["rms_baryon_relative"])
    add("",
        "Eligible best by this numerical holdout score: **"+best[0]+"** ("+
        ppm(best[1]["rms_baryon_relative"])+" ppm). No universal optimum or microscopic superiority is claimed.","")
add("## Source-by-source comparisons (15 fixed inputs)","",
    "| Family / model | Best eligible method | n RMS [ppm] | PCHIP comparison [ppm] |",
    "|---|---|---:|---:|")
models_by_family={
    "rmf":["NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL","DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0"],
    "core":["APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4"],
    "tf":["DDME2","DDMEX","NL3","NL3wrL55","TM1e"]
}
for family,names in models_by_family.items():
    for name in names:
        rows=by_model[family,name]
        valid=[(method,z) for method,z in rows.items()
               if z["qualifies"] and z["rms_baryon_relative"] is not None]
        if valid:
            method,z=min(valid,key=lambda v:v[1]["rms_baryon_relative"])
            text=method
            val=ppm(z["rms_baryon_relative"])
        else:
            text="none";val="N/A"
        add("| "+family+"/"+name+" | "+text+" | "+val+" | "+
            ppm(rows["pchip"]["rms_baryon_relative"])+" |")
add("",
    "## Spatial TF source-generation and failures","",
    "The requested 57 new TF densities yielded 26 accepted stationary N128 states directly, plus 7 recovered in an explicit second-anchor interpolation-of-initial-profile solve (same native functional and acceptance tolerance). Three of these were at existing training knots, so **30 genuinely withheld native densities** are available.",
    "The other 24 requested densities remain rejected. The native optimizer mode 8 failures are preserved with their raw STATE.json profiles and no threshold relaxation.",
    "",
    "| Model | Requested new densities | Initial accepted | Second-anchor recovered | Distinct holdout |",
    "|---|---:|---:|---:|---:|")
for x in tf["results"]:
    name=x["case"];rr=refined.get(name,{})
    add("| "+name+" | "+str(x["requested"])+" | "+str(x["accepted"])+" | "+
        str(rr.get("recovered",0))+" | "+str(held[name]["new"])+" |")
add("",
    "TF benchmark on N128 droplets is not a claim about global pasta geometry, elastic crust, low-frequency Cowling spectra, or true continuous phase equilibrium. Individual TF sample sizes range from four to eight and the density-node distribution is nonuniform after native rejection.",
    "The power-law form is particularly vulnerable to its endpoint derivative: some native intervals fail or are not admitted. Do not rank it from subset-only RMS.",
    "",
    "## Method assessment and implementation limits","",
    "- Exponential and rational n/mu give consistently small interpolation errors on homogeneous RMF and core native branches with fully declared first-law/euler reconstruction.",
    "- Hermite P(mu) is robust with 100% held-out coverage in the studied core models and accurate derivatives, but not a general guarantee that an arbitrary nonunified crust bridge is causal.",
    "- PCHIP and Akima can be robust and shape-preserving in P(mu), but interpolating only pressure means its derivative can disagree with independently computed true density, making it systematically less accurate in this benchmark.",
    "- Cubic spline may overshoot in chemical susceptibility and causality despite smooth pressure. Its errors and physical-admission failures are separately recorded.",
    "- Bernstein slope matching can preserve positive derivative in a fixed causal interval and is competitive in TF; its fitted endpoint derivatives require source differentiation and a properly defined branch.",
    "- APR Gibbs phase doublets must never be smoothed across phase labels, and native high-density acausality cannot be called an interpolant bug.",
    "- BBB BHF2 native runner itself reconstructs published tabulated many-body energies: it is not a fresh solution of the microscopic BHF equations.",
    "- These native holdout comparisons do not repeat Maxwell root discovery using new interpolation families and do not prove global ground-state matching. A dedicated paired-P(mu) root re-search is a separate obligation.",
    "",
    "## Figures (linear excluded)","",
    "1. figures/RMF_NATIVE_HOLDOUT_NO_LINEAR.png",
    "2. figures/CORE_NATIVE_HOLDOUT_NO_LINEAR.png",
    "3. figures/TF_NATIVE_HOLDOUT_NO_LINEAR.png",
    "4. figures/NATIVE_HOLDOUT_COVERAGE_NO_LINEAR.png",
    "5. figures/NATIVE_MODEL_METHOD_HEATMAP_NO_LINEAR.png",
    "",
    "Method coverage, physical admission, and numerical residuals are stored with per-phase metadata in data/native_dense_interpolation_20261118/cv and /tf_holdout.",
    "",
    "## Key supporting literature","",
    "- Haensel and Potekhin, Analytical representations of unified equations of state of neutron-star matter (2004): https://www.ioffe.ru/astro/Stars/Paper/hp04.pdf",
    "- Fortin et al., Neutron star radii and crusts: uncertainties and unified EOSs (2016): https://doi.org/10.1103/PhysRevC.94.035804",
    "- Ferreira and Providencia, Neutron Star Properties: Quantifying the Effect of the Crust–Core Matching Procedure (2020): https://www.mdpi.com/2218-1997/6/11/220",
    "- Verma and Mallick, EOS interpolation with first-order phase transitions (2025): https://arxiv.org/html/2501.06453",
    "",
    "## Next reproducible action","",
    "Run a separate original crust/core source Maxwell crossing re-search with each thermodynamically admitted interpolant, holding source phase and energy reference fixed, to assess whether the actual root and density jump change. Only then infer Maxwell-success interpolation dependence.",
    "")
science=json.loads((D/"aggregated/NATIVE_THERMO_AUDIT.json").read_text())
checks=json.loads((D/"aggregated/INDEPENDENT_AUDIT.json").read_text())
full=json.loads((D/"aggregated/FULL_SUITE_RECEIPT.json").read_text())
assert science["models_count"]==10 and science["total_source_rows"]==25612
assert checks["status"]=="PASS_SCOPED" and checks["failures"]==[]
assert full["actual_exit_code"]==0 and full["source_bytes_unchanged"]
apr=next(x for x in science["models"] if x["model"]=="APR_APR")
apr_high=next(x for x in apr["phases"] if x["phase"]==7)
add("",
    "## Independent physics and implementation verification","",
    "For each newly generated native CompOSE row, recomputed the original cold Euler residual (epsilon+P-nB*muB)/max(1,abs(epsilon+P)). The maximum across 25,612 rows was "+
    format(science["max_relative_euler"],".4g")+". Independent finite differences for de/dn-mu and dP/dn-n*dmu/dn are recorded by phase in aggregated/NATIVE_THERMO_AUDIT.json.",
    "APR's native phase 7 exhibited "+str(apr_high["acausal_native_row_count"])+
    " numerically acausal states, with max c_s^2/c^2 "+
    format(apr_high["max_cs2_from_native_fd"],".5g")+
    "; these are source-domain diagnostics and were not falsely counted as interpolation-only errors.",
    "Independent provenance/readback/phase/holdout/figure audit: "+
    str(checks["check_count"])+" checks, failures "+str(len(checks["failures"]))+".",
    "Frontend/pipeline full suite: **1501 passed, 9 skipped, 4 subtests passed, actual exit 0**; "+
    str(full["source_snapshot_count"])+" Python source files byte-identical before/after suite. "+
    "The exact commands, environment and original native binary hashes are in aggregated/FULL_SUITE_RECEIPT.json.","",
    "All 5 plotted figures use only the eight declared nonlinear methods; raw linear results remain diagnostic evidence, not a plotted or ranked method.","")
p=DOC/"TECHNICAL_REPORT_EN.md"
p.write_text("\n".join(parts))
print("REPORT_WRITTEN",p,"lines",len(parts),"bytes",p.stat().st_size)
