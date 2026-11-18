"""Produce a complete, origin-separated and physically scoped candidate inventory."""
from __future__ import annotations
from pathlib import Path
from collections import Counter,defaultdict
import json,csv,hashlib
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/all_origin_maxwell_css_20261118"
DOC=ROOT/"docs/all_origin_maxwell_css_20261118"
DOC.mkdir(parents=True,exist_ok=True)
read=lambda p:json.loads(Path(p).read_text())
data=read(D/"CAMPAIGN.json")
sources=read(D/"SOURCE_CATALOG.json")
whole=read(D/"whole_stars/CAMPAIGN.json")
source_by_id={s["source_id"]:s for s in sources["sources"]}
whole_by_name={s:read(D/"whole_stars"/s/"RESULT.json") for s in whole["source_model_ids"]}
variant_campaign=read(D/"HERMITED_INPUT_2X2/CAMPAIGN.json")
assert variant_campaign["status"]=="FINISHED_SCOPED"
assert variant_campaign["matchings"]==["maxwell","css"]
assert variant_campaign["interpolation_orders"]==["cubic","quintic"]
variant_by_id={}
for z in variant_campaign["outcomes"]:
    key=(z["tier"],z["crust"],z["core"])
    variant_by_id[key]=read(D/"HERMITED_INPUT_2X2"/z["tier"]/(z["crust"]+"__"+z["core"])/"RESULT.json")
assert len(variant_by_id)==7
native_campaign=read(D/"HERMITED_INPUT_2X2/NATIVE_MAXWELL/CAMPAIGN.json")
assert native_campaign["matching"]=="maxwell"
assert native_campaign["native_maxwell_accepted"]==14
native_by_id={}
for z in native_campaign["results"]:
    key=(z["tier"],z["crust"],z["core"],z["order"])
    native_by_id[key]=read(D/"HERMITED_INPUT_2X2/NATIVE_MAXWELL"/
        z["tier"]/(z["crust"]+"__"+z["core"])/z["order"]/"RESULT.json")
assert len(native_by_id)==14
assert data["status"]=="COMPLETED_WITH_EXPLICIT_ORIGIN_AND_SUPPORT_GATES"
assert len(data["cases"])==635
def write_csv(name,records,keys):
    target=D/name
    with target.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys,extrasaction="ignore")
        w.writeheader();w.writerows(records)
    return len(records)
def short_msg(z):
    for part in ("error","reason","type"):
        if z.get(part):return str(z[part])[:110]
    return ""
def reason(z):
    if z["both_status"]=="FINITE_ONLY_NOT_WHOLE":return z["whole"]["status"]
    if z["both_status"]=="BLOCKED_SOURCE":return z.get("reason","SOURCE_MISSING")
    return (z["css"].get("status") or "")+"/"+(z["maxwell"].get("status") or "")
rows=[]
dual=[]
full=[]
for c in data["cases"]:
    ca=source_by_id[c["crust_source_id"]];co=source_by_id[c["core_source_id"]]
    if ca["tier"]!=co["tier"] or ca["tier"]!=c["tier"]:
        raise ValueError("tier leakage in result")
    aux=whole_by_name.get(c["crust"]) if c["tier"]=="reconstructed" else None
    condition=bool(aux and aux["core"]==c["core"] and aux["both_conditional_whole"])
    stellar=aux if condition else None
    variants=variant_by_id.get((c["tier"],c["crust"],c["core"]),{}).get("interpolation_pair_results",{})
    mx_cubic=native_by_id.get((c["tier"],c["crust"],c["core"],"cubic"),{})
    mx_quintic=native_by_id.get((c["tier"],c["crust"],c["core"],"quintic"),{})
    item=dict(priority=c["priority"],tier=c["tier"],crust=c["crust"],core=c["core"],
        crust_source_id=c["crust_source_id"],core_source_id=c["core_source_id"],
        nuclear_functional_independent=not c["same_functional"],
        same_nuclear_functional_pair=c["same_functional"],
        same_functional_core_exists=ca["same_functional_core_exists"],
        dedicated_crust_only_input=ca["role"]=="crust",
        original_outer_present=ca.get("outer_present",False),
        physical_inner_and_outer_verified=bool(c["both_status"]=="BOTH_WHOLE_PASS"),
        conditional_reconstructed_BSk22_outer_TOV_both=condition,
        css_status=c["css"].get("status"),css_transition_pressure=c["css"].get("P_t"),
        css_density_jump=(c["css"].get("n_high") or 0)-(c["css"].get("n_low") or 0)
            if c["css"].get("status")=="PASS" else None,
        css_cs2=c["css"].get("cs2"),
        maxwell_status=c["maxwell"].get("status"),
        maxwell_transition_pressure=(
            c["maxwell"]["admitted"][0]["P_t"] if c["maxwell"].get("admitted")
            else None),
        maxwell_count=c["maxwell"].get("admitted_component_output_count",0),
        interpolation_cubic_maxwell_status=mx_cubic.get("status","NOT_EVALUATED"),
        interpolation_quintic_maxwell_status=mx_quintic.get("status","NOT_EVALUATED"),
        interpolation_cubic_css_status=variants.get("cubic+css",{}).get("status","NOT_EVALUATED"),
        interpolation_quintic_css_status=variants.get("quintic+css",{}).get("status","NOT_EVALUATED"),
        interpolation_cubic_maxwell_P_t=mx_cubic.get("transition_P"),
        interpolation_quintic_maxwell_P_t=mx_quintic.get("transition_P"),
        interpolation_cubic_maxwell_delta_n=mx_cubic.get("density_jump"),
        interpolation_quintic_maxwell_delta_n=mx_quintic.get("density_jump"),
        interpolation_cubic_css_P_t=variants.get("cubic+css",{}).get("P_t"),
        interpolation_quintic_css_P_t=variants.get("quintic+css",{}).get("P_t"),
        finite_pair_status=c["both_status"],unconditional_whole_status=c["whole"].get("status"),
        failure_category=reason(c))
    rows.append(item)
    if c["both_status"]=="FINITE_ONLY_NOT_WHOLE":
        item2=dict(item)
        item2.update(crust_nuclear_functional=ca["nuclear_functional"],
                     core_nuclear_functional=co["nuclear_functional"],
                     native_maxwell_evidence=c["maxwell"].get("source_ref","new-native-fortran"),
                     css_native_verified=c["css"]["status"]=="PASS",
                     independent_source_class="cross-functional" if not c["same_functional"] else "same-functional-control",
                     whole_status_note=c["whole"].get("note",c["whole"].get("status")),
                     outer_completion_contract="E reconstructed BSk22 + anchored surface + causal lower bridge" if condition else "NOT_PASSED")
        dual.append(item2)
        if condition:
            record=dict(tier=c["tier"],crust=c["crust"],core=c["core"],
                 outer="E reconstructed BSk22, conditional surface and causal exponential lower bridge",
                 source_native_crust_has_outer=False,same_functional_core_exists=ca["same_functional_core_exists"])
            for method in ("maxwell","css"):
                z=stellar["methods"][method]
                record[method+"_TOV_status"]=z["status"]
                record[method+"_R_km"]=z["admitted"]["stellar"]["R_km"]
                record[method+"_mass_msun"]=z["admitted"]["stellar"]["mass_msun"]
                record[method+"_n_join_fm3"]=z["admitted"]["join_n_fm3"]
                record[method+"_transition_p"]=z["admitted"]["transition_P"]
                record[method+"_density_jump"]=z["admitted"]["density_jump"]
            full.append(record)
rows.sort(key=lambda x:(x["priority"],x["tier"],x["crust"],x["core"]))
dual.sort(key=lambda x:(x["priority"],x["tier"],x["crust"],x["core"]))
full.sort(key=lambda x:(x["crust"],x["core"]))
column_keys=list(rows[0])
write_csv("ALL_CANDIDATE_PAIRS.csv",rows,column_keys)
write_csv("DUAL_FINITE_SUCCESS.csv",dual,list(dual[0]))
matrix=[dict(tier=x["tier"],crust=x["crust"],core=x["core"],
    same_functional_core_exists=x["same_functional_core_exists"],
    cubic_maxwell_P=x["interpolation_cubic_maxwell_P_t"],
    quintic_maxwell_P=x["interpolation_quintic_maxwell_P_t"],
    cubic_maxwell_delta_n=x["interpolation_cubic_maxwell_delta_n"],
    quintic_maxwell_delta_n=x["interpolation_quintic_maxwell_delta_n"],
    cubic_css_P=x["interpolation_cubic_css_P_t"],
    quintic_css_P=x["interpolation_quintic_css_P_t"],
    cubic_maxwell_status=x["interpolation_cubic_maxwell_status"],
    quintic_maxwell_status=x["interpolation_quintic_maxwell_status"],
    cubic_css_status=x["interpolation_cubic_css_status"],
    quintic_css_status=x["interpolation_quintic_css_status"],
    strict_whole_NS=x["physical_inner_and_outer_verified"],
    conditional_BSk22_outer_TOV_both=x["conditional_reconstructed_BSk22_outer_TOV_both"])
  for x in dual]
write_csv("MATCHING_2X2_VERIFIED.csv",matrix,list(matrix[0]))
write_csv("CONDITIONAL_WHOLE_TOV.csv",full,list(full[0]))
inventory=[]
paired_results={(x["tier"],x["crust"]):0 for x in rows}
for x in dual:
    paired_results[(x["tier"],x["crust"])]=paired_results.get((x["tier"],x["crust"]),0)+1
for x in sources["sources"]:
    if x["role"]!="crust":continue
    inventory.append(dict(priority=x["priority"],tier=x["tier"],source_id=x["source_id"],
        model=x["model"],source_class=x["provenance_class"],
        nuclear_functional=x.get("nuclear_functional"),
        same_functional_core_exists=x["same_functional_core_exists"],
        core_independent_crust_only=x.get("core_independent_model",False),
        original_outer_present=x.get("outer_present",False),
        declared_inner_and_outer=x.get("full_outer_inner",False),
        raw_status=x["status"],
        ready_for_both_matching=bool(x["status"]=="SOURCE_SNAPSHOT_AVAILABLE"),
        verified_dual_finite_partner_count=paired_results.get((x["tier"],x["model"]),0),
        artifact_path=x.get("source_path") or x.get("nuda_path"),
        external_source_raw_sha=x.get("raw_sha256"),
        important_limit="Sealed native source has no outer crust" if x["tier"]=="reconstructed" else
                        "No full projected thermodynamic EOS with a supported outer crust" if x["source_id"].startswith("external_nuda:") else "Selection may leave phase holes or finite support"))
inventory.sort(key=lambda x:(x["priority"],x["tier"],x["model"]))
write_csv("INDEPENDENT_CRUST_INVENTORY.csv",inventory,list(inventory[0]))
summary=dict(schema="ns_cowling.dual_matching_same_origin_exhaustive_scoped.v1",
    source_head="bc09c87010fccb280c455d765e6239d3910fb0db",
    scope="all materialized sealed EOS candidates from earlier cohorts, new XMLSLZ provider and 23 nucleardatapy raw-only crust sources, core-independent priority",
    number_sources=len(sources["sources"]),
    number_crust_sources=len(inventory),
    number_pairs_total=len(rows),
    matching_methods=["maxwell","css"],
    input_interpolation_methods=["cubic","quintic"],
    interpolation_matrix_2x2_total=len(variant_by_id)*4,
    actual_native_maxwell_interpolated_source_accepted=native_campaign["native_maxwell_accepted"],
    actual_native_css_interpolated_source_accepted=sum(
      variants[name]["status"]=="PASS_DERIVED_INTERPOLATED_SOURCE"
      for record in variant_by_id.values()
      for variants in (record["interpolation_pair_results"],)
      for name in ("cubic+css","quintic+css")),
    performed_native_css=sum(x["css_status"] in ("PASS","REJECTED","FAILED_CONSISTENCY") for x in rows),
    finite_dual_success=len(dual),
    local_hermite_cubic_quintic_crosscheck='7 source pairs each with one local admissible root from each Hermite order, see HERMITE_ROOT_CROSSCHECK.json',
    true_original_full_whole_success=sum(x["physical_inner_and_outer_verified"] for x in rows),
    conditional_reconstructed_whole_TOV_success=len(full),
    source_isolation=bool(all(x["tier"] in ("original","corrected","reconstructed") for x in rows)),
    status_counts=dict(Counter(x["finite_pair_status"] for x in rows)),
    by_tier=dict(Counter(x["tier"]+"|"+x["finite_pair_status"] for x in rows)),
    root_search_limit="reused prior production native Maxwell for 195 sealed legacy pairings; new pairs use vectorized cubic P(mu) preliminary roots and native Fortran, negative screen not a mathematical no-root certificate",
    strict_whole_definition="No unverified phase gaps and actual outer crust from same source; only source-checked full surface",
    conditional_whole_definition="Source inner + E pre-existing code-reconstructed BSk22 outer + anchored surface with explicit causal exponential connector, TOV 1.4M verified; not original fully unified micro EOS",
    external_prioritized_but_blocked="Negele–Vautherin has no P; Gogny and GMRS raw files do not include validated unified outer, equilibrium µ normalization/CompOSE format",
    nuda_wheel_sha256=sources["nuda_evidence"]["sha256"],
    nuda_package=sources["nuda_evidence"]["package"],
    catalogue_files=["ALL_CANDIDATE_PAIRS.csv","DUAL_FINITE_SUCCESS.csv","CONDITIONAL_WHOLE_TOV.csv","INDEPENDENT_CRUST_INVENTORY.csv","MATCHING_2X2_VERIFIED.csv"],
    dual_finite=[{k:x[k] for k in ("tier","crust","core","independent_source_class","same_functional_core_exists",
                                 "css_transition_pressure","maxwell_transition_pressure",
                                 "finite_pair_status","unconditional_whole_status","conditional_reconstructed_BSk22_outer_TOV_both")} for x in dual])
(D/"RESULTS_SUMMARY.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
lines=["# Maxwell ∩ CSS with cubic/quintic Hermite: source-origin independent census","",
 "## Source-protected result (2026-11-18)","",
 f'**{len(sources["sources"])} distinct source records, {len(rows)} same-origin pairing candidates**.',
 "Tier mapping: **original** = original provider CompOSE plus original nuda raw crust; **corrected** = energy-primary public thermodynamic reconstructions; **reconstructed** = independently generated local code outputs.",
 "A unique role-qualified source ID is retained for each; there is NO original × corrected or provider × reconstructed pairing.",
 f'All {len(rows)} candidate rows are available in ALL_CANDIDATE_PAIRS.csv, prioritized by independent crust.',
 f'Native finite Maxwell and native CSS both succeeded on **{len(dual)}** pairings. Matched output physical claims are scoped to their actually valid finite source ranges.',
 f'**Strict unconditional full inner+outer crust EOS passing both: {summary["true_original_full_whole_success"]}.**',
 f'**Conditional code-completed inner+outer stellar EOS passing both and M=1.4M⊙ TOV: {len(full)}.**',
 "",
 "**Matching 방법은 Maxwell와 CSS만 비교한다.** Cubic/quintic은 EOS 입력자료를 P(mu)로 표현하는 보간 차수로만 사용한다. 보간법 자체가 세 번째 matching construction이 아니다. 검증된 7조합의 2×2 자료는 HERMITED_INPUT_2X2/CAMPAIGN.json과 DUAL_FINITE_SUCCESS.csv에서 확인할 수 있다.", "",
 "## All finite dual-matching successes (7 pairs)",
 "",
 "| Source tier | Crust | Core | Same-functional core exists? | Maxwell P_t [MeV fm⁻³] | CSS P_t [MeV fm⁻³] | Cubic+Maxwell / Quintic+Maxwell | Cubic+CSS / Quintic+CSS | Whole status |",
 "|---|---|---|---|---:|---:|---|---|---|"]
for x in dual:
    line=(" | ".join([x["tier"],x["crust"],x["core"],x["same_functional_core_exists"],
                      f'{x["maxwell_transition_pressure"]:.6g}' if x["maxwell_transition_pressure"] else "historic validated",
                      f'{x["css_transition_pressure"]:.6g}',
                      " / ".join(x[k] for k in ("interpolation_cubic_maxwell_status","interpolation_quintic_maxwell_status")),
                      " / ".join(x[k] for k in ("interpolation_cubic_css_status","interpolation_quintic_css_status")),
                      "conditional BSk22+TOV" if x["conditional_reconstructed_BSk22_outer_TOV_both"] else x["unconditional_whole_status"]]))
    lines.append("| "+line+" |")
lines += ["","**The dual-success set has no pair sharing the same nuclear functional**. The 'same-functional core exists' column identifies the external/unified functional controls, not identity with the actual matched core.","",
  "**2×2 source-grounded comparisons:** 7 native admissible finite input pairs × (cubic/quintic P(mu) source representation) × (Maxwell/CSS matching) = 28 native matching evaluations. Maxwell was independently re-run by the Fortran transition generator on phase-restricted, Hermite-interpolated inputs (14/14 passed). CSS was likewise re-run with the native CSS generator (14/14 passed). The earlier root-only local Hermite scan remains an independent diagnostic. Original source knots were retained, readback was checked. These results do not prove complete unified microscopic EOSs.","",
  "### Exactly two matching methods × two EOS interpolation orders","",
 "| Crust/core | Cubic+Maxwell P_t | Quintic+Maxwell P_t | Cubic+CSS P_t | Quintic+CSS P_t |",
 "|---|---:|---:|---:|---:|"]
for x in matrix:
    lines.append(f'| {x["tier"]} / {x["crust"]}+{x["core"]} | {x["cubic_maxwell_P"]:.9g} | {x["quintic_maxwell_P"]:.9g} | {x["cubic_css_P"]:.9g} | {x["quintic_css_P"]:.9g} |')
lines += ["","Note: the two Maxwell values in each row are from native Fortran matched EOSs on Hermite-interpolated source representations, not merely a sign-change root. No cubic_F or quintic_F smooth connector is used as a matching method. Conditional TOV below was performed on source-matched constructions with E outer closure, not on all of the newly interpolated 2×2 variants.","",
 "## Conditional whole-star verification (not a unified-microscopic claim)","",
 "| Reconstructed crust+core | Maxwell R_1.4 [km] | CSS R_1.4 [km] | Same-source E outer and causal lower join? |",
 "|---|---:|---:|---|"]
for x in full:
    lines.append(f'| {x["crust"]} + {x["core"]} | {x["maxwell_R_km"]:.7f} | {x["css_R_km"]:.7f} | Yes, explicit conditional closure |')
lines += ["",
 "The five native SPG crust sources lack an outer crust. For these conditional whole models we retained their computed inner crust and appended the *previously reconstructed* BSk22 outer/anchored surface plus an explicit causal low-density connector, without altering source inner knots. Their mass was actually integrated to M≈1.4M⊙. This is not proof that the entire star is generated under one microscopic functional and not proof of global TF minimality.",
 "GPPVA_FSU2+GM1 failed the full-star consumer since inherited source phase boundaries cannot be smoothed across; it is retained as a finite-interface success only.",
 "Original provider XMLSLZ(TM1)+GM1 has an outer source but an invalid/unadmitted hole between accepted crust source components, thus full crust is not approved. It is a **finite Maxwell+CSS success, not a whole-neutron-star success**.",
 "",
 "## Prioritized independent crust candidates and constraints","",
 f'The inventory covers **{len(inventory)} crust source entries** with same-functional reference and original low-density support flags.',
 "External Negele–Vautherin (1973) raw table contains baryon density and E/A-m_n but not pressure P or a verified full EOS; original-table Maxwell cannot be calculated from it without explicitly reconstructing pressure and fixing integration constants. It stays BLOCKED, not a negative physics result.",
 "External Gogny D1S/D1M/D1MS and GMRS 2022 contain source n,P,E/A-m_n but do not have a complete verified same-source outer crust down to the surface nor independently verified total-energy μ and exact CompOSE slice. The original raw corpus is preserved. The external paired EOS is **not fabricated or counted**.",
 "The external BPS crust historical reference is outer-only, BBP is inner-only; without a full matched physics and normalization contract neither alone satisfies the user-mandated inner+outer condition.",
 "Previously user-uploaded MYN13_RHF_TF and OMEG1_EoS are inventoried as original ZIP-only sources with exact source-file SHA. Earlier direct ZIP audit found four mathematical crossings but no physically adopted crust/core Maxwell root. They have no sealed original CompOSE source slice on the remote workstation; therefore new cross-combinations are explicitly blocked, not called no-coexistence.",
 "All 13 provider-original SPG/GPPVA crust tables contain an attached BSk22 outer segment and source inner crust, although native original-source admissible inner components can have holes. Their corrected counterparts retained unchanged outer rows but have energy-primary reconstructed inner pressure; source selection/first-law assumptions differ and are never combined with originals.",
 "The five XMLSLZ original provider crust records are kept as original tier, with their own support provenance. Matched cores are original source of the same tier only.",
 "",
 "## Scientific decision rules and limitations","",
 "1. **오직 Maxwell와 CSS만 matching construction으로 사용한다.** Cubic/quintic Hermite는 이 두 matching에 공급하는 EOS의 P(mu) 입력 보간법이다. 7개 유한 성공 조합에서 (cubic,quintic)×(Maxwell,CSS)의 28개 조건을 별도로 계산했다. CSS는 Hermite spline 자체가 아니다.",
 "2. A finite Maxwell root + CSS PASS is not a complete NS EOS. Requiring an approved outer crust, complete inner phase chain and 1.4M⊙ TOV shrinks the viable set from 7 finite to 5 *conditional* and 0 strict source-native whole-star results.",
 "3. Generic origin-tier mixing prohibited. A newly created thermodynamic interpolation from an original source is a new reconstructed object, never an original-provider value. No pressure/root location/extrapolation is fabricated to rescue a rejected combination.",
 "4. The coarse cubic P(μ) screen identifies positive roots; for new unscreened pairings a screen-negative result is **not an exhaustive nonlinear root nonexistence proof** (tangential root, multiple narrow roots, invalid source branches). An old rigorous native Maxwell 195 pair campaign was reused with exact source identity.",
 "5. Physical assumption: single conserved baryon charge and local charge neutrality, beta-equilibrated cold barotropes; additional pasta phase competition, formation history, surface tension/interface tension, superfluidity, frozen composition, radial stability and full-GR oscillations remain unresolved.",
 "",
 "## Evidence and reproducibility","",
 "- data/all_origin_maxwell_css_20261118/SOURCE_CATALOG.json (source/provenance/input support tier)",
 f"- data/all_origin_maxwell_css_20261118/CAMPAIGN.json ({len(rows)} pair results and actual CSS/new Maxwell receipts)",
 "- data/all_origin_maxwell_css_20261118/whole_stars/CAMPAIGN.json and */RESULT.json (6 bounded attempts, 5 conditional dual TOV)",
 "- data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json (195 prior independent native Maxwell pairings)",
 "- nucleardatapy 1.0.2 original external distribution wheel SHA recorded in SOURCE_CATALOG, downloaded separately under scratch, not altered or silently republished",
 "",
 "## Next scientific boundary",
 "To qualify the full original/corrected crust models for a strict whole-star Maxwell+CSS study, source gaps must be resolved by author-supplied physically consistent low-density tables or independently solved source EOS with identical energy/charge conventions. The 5 conditional SPG+PK1 stars are usable phenomenological comparison baselines, not a validation of a shared microscopic crust-core functional.",
 ]
(DOC/"RESULTS_KO.md").write_text("\n".join(lines)+"\n")
print("GENERATED",len(rows),"PAIRS",len(dual),"DUAL_FINITE",len(full),"CONDITIONAL",len(inventory),"SOURCE_CRUST_ENTRIES")
