#!/usr/bin/env python3
"""Build a same-paper MIC/hemolysis benchmark for B2RP variants."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MIC_PATH = ROOT / "data/processed/dbaasp_a_baumannii_mic_parsed.csv"
H50_PATH = ROOT / "data/processed/dbaasp_human_rbc_h50_candidate_records.csv"
OUTPUT = ROOT / "results/candidate_prioritization/b2rp_published_lead_benchmark.csv"
REPORT = ROOT / "results/candidate_prioritization/b2rp_published_lead_benchmark.json"
PMID = "19793185"
STRAIN = re.compile(r"\bNM\d+\b", re.IGNORECASE)


def main() -> None:
    with MIC_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        mic = list(csv.DictReader(handle))
    with H50_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        h50 = list(csv.DictReader(handle))

    mic_by_id: dict[str, list[dict[str, str]]] = {}
    h50_by_id: dict[str, list[dict[str, str]]] = {}
    for row in mic:
        if PMID in (row.get("candidate_pubmed_ids_json") or "") and STRAIN.search(row.get("target_species_raw", "")):
            mic_by_id.setdefault(row["dbaasp_id"], []).append(row)
    for row in h50:
        if row.get("candidate_pubmed_id") == PMID:
            h50_by_id.setdefault(row["dbaasp_id"], []).append(row)

    fields = [
        "dbaasp_id", "peptide_name", "sequence", "number_of_named_isolates",
        "named_isolate_mic_values_raw", "mic_unit_raw", "highest_mic_value",
        "human_rbc_50_percent_hemolysis_raw", "h50_unit_raw",
        "apparent_h50_over_highest_mic", "priority_rank_within_published_panel",
        "ratio_interpretation", "evidence_pmid",
        "interpretation_limit",
    ]
    rows_out: list[dict[str, str]] = []
    for peptide_id in sorted(set(mic_by_id) & set(h50_by_id), key=int):
        mics = mic_by_id[peptide_id]
        h50s = h50_by_id[peptide_id]
        first = mics[0]
        point_mics = [r for r in mics if r.get("concentration_parse_status") == "exact"
                      and r.get("unit_raw") == "µM"]
        unit_h50 = [r for r in h50s if r.get("activity_measure_value_raw", "").strip() == "50% Hemolysis"
                    and r.get("unit_raw") == "µM"]
        if not point_mics or not unit_h50:
            continue
        max_mic = max(float(r["concentration_value"]) for r in point_mics)
        h50 = unit_h50[0]
        h50_raw = h50.get("concentration_raw", "")
        censored = h50.get("concentration_parse_status") == "censored"
        h50_value = h50.get("concentration_value", "")
        ratio = ""
        ratio_interpretation = "not calculated"
        priority_rank = ""
        if h50_value:
            ratio_number = float(h50_value) / max_mic
            operator = h50.get("concentration_operator", "")
            if censored and operator in (">", ">="):
                ratio = f">{ratio_number:.2f}"
                ratio_interpretation = "lower bound on apparent concentration ratio"
            elif censored and operator in ("<", "<="):
                ratio = f"<{ratio_number:.2f}"
                ratio_interpretation = "upper bound on apparent concentration ratio"
            else:
                ratio = f"{ratio_number:.2f}"
                ratio_interpretation = "apparent concentration ratio; assays remain protocol-specific"
        measurements = [f"{r.get('target_species_raw','')}: {r.get('concentration_raw','')} {r.get('unit_raw','')}"
                        for r in mics]
        rows_out.append({
            "dbaasp_id": peptide_id,
            "peptide_name": first.get("name", "").strip(),
            "sequence": first.get("sequence", ""),
            "number_of_named_isolates": str(len(mics)),
            "named_isolate_mic_values_raw": " | ".join(measurements),
            "mic_unit_raw": "µM",
            "highest_mic_value": str(max_mic),
            "human_rbc_50_percent_hemolysis_raw": h50_raw,
            "h50_unit_raw": h50.get("unit_raw", ""),
            "apparent_h50_over_highest_mic": ratio,
            "priority_rank_within_published_panel": priority_rank,
            "ratio_interpretation": ratio_interpretation,
            "evidence_pmid": PMID,
            "interpretation_limit": "Same cited study and matching concentration units; not a clinical therapeutic index or proof of mechanism.",
        })

    report = {
        "source_pmid": PMID,
        "paper_title": "Antimicrobial properties of brevinin-2-related peptide and its analogs: Efficacy against multidrug-resistant Acinetobacter baumannii",
        "variants_with_named_isolate_mic_and_human_rbc_h50": len(rows_out),
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(REPORT.relative_to(ROOT))],
        "lead_interpretation": "This is a published benchmark/control panel, not novel peptide discovery. The ratio is a rough same-study concentration comparison; it is not a validated safety index.",
    }
    # Transparent sorting rule: descending conservative ratio, then ascending
    # worst named-isolate MIC. This ranks evidence within this one paper only.
    def sort_key(row: dict[str, str]) -> tuple[float, float]:
        ratio_value = float(row["apparent_h50_over_highest_mic"].lstrip("><="))
        return (-ratio_value, float(row["highest_mic_value"]))

    rows_out.sort(key=sort_key)
    for rank, row in enumerate(rows_out, start=1):
        row["priority_rank_within_published_panel"] = str(rank)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows_out)
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
