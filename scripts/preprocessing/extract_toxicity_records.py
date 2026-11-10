#!/usr/bin/env python3
"""Extract DBAASP hemolysis/cytotoxicity observations without pooling assays."""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "data/raw/DBAASP/a_baumannii_details"
OUTPUT = ROOT / "data/processed/dbaasp_a_baumannii_toxicity_records.csv"
REPORT = ROOT / "data/processed/toxicity_data_report.json"

FIELDS = [
    "source", "dbaasp_id", "dbaasp_accession", "peptide_name", "sequence",
    "sequence_length", "n_terminus_json", "c_terminus_json", "unusual_amino_acids_json",
    "toxicity_record_id", "target_cell_name_raw", "target_cell_description_raw",
    "activity_measure_group_raw", "activity_measure_value_raw", "activity_numeric_raw",
    "concentration_raw", "unit_raw", "ionic_strength_raw", "salt_type_raw", "ph_raw",
    "note_raw", "reference_raw", "reference_resolution_status",
    "reference_article_index_candidate", "candidate_pubmed_id", "candidate_reference_title",
    "detail_file",
]


def candidate_article(detail: dict, ref: str) -> tuple[str, str, str, str]:
    """Resolve a numeric article-array index as a candidate, retaining uncertainty."""
    try:
        idx = int(ref)
    except (TypeError, ValueError):
        return "", "unresolved_non_numeric_reference", "", ""
    articles = detail.get("articles") or []
    if idx < 1 or idx > len(articles):
        return str(idx), "unresolved_out_of_range_reference", "", ""
    article = articles[idx - 1]
    pubmed = article.get("pubmed") or {}
    status = "candidate_1_based_article_index_mapping"
    return str(idx), status, str(pubmed.get("pubmedId") or ""), str(article.get("title") or "")


def main() -> None:
    detail_files = sorted(INPUT_DIR.glob("*.json"), key=lambda path: int(path.stem))
    rows: list[dict[str, str]] = []
    measure_counts: Counter[str] = Counter()
    cell_counts: Counter[str] = Counter()
    unit_counts: Counter[str] = Counter()
    peptide_ids: set[str] = set()

    for path in detail_files:
        detail = json.loads(path.read_text(encoding="utf-8"))
        dbaasp_id = str(detail.get("id", path.stem))
        peptide_ids.add(dbaasp_id)
        for activity in detail.get("hemoliticCytotoxicActivities") or []:
            group = (activity.get("activityMeasureForLysisGroup") or {}).get("name", "")
            cell = activity.get("targetCell") or {}
            unit = (activity.get("unit") or {}).get("name", "")
            ref = str(activity.get("reference") or "")
            index, status, pubmed, title = candidate_article(detail, ref)
            measure_counts[group or "(missing)"] += 1
            cell_counts[cell.get("name") or "(missing)"] += 1
            unit_counts[unit or "(missing)"] += 1
            rows.append({
                "source": "DBAASP",
                "dbaasp_id": dbaasp_id,
                "dbaasp_accession": str(detail.get("dbaaspId") or ""),
                "peptide_name": str(detail.get("name") or ""),
                "sequence": str(detail.get("sequence") or ""),
                "sequence_length": str(detail.get("sequenceLength") or ""),
                "n_terminus_json": json.dumps(detail.get("nTerminus") or detail.get("NTerminus"), ensure_ascii=False),
                "c_terminus_json": json.dumps(detail.get("cTerminus") or detail.get("CTerminus"), ensure_ascii=False),
                "unusual_amino_acids_json": json.dumps(detail.get("unusualAminoAcids") or [], ensure_ascii=False),
                "toxicity_record_id": str(activity.get("id") or ""),
                "target_cell_name_raw": str(cell.get("name") or ""),
                "target_cell_description_raw": str(cell.get("description") or ""),
                "activity_measure_group_raw": str(group),
                "activity_measure_value_raw": str(activity.get("activityMeasureForLysisValue") or ""),
                "activity_numeric_raw": str(activity.get("activity") or ""),
                "concentration_raw": str(activity.get("concentration") or ""),
                "unit_raw": str(unit),
                "ionic_strength_raw": str(activity.get("ionicStrength") or ""),
                "salt_type_raw": str(activity.get("saltType") or ""),
                "ph_raw": str(activity.get("ph") or ""),
                "note_raw": str(activity.get("note") or ""),
                "reference_raw": ref,
                "reference_resolution_status": status,
                "reference_article_index_candidate": index,
                "candidate_pubmed_id": pubmed,
                "candidate_reference_title": title,
                "detail_file": str(path.relative_to(ROOT)),
            })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "source": "DBAASP v4 detail responses downloaded for the A. baumannii target search",
        "detail_records_scanned": len(detail_files),
        "unique_dbaasp_peptide_ids_scanned": len(peptide_ids),
        "toxicity_observation_rows": len(rows),
        "unique_peptide_ids_with_toxicity_observations": len({r["dbaasp_id"] for r in rows}),
        "measurement_group_counts": dict(measure_counts),
        "target_cell_counts": dict(cell_counts.most_common()),
        "unit_raw_counts": dict(unit_counts),
        "reference_mapping_note": "Article index to PubMed mapping is retained as candidate pending record-level/source-paper validation.",
        "interpretation_note": "Hemolysis and cytotoxicity endpoints, cell types, concentration values, and units are not pooled or converted to binary safety labels.",
        "output": str(OUTPUT.relative_to(ROOT)),
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
