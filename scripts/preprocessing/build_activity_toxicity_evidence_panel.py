#!/usr/bin/env python3
"""Join DBAASP MIC and human-RBC H50 observations by peptide ID, without scoring."""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MIC_PATH = ROOT / "data/processed/dbaasp_a_baumannii_mic_parsed.csv"
H50_PATH = ROOT / "data/processed/dbaasp_human_rbc_h50_candidate_records.csv"
OUTPUT = ROOT / "data/processed/dbaasp_activity_toxicity_evidence_panel.csv"
REPORT = ROOT / "data/processed/activity_toxicity_overlap_report.json"

FIELDS = [
    "dbaasp_id", "dbaasp_accession", "peptide_name", "sequence", "sequence_qc",
    "mic_observation_count", "h50_observation_count", "has_exact_standard_aa_mic",
    "has_exact_standard_aa_h50", "mic_observations_json", "h50_observations_json",
]


def main() -> None:
    with MIC_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        mic_rows = [row for row in csv.DictReader(handle) if row.get("source") == "DBAASP"]
    with H50_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        h50_rows = list(csv.DictReader(handle))

    mic_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    h50_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mic_rows:
        mic_by_id[row["dbaasp_id"]].append(row)
    for row in h50_rows:
        h50_by_id[row["dbaasp_id"]].append(row)

    ids = sorted(mic_by_id, key=lambda value: int(value))
    out_rows: list[dict[str, str]] = []
    shared_units: Counter[str] = Counter()
    for peptide_id in ids:
        mics = mic_by_id[peptide_id]
        h50s = h50_by_id.get(peptide_id, [])
        first = mics[0]
        exact_mic = [r for r in mics if r.get("concentration_parse_status") == "exact"
                     and r.get("sequence_qc") == "valid_standard_aa"]
        exact_h50 = [r for r in h50s if r.get("concentration_parse_status") == "exact"
                     and r.get("sequence_qc") == "valid_standard_aa"]
        for mic in exact_mic:
            for h50 in exact_h50:
                if mic.get("unit_raw") and mic.get("unit_raw") == h50.get("unit_raw"):
                    shared_units[mic["unit_raw"]] += 1
        out_rows.append({
            "dbaasp_id": peptide_id,
            "dbaasp_accession": first.get("dbaasp_accession", ""),
            "peptide_name": first.get("name", ""),
            "sequence": first.get("sequence", ""),
            "sequence_qc": first.get("sequence_qc", ""),
            "mic_observation_count": str(len(mics)),
            "h50_observation_count": str(len(h50s)),
            "has_exact_standard_aa_mic": str(bool(exact_mic)).lower(),
            "has_exact_standard_aa_h50": str(bool(exact_h50)).lower(),
            # Retain all original columns, units, strain and reference context.
            "mic_observations_json": json.dumps(mics, ensure_ascii=False, separators=(",", ":")),
            "h50_observations_json": json.dumps(h50s, ensure_ascii=False, separators=(",", ":")),
        })

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(out_rows)

    overlap_ids = set(mic_by_id) & set(h50_by_id)
    exact_both_ids = {
        peptide_id for peptide_id in overlap_ids
        if any(r.get("concentration_parse_status") == "exact" and r.get("sequence_qc") == "valid_standard_aa"
               for r in mic_by_id[peptide_id])
        and any(r.get("concentration_parse_status") == "exact" and r.get("sequence_qc") == "valid_standard_aa"
                for r in h50_by_id[peptide_id])
    }
    report = {
        "mic_peptide_ids": len(mic_by_id),
        "human_rbc_h50_candidate_peptide_ids": len(h50_by_id),
        "peptide_ids_with_both_evidence_types": len(overlap_ids),
        "overlap_with_exact_standard_aa_mic_and_h50": len(exact_both_ids),
        "matched_exact_mic_h50_pairs_with_same_raw_unit": dict(shared_units),
        "output_rows": len(out_rows),
        "join_key": "DBAASP peptide record ID; records are not joined on sequence alone",
        "interpretation": "Evidence inventory only. Raw observations remain grouped by peptide but are not aggregated or converted into a therapeutic-index or safety score.",
        "output": str(OUTPUT.relative_to(ROOT)),
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
