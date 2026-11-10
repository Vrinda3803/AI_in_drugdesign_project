#!/usr/bin/env python3
"""Prepare a source-review cohort of reported human-RBC 50% hemolysis doses."""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/preprocessing"))
from parse_mic_observations import parse  # shared conservative numeric syntax parser

SOURCE = ROOT / "data/processed/dbaasp_a_baumannii_toxicity_records.csv"
OUTPUT = ROOT / "data/processed/dbaasp_human_rbc_h50_candidate_records.csv"
REPORT = ROOT / "data/processed/human_rbc_h50_candidate_report.json"
STANDARD_AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def main() -> None:
    with SOURCE.open("r", encoding="utf-8-sig", newline="") as handle:
        all_rows = list(csv.DictReader(handle))
    selected = [
        row for row in all_rows
        if row.get("target_cell_name_raw", "").strip().casefold() == "human erythrocytes"
        and row.get("activity_measure_value_raw", "").strip().casefold() == "50% hemolysis"
        and "hemolysis" in row.get("activity_measure_group_raw", "").casefold()
    ]

    output_rows: list[dict[str, str]] = []
    statuses: Counter[str] = Counter()
    units: Counter[str] = Counter()
    seqs: set[str] = set()
    valid_seqs: set[str] = set()
    for row in selected:
        parsed = parse(row.get("concentration_raw", ""))
        statuses[parsed["concentration_parse_status"]] += 1
        units[row.get("unit_raw") or "(missing)"] += 1
        seq = row.get("sequence", "").strip().upper()
        valid = bool(STANDARD_AA.fullmatch(seq))
        if seq:
            seqs.add(seq)
        if valid:
            valid_seqs.add(seq)
        out = dict(row)
        out.update(parsed)
        out["sequence_qc"] = "valid_standard_aa" if valid else "nonstandard_or_invalid"
        # Candidate for review only: not a validated biological safety label.
        out["endpoint_interpretation"] = "reported concentration at 50% human erythrocyte hemolysis; assay-specific"
        out["cohort_review_status"] = "candidate_requires_paper_and_assay_review"
        output_rows.append(out)

    fields = list(output_rows[0]) if output_rows else []
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)

    report = {
        "selection_rule": "Human erythrocytes; activityMeasureForLysisValue exactly 50% Hemolysis; measurement group contains Hemolysis.",
        "source_rows_scanned": len(all_rows),
        "candidate_observation_rows": len(selected),
        "unique_nonempty_sequences": len(seqs),
        "unique_standard_aa_sequences": len(valid_seqs),
        "concentration_parse_status_counts": dict(statuses),
        "unit_raw_counts": dict(units),
        "output": str(OUTPUT.relative_to(ROOT)),
        "interpretation": "A candidate assay-specific H50 cohort for source review. Rows do not define a universal toxicity threshold and are not labeled safe/unsafe.",
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
