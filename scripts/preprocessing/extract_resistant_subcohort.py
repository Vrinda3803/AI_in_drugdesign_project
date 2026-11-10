#!/usr/bin/env python3
"""Extract records whose DBAASP isolate note explicitly signals drug resistance.

This is a provenance-preserving curation aid, not a new activity-labeling rule:
the script selects text explicitly present in the isolate field and retains the
original MIC and source fields. Publication titles alone are deliberately not
used to label a record as MDR.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data/processed/dbaasp_a_baumannii_mic_parsed.csv"
OUTPUT = ROOT / "data/processed/dbaasp_a_baumannii_resistant_isolate_records.csv"
REPORT = ROOT / "data/processed/resistant_isolate_subcohort_report.json"

# Explicit terms only. Preserve the source wording; do not infer resistance from
# peptide names, an article title, or the species name.
RESISTANCE_NOTE = re.compile(
    r"\b(?:MDR|XDR|PDR|CRAB|carbapenem[- ]resistant|"
    r"multidrug[- ]resistant|extensively drug[- ]resistant)\b",
    re.IGNORECASE,
)


def main() -> None:
    with INPUT.open("r", encoding="utf-8-sig", newline="") as handle:
        source_rows = list(csv.DictReader(handle))
    selected = [
        row for row in source_rows
        if RESISTANCE_NOTE.search(row.get("strain_or_isolate_note_raw", "") or "")
    ]

    fields = list(source_rows[0]) if source_rows else []
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)

    report = {
        "input": str(INPUT.relative_to(ROOT)),
        "output": str(OUTPUT.relative_to(ROOT)),
        "selection_rule": (
            "strain_or_isolate_note_raw explicitly matches MDR, XDR, PDR, CRAB, "
            "carbapenem-resistant, multidrug-resistant, or extensively drug-resistant"
        ),
        "publication_title_used_for_selection": False,
        "source_rows_total": len(source_rows),
        "selected_activity_rows": len(selected),
        "selected_mic_rows": sum(row.get("activity_measure") == "MIC" for row in selected),
        "selected_exact_point_mic_rows": sum(
            row.get("activity_measure") == "MIC"
            and row.get("concentration_parse_status") == "exact"
            for row in selected
        ),
        "unique_nonempty_sequences": len({row["sequence"] for row in selected if row.get("sequence")}),
        "unique_valid_standard_aa_sequences": len({
            row["sequence"] for row in selected
            if row.get("sequence_qc") == "valid_standard_aa"
        }),
        "unit_raw_counts": dict(Counter(row.get("unit_raw", "") for row in selected)),
        "interpretation": (
            "An explicitly labeled resistant-isolate subcohort for audit and "
            "future stratified evaluation. It is not a standardized MDR-only "
            "benchmark; resistance terms and assay conditions require source-paper review."
        ),
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
