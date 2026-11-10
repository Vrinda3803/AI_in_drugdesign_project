#!/usr/bin/env python3
"""Extract DBAASP A. baumannii MIC records linked to reviewed MDR/CRAB papers."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data/processed/dbaasp_a_baumannii_mic_parsed.csv"
OUTPUT = ROOT / "data/processed/dbaasp_a_baumannii_reviewed_paper_records.csv"

# Paper-level claim verified from PubMed abstract; this does not imply that
# every single row is mapped to a particular resistant isolate.
REVIEWED = {
    "19793185": ("MDR", "PubMed abstract reports B2RP and analogues active against nosocomial MDR A. baumannii isolates."),
    "24071034": ("MDR", "PubMed abstract reports LL-37 fragments tested against MDR A. baumannii clinical isolates."),
    "40056800": ("CRAB", "PubMed abstract reports cathelicidin-derived peptides tested against CRAB clinical isolates."),
}


def main() -> None:
    with INPUT.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"No input records in {INPUT}")

    selected = []
    for row in rows:
        try:
            pmids = json.loads(row.get("candidate_pubmed_ids_json", "[]") or "[]")
        except json.JSONDecodeError:
            continue
        for pmid in pmids:
            pmid = str(pmid)
            if pmid not in REVIEWED:
                continue
            phenotype, basis = REVIEWED[pmid]
            out = dict(row)
            out["reviewed_pmid"] = pmid
            out["paper_level_resistance_phenotype"] = phenotype
            out["paper_abstract_evidence"] = basis
            out["row_level_isolate_assignment_reviewed"] = "no"
            selected.append(out)

    fields = list(rows[0]) + [
        "reviewed_pmid", "paper_level_resistance_phenotype",
        "paper_abstract_evidence", "row_level_isolate_assignment_reviewed",
    ]
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)
    print(f"Wrote {len(selected)} paper-linked MIC records to {OUTPUT.relative_to(ROOT)}")
    for pmid, (phenotype, _) in REVIEWED.items():
        matched = [r for r in selected if r["reviewed_pmid"] == pmid]
        print(f"{pmid} ({phenotype}): {len(matched)} rows")


if __name__ == "__main__":
    main()
