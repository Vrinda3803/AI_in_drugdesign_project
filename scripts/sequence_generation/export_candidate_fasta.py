"""Export the B2RP control and hypothesis panel for external peptide tools."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results" / "sequence_generation" / "b2rp_analogue_hypotheses.csv"
LEAD = ROOT / "results" / "candidate_prioritization" / "b2rp_published_lead_benchmark.csv"
OUT = ROOT / "results" / "candidate_prioritization" / "b2rp_candidate_panel.fasta"


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if LEAD.exists():
        with LEAD.open(encoding="utf-8-sig", newline="") as handle:
            lead = next(csv.DictReader(handle), None)
        if lead:
            rows.append({"candidate_label": "B2RP_top_published_lead_control",
                         "candidate_status": "published_positive_control",
                         "sequence": lead["sequence"]})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="ascii", newline="\n") as handle:
        for row in rows:
            handle.write(f">{row['candidate_label']}|{row['candidate_status']}\n")
            seq = row["sequence"]
            for start in range(0, len(seq), 60):
                handle.write(seq[start:start + 60] + "\n")
    print(f"Exported {len(rows)} sequences to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
