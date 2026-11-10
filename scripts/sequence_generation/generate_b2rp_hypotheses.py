#!/usr/bin/env python3
"""Generate auditable B2RP analog hypotheses from published single-site SAR."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FEATURES = ROOT / "data/processed/sequence_features.csv"
OUTPUT_DIR = ROOT / "results/sequence_generation"
OUTPUT = OUTPUT_DIR / "b2rp_analogue_hypotheses.csv"
REPORT = OUTPUT_DIR / "b2rp_analogue_generation_report.json"
SEED = "GIWDTIKSMGKVFAGKILQNL"
SEED_NAME = "B2RP parent; source paper reports C-terminal amide"

VARIANTS = [
    {"label": "B2RP_parent", "changes": {}, "status": "reference_parent",
     "rationale": "Parent control from the source study."},
    {"label": "D4K", "changes": {4: "K"}, "status": "published_control",
     "rationale": "The source paper reports increased activity against MDR A. baumannii isolates."},
    {"label": "L18K", "changes": {18: "K"}, "status": "published_control",
     "rationale": "The source paper reports reduced potency and reduced hemolysis; retain as an SAR control."},
    {"label": "D4K_L18K", "changes": {4: "K", 18: "K"}, "status": "published_control",
     "rationale": "The source paper reports retained activity with low hemolytic activity for this combination."},
    {"label": "D4R_hypothesis", "changes": {4: "R"}, "status": "generated_unvalidated_hypothesis",
     "rationale": "Arginine is a conservative basic substitution at the source-supported D4K position; activity is unknown."},
    {"label": "L18R_hypothesis", "changes": {18: "R"}, "status": "generated_unvalidated_hypothesis",
     "rationale": "Arginine tests a basic replacement at the source-supported L18K position; activity and hemolysis are unknown."},
    {"label": "D4K_L18R_hypothesis", "changes": {4: "K", 18: "R"}, "status": "generated_unvalidated_hypothesis",
     "rationale": "Combines the source-supported D4K change with an untested basic alternative at position 18; no activity prediction."},
    {"label": "D4R_L18K_hypothesis", "changes": {4: "R", 18: "K"}, "status": "generated_unvalidated_hypothesis",
     "rationale": "Tests an unvalidated basic alternative at position 4 with the source-supported L18K change."},
    {"label": "D4R_L18R_hypothesis", "changes": {4: "R", 18: "R"}, "status": "generated_unvalidated_hypothesis",
     "rationale": "Tests basic substitutions at both source-informed positions; activity and hemolysis are unknown."},
]


def levenshtein(a: str, b: str) -> int:
    previous = list(range(len(b) + 1))
    for i, char_a in enumerate(a, start=1):
        current = [i]
        for j, char_b in enumerate(b, start=1):
            current.append(min(current[-1] + 1, previous[j] + 1,
                               previous[j - 1] + (char_a != char_b)))
        previous = current
    return previous[-1]


def main() -> None:
    with FEATURES.open("r", encoding="utf-8-sig", newline="") as handle:
        known_rows = list(csv.DictReader(handle))
    known_sequences = {row["sequence"]: row for row in known_rows if row.get("sequence")}
    output_rows = []
    for variant in VARIANTS:
        sequence = list(SEED)
        mutation_text = []
        for position, target in sorted(variant["changes"].items()):
            original = sequence[position - 1]
            sequence[position - 1] = target
            mutation_text.append(f"{original}{position}{target}")
        sequence = "".join(sequence)
        exact = known_sequences.get(sequence)
        nearest_sequence, nearest_row = min(
            known_sequences.items(),
            key=lambda item: (levenshtein(sequence, item[0]), item[0]),
        )
        distance = levenshtein(sequence, nearest_sequence)
        output_rows.append({
            "candidate_label": variant["label"],
            "sequence": sequence,
            "length": str(len(sequence)),
            "mutations_from_B2RP": ";".join(mutation_text) or "none",
            "candidate_status": variant["status"],
            "local_exact_match": "yes" if exact else "no",
            "local_exact_match_sources": exact.get("sources_and_ids_json", "") if exact else "",
            "nearest_local_sequence": nearest_sequence,
            "nearest_local_edit_distance": str(distance),
            "nearest_local_edit_identity": f"{(max(len(sequence), len(nearest_sequence)) - distance) / max(len(sequence), len(nearest_sequence)):.4f}",
            "nearest_local_sources": nearest_row.get("sources_and_ids_json", ""),
            "terminal_modification_assumption": "C-terminal amide copied from published B2RP parent; confirm before synthesis.",
            "design_rationale": variant["rationale"],
            "predicted_activity": "not predicted",
            "predicted_toxicity": "not predicted",
        })

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    novel_local = [r for r in output_rows if r["candidate_status"] == "generated_unvalidated_hypothesis"
                   and r["local_exact_match"] == "no"]
    report = {
        "seed_sequence": SEED,
        "seed_name": SEED_NAME,
        "source_sar_pmid": "19793185",
        "known_sequence_panel_size": len(known_sequences),
        "candidate_sequences_total_including_controls": len(output_rows),
        "generated_hypotheses": sum(r["candidate_status"] == "generated_unvalidated_hypothesis" for r in output_rows),
        "hypotheses_not_found_in_local_panel": len(novel_local),
        "hypothesis_labels": [r["candidate_label"] for r in novel_local],
        "design_rule": "Use only source-supported positions 4 and 18; explore Arg in place of Lys as an explicit, unvalidated basic-residue hypothesis. Do not use the weak MIC regressor to score candidates.",
        "novelty_scope": "Exact sequence absence from the local 981-sequence feature panel only; no global literature/patent novelty claim.",
        "prediction_scope": "No MIC, hemolysis, safety, uptake, target-binding, or efficacy predictions are made.",
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(REPORT.relative_to(ROOT))],
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
