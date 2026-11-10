#!/usr/bin/env python3
"""Create a transparent, assay-linked Pareto screen for measured leads.

Generated hypotheses are retained in the table but stay explicitly unscored
until measurements or validated predictor outputs exist.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = ROOT / "results/candidate_prioritization/b2rp_published_lead_benchmark.csv"
HYPOTHESES = ROOT / "results/sequence_generation/b2rp_analogue_hypotheses.csv"
OUTPUT = ROOT / "results/candidate_prioritization/multiobjective_candidate_priority.csv"
REPORT = ROOT / "results/candidate_prioritization/multiobjective_candidate_priority.json"
CONCENTRATION = re.compile(r"^\s*([<>]=?)?\s*([0-9]+(?:\.[0-9]+)?)\s*$")


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def concentration_lower_bound(raw: str) -> tuple[float | None, str]:
    match = CONCENTRATION.fullmatch(raw or "")
    if not match:
        return None, "not_parseable"
    operator, value_text = match.groups()
    value = float(value_text)
    if operator in ("<", "<="):
        return None, "upper_censored_no_lower_bound"
    if operator in (">", ">="):
        return value, f"{operator}{value_text}_reported_lower_bound"
    return value, "exact_reported_value"


def dominates(left: dict, right: dict) -> bool:
    """Conservative point screen: lower MIC and higher reported H50 lower bound."""
    no_worse = (left["mic_value"] <= right["mic_value"] and
                left["h50_lower_bound"] >= right["h50_lower_bound"])
    strictly_better = (left["mic_value"] < right["mic_value"] or
                       left["h50_lower_bound"] > right["h50_lower_bound"])
    return no_worse and strictly_better


def assign_pareto_fronts(rows: list[dict]) -> None:
    remaining = [row for row in rows if row["evidence_status"] == "paired_experimental_measurements"]
    tier = 1
    while remaining:
        front = [candidate for candidate in remaining if not any(
            other is not candidate and dominates(other, candidate) for other in remaining
        )]
        if not front:
            raise RuntimeError("Pareto front assignment did not converge.")
        for candidate in front:
            candidate["pareto_front"] = tier
            candidate["priority_status"] = f"Measured Pareto front {tier}"
        front_ids = {id(row) for row in front}
        remaining = [row for row in remaining if id(row) not in front_ids]
        tier += 1


def main() -> None:
    measured: list[dict] = []
    for source in load_csv(BENCHMARK):
        h50, censor_status = concentration_lower_bound(source.get("human_rbc_50_percent_hemolysis_raw", ""))
        eligible = (
            source.get("mic_unit_raw") == source.get("h50_unit_raw") == "µM"
            and h50 is not None
        )
        measured.append({
            "candidate_id": f"DBAASP:{source['dbaasp_id']}",
            "candidate_name": source.get("peptide_name", ""),
            "sequence": source.get("sequence", ""),
            "candidate_class": "published_control",
            "evidence_status": "paired_experimental_measurements" if eligible else "ineligible_unit_or_censoring",
            "mic_value": float(source["highest_mic_value"]) if eligible else None,
            "h50_lower_bound": h50 if eligible else None,
            "activity_objective": "worst_named_isolate_MIC",
            "activity_value": source.get("highest_mic_value", ""),
            "activity_unit": source.get("mic_unit_raw", ""),
            "safety_objective": "human_RBC_50_percent_hemolysis_lower_bound",
            "safety_reported_value": source.get("human_rbc_50_percent_hemolysis_raw", ""),
            "safety_lower_bound": h50 if eligible else "",
            "safety_unit": source.get("h50_unit_raw", ""),
            "safety_censoring_status": censor_status,
            "pareto_front": "" if eligible else "",
            "priority_status": "Awaiting eligible matched measurements" if not eligible else "Pareto tier pending",
            "evidence_pmid": source.get("evidence_pmid", ""),
            "interpretation": (
                "Same-paper observed MIC and human-RBC H50 evidence; this Pareto screen is not a clinical therapeutic index."
                if eligible else "Not ranked because units or censoring do not support the defined comparison."
            ),
        })

    assign_pareto_fronts(measured)
    hypotheses = []
    for source in load_csv(HYPOTHESES):
        if source.get("candidate_status") != "generated_unvalidated_hypothesis":
            continue
        hypotheses.append({
            "candidate_id": source.get("candidate_label", ""),
            "candidate_name": source.get("candidate_label", ""),
            "sequence": source.get("sequence", ""),
            "candidate_class": source.get("candidate_status", "generated_hypothesis"),
            "evidence_status": "no_paired_experimental_measurements",
            "activity_objective": "worst_named_isolate_MIC",
            "activity_value": "",
            "activity_unit": "",
            "safety_objective": "human_RBC_50_percent_hemolysis_lower_bound",
            "safety_reported_value": "",
            "safety_lower_bound": "",
            "safety_unit": "",
            "safety_censoring_status": "not_measured",
            "pareto_front": "",
            "priority_status": "Unscored hypothesis; experimental evidence required",
            "evidence_pmid": "",
            "interpretation": "No activity, toxicity, hemolysis, or docking value is imputed for this sequence.",
        })

    rows = measured + hypotheses
    for row in rows:
        row.pop("mic_value", None)
        row.pop("h50_lower_bound", None)
    fields = [
        "candidate_id", "candidate_name", "sequence", "candidate_class", "evidence_status",
        "activity_objective", "activity_value", "activity_unit", "safety_objective",
        "safety_reported_value", "safety_lower_bound", "safety_unit", "safety_censoring_status",
        "pareto_front", "priority_status", "evidence_pmid", "interpretation",
    ]
    rows.sort(key=lambda row: (
        0 if row["evidence_status"] == "paired_experimental_measurements" else 1,
        int(row["pareto_front"]) if row["pareto_front"] else 999,
        row["candidate_name"].casefold(),
    ))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "method": "Pareto tiering of paired same-paper observed activity and human-RBC hemolysis evidence",
        "activity_objective": "lower worst MIC across five named NM isolates is better; µM only",
        "safety_objective": "higher reported human-RBC 50% hemolysis concentration is better; µM only",
        "censoring_rule": "For right-censored H50 (e.g. >200 µM), use the reported bound as a conservative value and preserve the censoring label; do not treat it as the exact H50.",
        "dominance_rule": "A measured candidate is no worse on both objectives and strictly better on at least one.",
        "unscored_rule": "Generated hypotheses without paired experimental measurements remain unscored; no model or docking value is imputed.",
        "measured_paired_candidates": sum(r["evidence_status"] == "paired_experimental_measurements" for r in rows),
        "generated_hypotheses_unscored": sum(r["evidence_status"] == "no_paired_experimental_measurements" for r in rows),
        "output": str(OUTPUT.relative_to(ROOT)),
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
