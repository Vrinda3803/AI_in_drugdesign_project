"""Validate and import official ToxinPred3/HemoPI2 CSV results for the UI.

The tool result file is never submitted anywhere by this script. It reads a
user-saved result CSV, maps tool IDs to the project's candidate FASTA, verifies
the returned sequence, and appends normalized provenance rows for the UI.
"""
from __future__ import annotations

import argparse
import csv
import re
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "results/candidate_prioritization/b2rp_candidate_panel.fasta"
OUT = ROOT / "results/candidate_prioritization/external_safety_predictions.csv"
FIELDS = ["candidate_id", "sequence", "predictor", "model", "endpoint", "prediction",
          "score_name", "score", "score_unit", "ppv", "threshold", "run_date"]


def normalize(value: str) -> str:
    value = value.casefold().replace("μ", "u").replace("µ", "u")
    return re.sub(r"[^a-z0-9]", "", value)


def read_panel() -> dict[str, str]:
    mapping: dict[str, str] = {}
    current = None
    with PANEL.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split("|", 1)[0]
                if current in mapping:
                    raise ValueError(f"Duplicate candidate ID in FASTA: {current}")
                mapping[current] = ""
            elif current:
                mapping[current] += line.upper()
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", choices=["toxinpred3", "hemopi2"], required=True)
    parser.add_argument("--endpoint", choices=["toxicity", "hemolysis-classification", "hemolysis-hc50"], required=True)
    parser.add_argument("--input", type=Path, required=True, help="Saved CSV downloaded/exported from the selected predictor")
    parser.add_argument("--model", required=True, help="Exact model name selected in the predictor")
    parser.add_argument("--threshold", default="", help="Selected decision threshold, if applicable")
    parser.add_argument("--run-date", default=date.today().isoformat())
    args = parser.parse_args()

    expected_endpoint_tool = {"toxicity": "toxinpred3", "hemolysis-classification": "hemopi2", "hemolysis-hc50": "hemopi2"}
    if expected_endpoint_tool[args.endpoint] != args.tool:
        raise SystemExit(f"Endpoint {args.endpoint} requires tool {expected_endpoint_tool[args.endpoint]}.")
    panel = read_panel()
    records: list[dict[str, str]] = []
    with args.input.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise SystemExit("Input CSV has no header row.")
        columns = {normalize(name): name for name in reader.fieldnames}

        def field(row: dict[str, str], *aliases: str) -> str:
            for alias in aliases:
                name = columns.get(normalize(alias))
                if name and (row.get(name) or "").strip():
                    return row[name].strip()
            return ""

        def matched_field_name(row: dict[str, str], *aliases: str) -> str:
            for alias in aliases:
                name = columns.get(normalize(alias))
                if name and (row.get(name) or "").strip():
                    return name
            return ""

        if args.endpoint == "toxicity":
            score_aliases = ("Hybrid Score", "ML Score", "Score") if "hybrid" in args.model.casefold() else ("ML Score", "Hybrid Score", "Score")
            predictor, endpoint = "ToxinPred3", "peptide toxicity classification"
            score_unit = "model score"
        elif args.endpoint == "hemolysis-classification":
            score_aliases = ("Hybrid Score", "ESM Score", "ML Score", "Score") if "hybrid" in args.model.casefold() else ("ESM Score", "ML Score", "Hybrid Score", "Score")
            predictor, endpoint = "HemoPI2", "erythrocyte hemolysis classification"
            score_unit = "model score"
        else:
            score_aliases = ("HC50(μM)", "HC50(µM)", "HC50 (μM)", "HC50", "Score")
            predictor, endpoint = "HemoPI2", "predicted hemolysis concentration"
            score_unit = "µM"

        bad_rows = []
        for line_number, row in enumerate(reader, 2):
            raw_id = field(row, "SeqID", "ID", "Subject", "Name")
            candidate_id = raw_id.lstrip("> ").split("|", 1)[0]
            sequence = field(row, "Sequence", "Seq").upper()
            if candidate_id not in panel:
                bad_rows.append(f"line {line_number}: unrecognized candidate ID {raw_id!r}")
                continue
            if not sequence:
                bad_rows.append(f"line {line_number}: no returned sequence for {candidate_id}")
                continue
            if sequence != panel[candidate_id]:
                bad_rows.append(f"line {line_number}: sequence mismatch for {candidate_id}")
                continue
            prediction = field(row, "Prediction")
            score = field(row, *score_aliases)
            score_name = matched_field_name(row, *score_aliases)
            if args.endpoint == "hemolysis-hc50" and score:
                score_name = "HC50"
            if not prediction and not score:
                bad_rows.append(f"line {line_number}: no recognized prediction or score for {candidate_id}")
                continue
            records.append({"candidate_id": candidate_id, "sequence": sequence, "predictor": predictor,
                            "model": args.model, "endpoint": endpoint, "prediction": prediction,
                            "score_name": score_name, "score": score, "score_unit": score_unit,
                            "ppv": field(row, "PPV") if args.endpoint == "toxicity" else "",
                            "threshold": args.threshold,
                            "run_date": args.run_date})

    if bad_rows:
        raise SystemExit("Import stopped; no rows written. Resolve these issues first:\n" + "\n".join(bad_rows))
    if not records:
        raise SystemExit("No result rows were imported; check the selected CSV and model output.")

    existing = []
    if OUT.exists():
        with OUT.open(encoding="utf-8-sig", newline="") as handle:
            existing = list(csv.DictReader(handle))
    key = lambda row: (row["candidate_id"], row["predictor"], row["model"], row["endpoint"], row["run_date"])
    incoming = {key(row): row for row in records}
    merged = [row for row in existing if key(row) not in incoming] + list(incoming.values())
    with OUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)
    print(f"Imported {len(records)} rows into {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
