#!/usr/bin/env python3
"""Build the static browser UI's local peptide lookup JSON."""
from __future__ import annotations

import csv
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UI_DATA = ROOT / "ui/data/peptides.json"
UI_DATA_JS = ROOT / "ui/data/peptides.js"
DETAILS = ROOT / "data/raw/DBAASP/a_baumannii_details"
ABAMP = ROOT / "data/processed/abampdb_records_qc.csv"
GENERATED = ROOT / "results/sequence_generation/b2rp_analogue_hypotheses.csv"
PREDICTIONS = ROOT / "results/candidate_prioritization/external_safety_predictions.csv"
PANEL = ROOT / "results/candidate_prioritization/b2rp_candidate_panel.fasta"
DOCKING_RUNS = ROOT / "results/docking/cabsdock/observed_runs.csv"
DOCKING_WEB_JOBS = ROOT / "results/docking/cabsdock/web_jobs.csv"


def read_fasta_sequences(path: Path) -> set[str]:
    sequences: set[str] = set()
    current: list[str] = []
    if not path.exists():
        return sequences
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    sequences.add("".join(current))
                current = []
            elif line:
                current.append(line)
    if current:
        sequences.add("".join(current))
    return sequences


def docking_fields(
    sequence: str,
    planned_sequences: set[str],
    observed_runs: dict[str, list[dict[str, str]]],
    web_jobs: dict[str, dict[str, str]],
) -> dict[str, str | None]:
    runs = observed_runs.get(sequence, [])
    completed = [row for row in runs if row.get("status") == "complete"]
    partial = [row for row in runs if row.get("status") == "partial"]
    if completed:
        receptor_states = sorted({Path(row["receptor_file"]).stem for row in completed})
        models = sum(int(row.get("ranked_pose_model_count") or 0) for row in completed)
        return {
            "docking_score": None,
            "docking_status": "Pose outputs available",
            "docking_note": (
                f"Ranked pose models are inventoried for {len(completed)} completed run(s) "
                f"across {', '.join(receptor_states)} ({models} model files). "
                "No numeric affinity score is reported. Docking uses one-letter sequence and "
                "does not encode terminal modifications."
            ),
        }
    if partial:
        return {
            "docking_score": None,
            "docking_status": "Partial outputs; review logs",
            "docking_note": (
                f"Artifacts are present for {len(partial)} incomplete run(s); no completed docking "
                "result or numeric score is available. Inspect run logs and output files."
            ),
        }
    web_job = web_jobs.get(sequence)
    if web_job and web_job.get("status") == "completed_results_available":
        density = web_job.get("top_cluster_density") or ""
        density_note = (
            f" The top cluster density reported by CABS-dock is {density}; "
            "this is a method-specific cluster metric, not a docking affinity or binding free energy."
            if density else
            " Cluster details and models are available on the linked CABS-dock results page; "
            "no single docking affinity score is reported."
        )
        return {
            "docking_score": None,
            "docking_cluster_density": density or None,
            "docking_status": "CABS-dock models and clusters available",
            "docking_note": (
                f"CABS-dock returned {web_job.get('returned_model_count') or 'representative'} "
                f"representative models and {web_job.get('returned_cluster_count') or 'cluster'} "
                "clusters for apo FtsZ 22HR chain A (50 cycles)." + density_note +
                " A predicted pose does not establish target binding or antibacterial activity."
            ),
            "docking_job_url": web_job.get("job_url", ""),
            "docking_job_id": web_job.get("candidate_id", ""),
        }
    if web_job:
        return {
            "docking_score": None,
            "docking_cluster_density": None,
            "docking_status": "CABS-dock pending or running",
            "docking_note": (
                "The peptide-protein docking job was accepted by CABS-dock for apo FtsZ 22HR chain A "
                "(50 simulation cycles). Results are still pending. No score is available yet; "
                "cluster or energy outputs are not binding free energies."
            ),
            "docking_job_url": web_job.get("job_url", ""),
            "docking_job_id": web_job.get("candidate_id", ""),
        }
    if sequence in planned_sequences:
        return {
            "docking_score": None,
            "docking_status": "Sequence in plan; not run",
            "docking_note": (
                "This exact one-letter sequence is in the prepared 60-run CABSdock panel, "
                "but no poses or scores exist yet. Terminal modifications are not encoded."
            ),
        }
    return {
        "docking_score": None,
        "docking_status": "Not in plan; not run",
        "docking_note": (
            "Only the ten-sequence candidate panel has a prepared docking plan. "
            "This sequence has no docking run or score."
        ),
    }


def article_candidate(detail: dict, ref: str) -> dict[str, str]:
    try:
        index = int(ref)
    except (TypeError, ValueError):
        return {"status": "unresolved_reference", "pubmed": "", "title": ""}
    articles = detail.get("articles") or []
    if index < 1 or index > len(articles):
        return {"status": "unresolved_reference", "pubmed": "", "title": ""}
    article = articles[index - 1]
    pubmed = article.get("pubmed") or {}
    return {
        "status": "candidate_article_index_mapping",
        "pubmed": str(pubmed.get("pubmedId") or ""),
        "title": str(article.get("title") or ""),
    }


def main() -> None:
    planned_sequences = read_fasta_sequences(PANEL)
    observed_runs: dict[str, list[dict[str, str]]] = {}
    if DOCKING_RUNS.is_file():
        with DOCKING_RUNS.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                sequence = (row.get("sequence") or "").strip()
                if sequence:
                    observed_runs.setdefault(sequence, []).append(row)
    web_jobs: dict[str, dict[str, str]] = {}
    if DOCKING_WEB_JOBS.is_file():
        with DOCKING_WEB_JOBS.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                sequence = (row.get("sequence") or "").strip()
                if sequence:
                    web_jobs[sequence] = row
    predictions_by_id: dict[str, list[dict[str, str]]] = {}
    if PREDICTIONS.exists():
        with PREDICTIONS.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                candidate_id = (row.get("candidate_id") or "").strip()
                if candidate_id:
                    predictions_by_id.setdefault(candidate_id, []).append({
                        "sequence": (row.get("sequence") or "").strip(),
                        "predictor": (row.get("predictor") or "Not specified").strip(),
                        "model": (row.get("model") or "Not specified").strip(),
                        "endpoint": (row.get("endpoint") or "Not specified").strip(),
                        "prediction": (row.get("prediction") or "Not specified").strip(),
                        "score_name": (row.get("score_name") or "Score").strip(),
                        "score": (row.get("score") or "").strip(),
                        "score_unit": (row.get("score_unit") or "").strip(),
                        "ppv": (row.get("ppv") or "").strip(),
                        "threshold": (row.get("threshold") or "").strip(),
                        "run_date": (row.get("run_date") or "").strip(),
                    })

    abamp_by_seq: dict[str, list[str]] = {}
    with ABAMP.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            seq = (row.get("Sequence") or "").strip()
            if seq:
                abamp_by_seq.setdefault(seq, []).append(row.get("Name") or row.get("AMP") or "AbAMPdb entry")

    catalog = []
    for path in sorted(DETAILS.glob("*.json"), key=lambda p: int(p.stem)):
        detail = json.loads(path.read_text(encoding="utf-8"))
        # Keep source casing: lowercase letters can encode stereochemistry or
        # source-specific residue notation and must not be silently normalized.
        sequence = str(detail.get("sequence") or "").strip()
        if not sequence:
            continue
        observations = []
        for item in detail.get("hemoliticCytotoxicActivities") or []:
            cell = item.get("targetCell") or {}
            group = item.get("activityMeasureForLysisGroup") or {}
            ref = str(item.get("reference") or "")
            paper = article_candidate(detail, ref)
            observations.append({
                "target_cell": str(cell.get("name") or "Not specified"),
                "endpoint_group": str(group.get("name") or "Endpoint not specified"),
                "reported_result": str(item.get("activityMeasureForLysisValue") or "Not specified"),
                "concentration": str(item.get("concentration") or "Not specified"),
                "unit": str((item.get("unit") or {}).get("name") or "Not specified"),
                "pubmed_id": paper["pubmed"],
                "paper_title": paper["title"],
                "reference_status": paper["status"],
            })
        activity_observations = []
        for item in detail.get("targetActivities") or []:
            measure = item.get("activityMeasureGroup") or {}
            target = item.get("targetSpecies") or {}
            if str(measure.get("name") or "").strip().casefold() != "mic":
                continue
            if not str(target.get("name") or "").casefold().startswith("acinetobacter baumannii"):
                continue
            paper = article_candidate(detail, str(item.get("reference") or ""))
            medium = item.get("medium") or {}
            activity_observations.append({
                "target": str(target.get("name") or "Not specified"),
                "reported_result": str(item.get("activityMeasureValue") or "MIC"),
                "concentration": str(item.get("concentration") or "Not specified"),
                "unit": str((item.get("unit") or {}).get("name") or "Not specified"),
                "medium": str(medium.get("name") or "Not specified"),
                "cfu": str(item.get("cfu") or "Not specified"),
                "pubmed_id": paper["pubmed"],
                "paper_title": paper["title"],
                "reference_status": paper["status"],
            })
        aliases = sorted(set(abamp_by_seq.get(sequence, [])))
        catalog.append({
            "id": f"DBAASP:{detail.get('id', path.stem)}",
            "record_id": str(detail.get("id", path.stem)),
            "name": str(detail.get("name") or "Unnamed peptide").strip(),
            "aliases": aliases,
            "sequence": sequence,
            "source": "DBAASP",
            "additional_sources": ["AbAMPdb"] if aliases else [],
            "modification_note": {
                "n_terminus": (detail.get("nTerminus") or detail.get("NTerminus") or {}).get("name", ""),
                "c_terminus": (detail.get("cTerminus") or detail.get("CTerminus") or {}).get("name", ""),
                "unusual_amino_acids": detail.get("unusualAminoAcids") or [],
            },
            "toxicity_observations": observations,
            "activity_observations": activity_observations,
            "toxicity_predictions": [],
            **docking_fields(sequence, planned_sequences, observed_runs, web_jobs),
            "docking_target": "FtsZ from Acinetobacter baumannii",
        })

    # Add AbAMPdb-only rows so its source names and sequences are searchable.
    dbaasp_sequences = {entry["sequence"] for entry in catalog}
    with ABAMP.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            sequence = (row.get("Sequence") or "").strip()
            if not sequence or sequence in dbaasp_sequences:
                continue
            catalog.append({
                "id": f"AbAMPdb:{row.get('source_record_id') or row.get('AMP') or row.get('Name')}",
                "record_id": str(row.get("source_record_id") or row.get("AMP") or ""),
                "name": str(row.get("Name") or row.get("AMP") or "Unnamed AbAMPdb peptide"),
                "aliases": [],
                "sequence": sequence,
                "source": "AbAMPdb",
                "additional_sources": [],
                "modification_note": {"n_terminus": "", "c_terminus": "", "unusual_amino_acids": []},
                "toxicity_observations": [],
                "activity_observations": [],
                "toxicity_predictions": [],
                **docking_fields(sequence, planned_sequences, observed_runs, web_jobs),
                "docking_target": "FtsZ from Acinetobacter baumannii",
            })

    if GENERATED.exists():
        with GENERATED.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                candidate_status = row["candidate_status"]
                if candidate_status == "generated_unvalidated_hypothesis":
                    source = "Source-informed hypothesis"
                elif candidate_status == "reference_parent":
                    source = "Reference control"
                else:
                    source = "Published control"
                catalog.append({
                    "id": f"Generated:{row['candidate_label']}",
                    "record_id": row["candidate_label"],
                    "name": row["candidate_label"],
                    "aliases": [],
                    "sequence": row["sequence"],
                    "source": source,
                    "additional_sources": [],
                    "modification_note": {"n_terminus": "", "c_terminus": "Amidation assumed from B2RP parent; confirm", "unusual_amino_acids": []},
                    "toxicity_observations": [],
                    "activity_observations": [],
                    "toxicity_predictions": [],
                    **docking_fields(row["sequence"], planned_sequences, observed_runs, web_jobs),
                    "docking_target": "FtsZ from Acinetobacter baumannii",
                    "candidate_status": candidate_status,
                })

    by_sequence: dict[str, list[dict[str, str]]] = {}
    for prediction_rows in predictions_by_id.values():
        for prediction in prediction_rows:
            if prediction["sequence"]:
                by_sequence.setdefault(prediction["sequence"], []).append(prediction)
    known_sequences = {entry["sequence"] for entry in catalog}
    unmatched_prediction_rows = sum(
        len(prediction_rows) for candidate_id, prediction_rows in predictions_by_id.items()
        if not any((candidate_id in {entry["id"], entry["record_id"]} and row["sequence"] == entry["sequence"])
                   or row["sequence"] in known_sequences
                   for entry in catalog for row in prediction_rows)
    )
    for entry in catalog:
        imported = predictions_by_id.get(entry["id"], [])
        if not imported:
            imported = predictions_by_id.get(entry["record_id"], [])
        imported = [row for row in imported if row["sequence"] == entry["sequence"]]
        if not imported:
            imported = by_sequence.get(entry["sequence"], [])
        entry["toxicity_predictions"] = []
        for prediction in imported:
            if prediction["sequence"] and prediction["sequence"] == entry["sequence"]:
                entry["toxicity_predictions"].append({k: v for k, v in prediction.items() if k != "sequence"})
            else:
                unmatched_prediction_rows += 1
    catalog.sort(key=lambda entry: (entry["name"].casefold(), entry["id"]))
    result = {
        "snapshot_date": date.today().isoformat(),
        "record_count": len(catalog),
        "docking_note": "No peptide docking run is available. A docking score is never imputed.",
        "toxicity_note": "Observations are source-reported assays, not a universal safe/unsafe label. No observation means unknown.",
        "records": catalog,
    }
    UI_DATA.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(result, ensure_ascii=False, separators=(",", ":"))
    UI_DATA.write_text(serialized, encoding="utf-8")
    UI_DATA_JS.write_text(f"window.PEPTIDE_DATA = {serialized};\n", encoding="utf-8")
    print(f"Built UI data for {len(catalog)} records: {UI_DATA.relative_to(ROOT)}; unmatched/mismatched prediction rows: {unmatched_prediction_rows}")


if __name__ == "__main__":
    main()
