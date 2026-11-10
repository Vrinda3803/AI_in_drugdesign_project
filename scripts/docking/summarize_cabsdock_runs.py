#!/usr/bin/env python3
"""Inventory CABSdock outputs against the prepared run matrix.

This collector records artifacts and completion evidence. It does not extract
or invent binding-affinity values from pose ranks, cluster labels, plots, or
PDB metadata.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN_ROOT = ROOT / "results/docking/cabsdock"
MATRIX = RUN_ROOT / "run_matrix.csv"
OUTPUTS = RUN_ROOT / "runs"
RUN_CSV = RUN_ROOT / "observed_runs.csv"
SUMMARY_JSON = RUN_ROOT / "observed_summary.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inventory(row: dict[str, str]) -> dict[str, str | int]:
    work = OUTPUTS / row["run_id"]
    models_dir = work / "output_pdbs"
    model_paths = sorted(
        (path for path in models_dir.glob("model_*.pdb")
         if path.stem.removeprefix("model_").isdigit()),
        key=lambda path: int(path.stem.removeprefix("model_")),
    ) if models_dir.is_dir() else []
    cluster_path = work / "output_data/cluster_names.txt"
    log_path = work / "CABS.log"
    log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
    success_marker = "simulation completed successfully" in log_text.casefold()
    cluster_present = cluster_path.is_file() and cluster_path.stat().st_size > 0
    models_present = len(model_paths) > 0
    status = "complete" if success_marker and cluster_present and models_present else (
        "partial" if (log_path.exists() or cluster_path.exists() or models_present) else "not_run"
    )
    return {
        **row,
        "status": status,
        "success_marker_in_log": int(success_marker),
        "ranked_pose_model_count": len(model_paths),
        "model_1_present": int((models_dir / "model_1.pdb").is_file()),
        "cluster_density_file_present": int(cluster_present),
        "cluster_density_file_sha256": sha256(cluster_path) if cluster_present else "",
        "run_log_present": int(log_path.is_file()),
        "run_log_sha256": sha256(log_path) if log_path.is_file() else "",
        "work_directory": str(work.relative_to(ROOT)) if work.exists() else "",
    }


def main() -> None:
    if not MATRIX.is_file():
        raise SystemExit(f"Run matrix not found: {MATRIX}; prepare it first.")
    with MATRIX.open("r", encoding="utf-8-sig", newline="") as stream:
        planned = list(csv.DictReader(stream))
    observed = [inventory(row) for row in planned]
    fields = list(observed[0]) if observed else []
    with RUN_CSV.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(observed)

    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in observed:
        grouped[(str(row["candidate_id"]), str(row["receptor_file"]))].append(row)
    candidate_summary = []
    for (candidate_id, receptor_file), rows in sorted(grouped.items()):
        complete = sum(row["status"] == "complete" for row in rows)
        candidate_summary.append({
            "candidate_id": candidate_id,
            "sequence": rows[0]["sequence"],
            "receptor_file": receptor_file,
            "planned_seed_count": len(rows),
            "complete_run_count": complete,
            "partial_run_count": sum(row["status"] == "partial" for row in rows),
            "not_run_count": sum(row["status"] == "not_run" for row in rows),
            "total_ranked_pose_models": sum(int(row["ranked_pose_model_count"]) for row in rows),
        })
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "planned_run_count": len(planned),
        "complete_run_count": sum(row["status"] == "complete" for row in observed),
        "partial_run_count": sum(row["status"] == "partial" for row in observed),
        "not_run_count": sum(row["status"] == "not_run" for row in observed),
        "score_interpretation": (
            "No numeric affinity or docking score is calculated here. Pose model indices are ranks; "
            "cluster-density artifacts are retained by path and checksum for review."
        ),
        "completion_evidence": (
            "CABS.log contains 'Simulation completed successfully', at least one output_pdbs/model_N.pdb "
            "exists, and output_data/cluster_names.txt is nonempty."
        ),
        "official_output_reference": "https://cabsflex.lcbio.pl/peptide-protein-docking.html",
        "run_inventory_csv": str(RUN_CSV.relative_to(ROOT)),
        "candidate_receptor_summary": candidate_summary,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in (
        "planned_run_count", "complete_run_count", "partial_run_count", "not_run_count",
        "run_inventory_csv"
    )}, indent=2))


if __name__ == "__main__":
    main()
