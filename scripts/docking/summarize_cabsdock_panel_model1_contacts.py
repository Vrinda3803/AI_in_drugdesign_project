#!/usr/bin/env python3
"""Compare geometric contacts in the representative model for each CABS job."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from analyze_cabsdock_pose_contacts import ROOT, chain_sequence, parse_pdb, summarize_pose

JOBS = ROOT / "results/docking/cabsdock/web_jobs.csv"
MODELS = ROOT / "results/docking/cabsdock/web_models"
REFERENCE = ROOT / "data/external/ftsZ/prepared/22HR_chainA_apo.pdb"
OUTPUT_CSV = ROOT / "results/docking/cabsdock/web_panel_model1_contact_summary.csv"
OUTPUT_JSON = ROOT / "results/docking/cabsdock/web_panel_model1_contact_summary.json"


def main() -> None:
    with JOBS.open("r", encoding="utf-8-sig", newline="") as handle:
        jobs = list(csv.DictReader(handle))
    target_chains = parse_pdb(REFERENCE)
    target_sequence = chain_sequence(target_chains["A"])
    rows: list[dict] = []
    for job in jobs:
        job_id = job["job_url"].rstrip("/").split("/")[-1]
        job_dir = MODELS / job_id
        candidates = [job_dir / "model_1.pdb", *sorted(job_dir.glob("*_model_1.pdb"))]
        model = next((path for path in candidates if path.exists()), None)
        if model is None:
            raise FileNotFoundError(f"Representative model_1 not found for {job['candidate_id']} at {job_dir}")
        row = summarize_pose(model, target_sequence, job["sequence"])
        row.update({
            "candidate_id": job["candidate_id"],
            "job_id": job_id,
            "job_url": job["job_url"],
            "CABS_top_cluster_density": job.get("top_cluster_density", ""),
            "CABS_top_cluster_size": job.get("top_cluster_size", ""),
        })
        rows.append(row)
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with OUTPUT_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "method": "Heavy-atom geometric contact summary for model_1 (most representative model of each CABS-dock best cluster).",
        "models": len(rows),
        "receptor": "A. baumannii FtsZ 22HR chain A, apo monomer",
        "cutoffs_angstrom": [4, 5, 8],
        "interpretation_limit": "Geometric counts and minimum distances do not represent docking energy, binding affinity, or antibacterial activity. Candidate cross-comparisons are exploratory single-run summaries.",
        "csv": str(OUTPUT_CSV.relative_to(ROOT)),
        "results": rows,
    }
    OUTPUT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"candidate_models": len(rows), "csv": str(OUTPUT_CSV), "json": str(OUTPUT_JSON)}, indent=2))


if __name__ == "__main__":
    main()
