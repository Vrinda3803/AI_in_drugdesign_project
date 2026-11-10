"""Create a reproducible, multi-seed CABSdock command matrix for the panel.

This prepares commands only. It does not install or run CABSdock and produces
no poses or scores.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "results/candidate_prioritization/b2rp_candidate_panel.fasta"
OUT_DIR = ROOT / "results/docking/cabsdock"
SCRIPT = ROOT / "scripts/docking/run_cabsdock_panel.sh"
MATRIX = OUT_DIR / "run_matrix.csv"
MANIFEST = OUT_DIR / "run_plan.json"
SEEDS = [104729, 130363, 155921]
RECEPTORS = ["22HR_chainA_apo.pdb", "22HQ_chainA_GDPbound.pdb"]
SAMPLING = {"mc_annealing": 20, "mc_cycles": 50, "mc_steps": 200}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier = None
    sequence: list[str] = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence)))
                identifier = line[1:].split("|", 1)[0]
                sequence = []
            elif identifier is not None:
                sequence.append(line.upper())
    if identifier is not None:
        records.append((identifier, "".join(sequence)))
    return records


def main() -> None:
    records = read_fasta(PANEL)
    panel_hash = hashlib.sha256(PANEL.read_bytes()).hexdigest()
    receptor_hashes = {
        name: hashlib.sha256((ROOT / "data/external/ftsZ/prepared" / name).read_bytes()).hexdigest()
        for name in RECEPTORS
    }
    rows = []
    for candidate_id, sequence in records:
        for receptor in RECEPTORS:
            for seed in SEEDS:
                run_id = f"{candidate_id}__{Path(receptor).stem}__seed{seed}"
                rows.append({"run_id": run_id, "candidate_id": candidate_id, "sequence": sequence,
                             "receptor_file": receptor, "receptor_sha256": receptor_hashes[receptor],
                             "candidate_panel_sha256": panel_hash, "seed": seed, **SAMPLING})

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with MATRIX.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    shell = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"',
        'RECEPTOR_DIR="$PROJECT_ROOT/data/external/ftsZ/prepared"',
        'RUN_ROOT="$PROJECT_ROOT/results/docking/cabsdock/runs"',
        'command -v CABSdock >/dev/null 2>&1 || { echo "CABSdock is not installed or not on PATH" >&2; exit 127; }',
        'mkdir -p "$RUN_ROOT"',
        'echo "CABSdock path: $(command -v CABSdock)"',
        'echo "Run plan: $PROJECT_ROOT/results/docking/cabsdock/run_matrix.csv"',
        "",
    ]
    for row in rows:
        receptor = f'"$RECEPTOR_DIR/{row["receptor_file"]}"'
        run_dir = f'"$RUN_ROOT/{row["run_id"]}"'
        shell.extend([
            f'if [[ -e {run_dir} ]]; then echo "Output already exists: {row["run_id"]}" >&2; exit 2; fi',
            f'echo "Starting {row["run_id"]}"',
            f'CABSdock -i {receptor} -p "{row["sequence"]}" -a {row["mc_annealing"]} -y {row["mc_cycles"]} -s {row["mc_steps"]} -z {row["seed"]} -M -C -S -o A --csv-output A --json-output --log --no-progress-bar -w {run_dir}',
            "",
        ])
    shell.append('echo "All CABSdock runs finished. Inspect outputs and logs before importing any score."')
    SCRIPT.write_text("\n".join(shell) + "\n", encoding="utf-8", newline="\n")
    plan = {
        "status": "commands_prepared_not_executed",
        "candidate_count": len(records), "receptor_states": RECEPTORS,
        "seeds": SEEDS, "sampling": SAMPLING,
        "candidate_panel_sha256": panel_hash, "receptor_sha256": receptor_hashes,
        "run_count": len(rows), "run_matrix": str(MATRIX.relative_to(ROOT)),
        "shell_script": str(SCRIPT.relative_to(ROOT)),
        "parameters_source": "Official CABS-flex standalone 3 examples and options reference; workflow remains exploratory and unvalidated on an AMP-FtsZ reference complex.",
        "interpretation": "CABSdock scores are method-specific pose rankings, not binding free energy or measured affinity.",
    }
    MANIFEST.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
