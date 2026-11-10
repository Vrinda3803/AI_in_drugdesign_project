"""Audit generated LSTM sequences against local AMP data; not an activity predictor."""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
FASTA = OUT / "lstm_hypotheses.fasta"
PANEL = ROOT / "data/processed/sequence_provenance_panel.csv"
AA = "ACDEFGHIKLMNPQRSTVWY"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                rows.append((name, "".join(seq).upper()))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name is not None:
        rows.append((name, "".join(seq).upper()))
    return rows


def edit_distance(a: str, b: str) -> int:
    """Levenshtein distance with O(min(len(a), len(b))) memory."""
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, aa in enumerate(a, 1):
        curr = [i]
        for j, bb in enumerate(b, 1):
            curr.append(min(curr[-1] + 1, prev[j] + 1, prev[j - 1] + (aa != bb)))
        prev = curr
    return prev[-1]


def main() -> None:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise SystemExit("Install Biopython in .venv to calculate instability index.") from exc
    import sys
    sys.path.insert(0, str(ROOT))
    from scripts.features.engineer_sequence_features import (
        EISENBERG, FREE_AA_MASS, KD, WATER, charge, hydrophobic_moment,
        isoelectric_point,
    )

    training: set[str] = set()
    with PANEL.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            seq = (row.get("sequence") or "").strip().upper()
            if 8 <= len(seq) <= 80 and set(seq) <= set(AA):
                training.add(seq)
    if not training:
        raise SystemExit(f"No valid training sequences found in {PANEL}")

    rows: list[dict] = []
    sequences = read_fasta(FASTA)
    invalid = []
    for identifier, seq in sequences:
        if not seq or not set(seq) <= set(AA):
            invalid.append(identifier)
            continue
        closest, distance = min(
            ((known, edit_distance(seq, known)) for known in training),
            key=lambda pair: pair[1] / max(len(seq), len(pair[0])),
        )
        counts = Counter(seq)
        rows.append({
            "candidate_id": identifier.split("|")[0],
            "sequence": seq,
            "length_aa": len(seq),
            "molecular_weight_free_termini_g_mol": round(
                sum(FREE_AA_MASS[aa] - WATER for aa in seq) + WATER, 3
            ),
            "net_charge_pH7_free_termini": round(charge(seq, 7.0), 3),
            "theoretical_pI_free_termini": round(isoelectric_point(seq), 3),
            "gravy_kyte_doolittle": round(sum(KD[aa] for aa in seq) / len(seq), 3),
            "hydrophobic_moment_eisenberg_100deg": round(hydrophobic_moment(seq, EISENBERG), 3),
            "instability_index_biopython": round(ProteinAnalysis(seq).instability_index(), 3),
            "basic_residue_fraction_KR": round((counts["K"] + counts["R"]) / len(seq), 3),
            "largest_single_residue_fraction": round(max(counts.values()) / len(seq), 3),
            "closest_local_training_sequence": closest,
            "nearest_levenshtein_identity_estimate": round(
                1 - distance / max(len(seq), len(closest)), 4
            ),
            "exact_training_match": seq in training,
            "interpretation": "Sequence-derived descriptors and local similarity audit only; no activity or safety inference.",
        })

    fields = list(rows[0]) if rows else []
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "lstm_candidate_audit.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "input_fasta": str(FASTA),
        "training_panel": str(PANEL),
        "candidates_audited": len(rows),
        "training_sequences_compared": len(training),
        "invalid_sequences": invalid,
        "exact_training_matches": sum(row["exact_training_match"] for row in rows),
        "nearest_identity_estimate_definition": "1 - Levenshtein edit distance / max(candidate length, training sequence length); approximate global sequence similarity, not a biological alignment score.",
        "feature_assumptions": "Linear peptide with unmodified free N/C termini; theoretical descriptors do not account for PTMs, cyclization, D-amino acids, or noncanonical residues.",
        "limitations": [
            "Descriptors do not predict antibacterial activity, toxicity, hemolysis, target binding, stability in vivo, or manufacturability.",
            "Local similarity is not a comprehensive novelty search against external patent or sequence databases.",
            "The LSTM validation loss is higher than its training loss, consistent with overfitting.",
        ],
    }
    (OUT / "lstm_candidate_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
