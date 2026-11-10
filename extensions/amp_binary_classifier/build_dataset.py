"""Build AMP/putative-non-AMP descriptors without modifying the main pipeline."""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "output"
NEGATIVE_CSV = Path(__file__).resolve().parent / "data" / "non_amp_sequences.csv"
AA = "ACDEFGHIKLMNPQRSTVWY"

def instability_index(seq: str) -> float:
    """Return ProtParam-style instability index (requires complete DIWV table)."""
    # Import the project's installed Biopython only if available; otherwise the
    # extension reports that this feature cannot be computed, not a fake score.
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise RuntimeError("Instability index requires Biopython. Install it in .venv before building.") from exc
    return float(ProteinAnalysis(seq).instability_index())


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def positive_sequences() -> dict[str, str]:
    seqs: dict[str, str] = {}
    for path in (ROOT / "data/processed/sequence_provenance_panel.csv",):
        for row in read_csv(path):
            seq = (row.get("sequence") or "").strip().upper()
            if 5 <= len(seq) <= 100 and set(seq) <= set(AA):
                seqs[seq] = "local_AbAMPdb_DBAASP_sequence_panel"
    return seqs


def descriptors(seq: str) -> dict[str, float | int]:
    counts = Counter(seq)
    # These simple descriptors complement the richer main-project table and
    # make the extension usable independently from its exact feature schema.
    sys.path.insert(0, str(ROOT))
    from scripts.features.engineer_sequence_features import charge, isoelectric_point, hydrophobic_moment, KD, EISENBERG
    return {
        "length": len(seq),
        "molecular_weight_g_mol": sum({"A":89.0932,"C":121.1590,"D":133.1032,"E":147.1299,"F":165.1900,"G":75.0666,"H":155.1552,"I":131.1736,"K":146.1882,"L":131.1736,"M":149.2124,"N":132.1184,"P":115.1310,"Q":146.1451,"R":174.2017,"S":105.0926,"T":119.1192,"V":117.1469,"W":204.2262,"Y":181.1894}[a] for a in seq) - (len(seq)-1)*18.01528,
        "net_charge_pH7": charge(seq, 7.0),
        "isoelectric_point": isoelectric_point(seq),
        "hydrophobicity_gravy": sum(KD[a] for a in seq) / len(seq),
        "hydrophobic_moment": hydrophobic_moment(seq, EISENBERG),
        "instability_index": instability_index(seq),
        **{f"AAC_{a}": counts.get(a, 0) / len(seq) for a in AA},
    }


def main() -> None:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis  # noqa: F401
    except ImportError as exc:
        raise SystemExit("Biopython is required for instability index. Install extensions/amp_binary_classifier/requirements.txt in .venv.") from exc
    if not NEGATIVE_CSV.exists():
        raise SystemExit(f"Missing curated negative set: {NEGATIVE_CSV}\nSee README.md for schema and labeling cautions.")
    pos = positive_sequences()
    negatives = read_csv(NEGATIVE_CSV)
    combined: dict[str, dict] = {}
    rejected, conflicts = [], []
    for seq, source in pos.items():
        combined[seq] = {"sequence": seq, "label": 1, "source": source, "source_id": "local_panel"}
    for i, row in enumerate(negatives, 2):
        seq = (row.get("sequence") or "").strip().upper()
        if row.get("label", "0").strip() not in {"0", "non_amp", "negative"}:
            rejected.append({"row": i, "reason": "negative label must be 0/non_amp/negative"}); continue
        if not (5 <= len(seq) <= 100) or not set(seq) <= set(AA):
            rejected.append({"row": i, "reason": "sequence must be 5-100 standard amino acids"}); continue
        if seq in combined:
            conflicts.append({"sequence": seq, "negative_source": row.get("source", ""), "action": "excluded from negative class"}); continue
        combined[seq] = {"sequence": seq, "label": 0, "source": row.get("source", ""), "source_id": row.get("source_id", "")}
    rows = []
    for item in combined.values():
        try:
            rows.append({**item, **descriptors(item["sequence"])})
        except Exception as exc:
            rejected.append({"sequence": item["sequence"], "reason": f"feature_error: {exc}"})
    fields = list(rows[0]) if rows else []
    write_csv(OUT / "amp_binary_dataset.csv", rows, fields)
    report = {
        "label_definition": "1=sequence present in local AbAMPdb/DBAASP panel; 0=UniProt-derived, length-matched putative-negative fragment",
        "class_counts": {str(label): sum(r["label"] == label for r in rows) for label in (0, 1)},
        "negative_input_path": str(NEGATIVE_CSV), "rejected_rows": rejected, "positive_negative_exact_conflicts": conflicts,
        "limitations": ["Negative fragments are weak labels and may have antimicrobial activity.", "Lengths are sampled to match positives, but biological source and sequence-composition bias may inflate scores.", "Grouping is only a leakage control, not independent external validation."],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "dataset_qc.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
