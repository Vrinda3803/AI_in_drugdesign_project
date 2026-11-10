"""Compute transparent, sequence-only features for unique standard-AA peptides."""
from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "data" / "processed" / "sequence_provenance_panel.csv"
OUTPUT = ROOT / "data" / "processed" / "sequence_features.csv"
AA = "ACDEFGHIKLMNPQRSTVWY"

# IUPAC average free-amino-acid masses (Da); subtract H2O per residue then add
# one H2O for an unmodified, free-N/free-C terminal linear peptide.
FREE_AA_MASS = dict(zip(AA, [89.0932,121.1590,133.1032,147.1299,165.1900,75.0666,155.1552,131.1736,146.1882,131.1736,149.2124,132.1184,115.1310,146.1451,174.2017,105.0926,119.1192,204.2262,181.1894,117.1469]))
WATER = 18.01528

# Kyte-Doolittle scale (GRAVY).
KD = dict(zip(AA, [1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,-0.9,-1.3,4.2]))

# Eisenberg consensus hydrophobicity scale, normalized values.
EISENBERG = dict(zip(AA, [.62,.29,-.90,-.74,1.19,.48,-.40,1.38,-.80,1.06,.64,-.85,.12,-.78,-.74,-.18,-.05,.81,.26,1.08]))

# Bjellqvist pK values used by Biopython's IsoelectricPoint implementation.
POS_PKS = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
NEG_PKS = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
NTERM_PKS = {"A":7.59,"M":7.0,"S":6.93,"P":8.36,"T":6.82,"V":7.44,"E":7.7}
CTERM_PKS = {"D":4.55,"E":4.75}


def charge(sequence: str, ph: float) -> float:
    positive = POS_PKS.copy()
    negative = NEG_PKS.copy()
    if sequence[0] in NTERM_PKS:
        positive["Nterm"] = NTERM_PKS[sequence[0]]
    if sequence[-1] in CTERM_PKS:
        negative["Cterm"] = CTERM_PKS[sequence[-1]]
    count = Counter(sequence)
    pos = sum(count.get(aa, 0) / (10 ** (ph - pk) + 1) for aa, pk in positive.items() if aa != "Nterm")
    pos += 1 / (10 ** (ph - positive["Nterm"]) + 1)
    neg = sum(count.get(aa, 0) / (10 ** (pk - ph) + 1) for aa, pk in negative.items() if aa != "Cterm")
    neg += 1 / (10 ** (negative["Cterm"] - ph) + 1)
    return pos - neg


def isoelectric_point(sequence: str) -> float:
    low, high = 4.05, 12.0
    for _ in range(40):
        mid = (low + high) / 2
        if charge(sequence, mid) > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def hydrophobic_moment(sequence: str, scale: dict[str, float], angle_deg: float = 100.0) -> float:
    angle = math.radians(angle_deg)
    x = sum(scale[aa] * math.cos(i * angle) for i, aa in enumerate(sequence))
    y = sum(scale[aa] * math.sin(i * angle) for i, aa in enumerate(sequence))
    return math.hypot(x, y) / len(sequence)


def maximum_11mer_moment(sequence: str) -> float:
    width = min(11, len(sequence))
    return max(hydrophobic_moment(sequence[i:i + width], EISENBERG) for i in range(len(sequence) - width + 1))


def main() -> None:
    with PANEL.open(encoding="utf-8-sig", newline="") as handle:
        panel = list(csv.DictReader(handle))
    rows: list[dict] = []
    for item in panel:
        seq = item["sequence"]
        if not seq or any(aa not in AA for aa in seq):
            continue
        count = Counter(seq)
        row = {"sequence": seq, "length": len(seq), "sources_and_ids_json": item["sources_and_ids_json"],
               "source_record_count": item["source_record_count"],
               "molecular_weight_average_unmodified_g_mol": sum(FREE_AA_MASS[aa] - WATER for aa in seq) + WATER,
               "net_charge_unmodified_pH7": charge(seq, 7.0), "isoelectric_point_unmodified": isoelectric_point(seq),
               "gravy_kyte_doolittle": sum(KD[aa] for aa in seq) / len(seq),
               "hydrophobic_moment_eisenberg_100deg_full_sequence": hydrophobic_moment(seq, EISENBERG),
               "hydrophobic_moment_eisenberg_100deg_max_11mer": maximum_11mer_moment(seq),
               "hydrophobic_moment_window_length": min(11, len(seq)),
               "descriptor_assumption": "linear peptide, free unmodified N/C termini; theoretical sequence descriptors; no PTM/cyclization correction"}
        row.update({f"AAC_{aa}": count.get(aa, 0) / len(seq) for aa in AA})
        rows.append(row)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "feature_rows": len(rows),
        "feature_columns": fields,
        "molecular_weight": "Average free amino acid masses minus water per peptide bond plus one water; linear peptide with unmodified free termini.",
        "net_charge_and_pI": "Henderson-Hasselbalch charge calculation and bisection pI using Bjellqvist pKa values, matching the documented Biopython IsoelectricPoint method; free termini assumed.",
        "hydrophobicity": "Mean Kyte-Doolittle scale (GRAVY).",
        "hydrophobic_moment": "Eisenberg consensus scale; alpha-helix angle 100 degrees; full-sequence mean vector magnitude and maximum sliding window up to 11 residues.",
        "modifications": "Sequence-derived theoretical descriptors do not account for terminal chemistry, cyclization, D-amino acids, or other noncanonical modifications.",
        "source_refs": ["https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html", "https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/IsoelectricPoint.py", "https://doi.org/10.1016/0022-2836(84)90309-7"]
    }
    (OUTPUT.parent / "feature_engineering_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({"feature_rows": len(rows), "feature_columns": len(fields)}, ensure_ascii=True))


if __name__ == "__main__":
    main()
