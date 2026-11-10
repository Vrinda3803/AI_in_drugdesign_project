#!/usr/bin/env python3
"""Summarize chain-A/chain-D heavy-atom contacts in downloaded CABS-dock PDBs.

This is a geometric contact analysis only. It does not calculate an affinity,
free energy, or biological activity score.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "results/docking/cabsdock/web_models/3682c266f3d6adc"
DEFAULT_REFERENCE = ROOT / "data/external/ftsZ/prepared/22HR_chainA_apo.pdb"
ONE_TO_THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}
THREE_TO_ONE = {three: one for one, three in ONE_TO_THREE.items()}
MODEL_NAME = re.compile(r"(?:_model_|^model_)(\d+)\.pdb$", re.IGNORECASE)
EXPECTED_SEQUENCE = "GIWKTIKSMGKVFAGAIKQNL"


def parse_pdb(path: Path) -> dict[str, dict]:
    chains: dict[str, dict] = {}
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        chain = line[21].strip() or "_"
        residue_key = (line[22:26].strip(), line[26:27].strip())
        residue_name = line[17:20].strip().upper()
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        atom_name = line[12:16].strip().upper()
        if not element:
            element = re.sub(r"^[0-9]+", "", atom_name)[:1]
        chain_data = chains.setdefault(chain, {"residues": OrderedDict(), "atoms": []})
        chain_data["residues"].setdefault(residue_key, residue_name)
        if element not in {"H", "D"}:
            chain_data["atoms"].append((xyz, residue_key, residue_name, atom_name))
    return chains


def chain_sequence(chain: dict) -> str:
    return "".join(THREE_TO_ONE.get(name, "X") for name in chain["residues"].values())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize_pose(path: Path, target_sequence: str, expected_sequence: str = EXPECTED_SEQUENCE) -> dict:
    match = MODEL_NAME.search(path.name)
    if not match:
        raise ValueError(f"Filename must end with _model_N.pdb: {path.name}")
    chains = parse_pdb(path)
    if "A" not in chains or "D" not in chains:
        raise ValueError(f"{path.name} must contain receptor chain A and peptide chain D")
    receptor_sequence = chain_sequence(chains["A"])
    peptide_sequence = chain_sequence(chains["D"])
    if receptor_sequence != target_sequence:
        raise ValueError(f"Target chain A in {path.name} does not match the prepared 22HR chain A")
    if peptide_sequence != expected_sequence:
        raise ValueError(f"Peptide chain D in {path.name} does not match the expected panel sequence")

    atoms_a = chains["A"]["atoms"]
    atoms_d = chains["D"]["atoms"]
    min_distance = math.inf
    atom_pair_counts = {4.0: 0, 5.0: 0, 8.0: 0}
    residue_pairs = {cutoff: set() for cutoff in atom_pair_counts}
    target_residues = {cutoff: set() for cutoff in atom_pair_counts}
    peptide_residues = {cutoff: set() for cutoff in atom_pair_counts}
    for (xa, ya, za), key_a, name_a, _ in atoms_a:
        for (xd, yd, zd), key_d, name_d, _ in atoms_d:
            dx, dy, dz = xa - xd, ya - yd, za - zd
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < min_distance * min_distance:
                min_distance = math.sqrt(d2)
            for cutoff in atom_pair_counts:
                if d2 <= cutoff * cutoff:
                    atom_pair_counts[cutoff] += 1
                    residue_pairs[cutoff].add((key_a, name_a, key_d, name_d))
                    target_residues[cutoff].add(key_a)
                    peptide_residues[cutoff].add(key_d)

    def residue_list(keys: set[tuple[str, str]]) -> str:
        return ";".join(f"{number}{insertion}" for number, insertion in sorted(keys, key=lambda x: (int(x[0]), x[1])))

    return {
        "model_rank": int(match.group(1)),
        "file": path.name,
        "sha256": sha256(path),
        "peptide_sequence": peptide_sequence,
        "target_sequence_length": len(receptor_sequence),
        "receptor_heavy_atom_count": len(atoms_a),
        "peptide_heavy_atom_count": len(atoms_d),
        "minimum_heavy_atom_distance": f"{min_distance:.3f}",
        "atom_contacts_le_4A": atom_pair_counts[4.0],
        "residue_contacts_le_4A": len(residue_pairs[4.0]),
        "target_residue_contacts_le_4A": len(target_residues[4.0]),
        "peptide_residue_contacts_le_4A": len(peptide_residues[4.0]),
        "atom_contacts_le_5A": atom_pair_counts[5.0],
        "residue_contacts_le_5A": len(residue_pairs[5.0]),
        "target_residue_contacts_le_5A": len(target_residues[5.0]),
        "peptide_residue_contacts_le_5A": len(peptide_residues[5.0]),
        "target_interface_residues_le_5A": residue_list(target_residues[5.0]),
        "peptide_interface_residues_le_5A": residue_list(peptide_residues[5.0]),
        "atom_contacts_le_8A": atom_pair_counts[8.0],
        "target_residue_contacts_le_8A": len(target_residues[8.0]),
        "peptide_residue_contacts_le_8A": len(peptide_residues[8.0]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--reference-pdb", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/docking/cabsdock")
    parser.add_argument("--job-id", default="3682c266f3d6adc")
    parser.add_argument("--candidate-id", default="B2RP_top_published_lead_control")
    parser.add_argument("--expected-sequence", default=EXPECTED_SEQUENCE)
    parser.add_argument("--expected-model-count", type=int, default=10)
    args = parser.parse_args()
    target_chains = parse_pdb(args.reference_pdb)
    if "A" not in target_chains:
        raise ValueError(f"Reference PDB has no chain A: {args.reference_pdb}")
    target_sequence = chain_sequence(target_chains["A"])
    pose_files = sorted(
        (path for path in args.input_dir.glob("3682c266f3d6adc_model_*.pdb") if MODEL_NAME.search(path.name)),
        key=lambda path: int(MODEL_NAME.search(path.name).group(1)),
    )
    if not pose_files:
        raise ValueError(f"No model PDB files found in {args.input_dir}")
    if args.expected_model_count and len(pose_files) != args.expected_model_count:
        raise ValueError(f"Expected {args.expected_model_count} model PDBs, found {len(pose_files)} in {args.input_dir}")
    rows = [summarize_pose(path, target_sequence, args.expected_sequence) for path in pose_files]
    for row in rows:
        row["candidate_id"] = args.candidate_id
        row["job_id"] = args.job_id
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / f"web_model_contact_summary_{args.job_id}.csv"
    json_path = args.output_dir / f"web_model_contact_summary_{args.job_id}.json"
    fields = list(rows[0])
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "job_id": args.job_id,
        "job_url": f"https://lcbio.pl/cabsdock/job/{args.job_id}/",
        "candidate_id": args.candidate_id,
        "peptide_sequence": args.expected_sequence,
        "receptor": "A. baumannii FtsZ PDB 22HR chain A, apo",
        "source_file_count": len(rows),
        "reference_target_sequence_length": len(target_sequence),
        "analysis": "Heavy-atom geometric contact counts and interface residue counts at 4, 5, and 8 Angstrom cutoffs; minimum inter-chain heavy-atom distance.",
        "interpretation_limit": "Contact counts and distances summarize a modeled pose only. They are not docking energies, binding affinities, target engagement, or antibacterial activity. CABS-dock model rank is occurrence in the docking trajectory, not affinity rank.",
        "csv": str(csv_path.relative_to(ROOT)),
        "models_directory": str(args.input_dir.relative_to(ROOT)),
        "models": rows,
    }
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"models": len(rows), "csv": str(csv_path), "json": str(json_path),
                      "contact_pairs_5A": [row["atom_contacts_le_5A"] for row in rows]}, indent=2))


if __name__ == "__main__":
    main()
