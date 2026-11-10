"""Prepare single-chain protein-only receptor files from downloaded structures."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

HERE = Path(__file__).resolve().parent
STRUCTURES = HERE / "structures"
OUT = STRUCTURES / "prepared"
INPUTS = {
    "A_baumannii_apo_22HR": (STRUCTURES / "A_baumannii_apo_22HR.cif", "cif"),
    "E_coli_FtsZ_6LL6": (STRUCTURES / "E_coli_FtsZ_6LL6.cif", "cif"),
    "S_aureus_apo_3VPA": (STRUCTURES / "S_aureus_apo_3VPA.cif", "cif"),
    "E_coli_FtsZ_AFDB": (STRUCTURES / "E_coli_FtsZ_AFDB.pdb", "pdb"),
    "S_aureus_FtsZ_AFDB": (STRUCTURES / "S_aureus_FtsZ_AFDB.pdb", "pdb"),
}


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def seq_for_chain(chain):
    residues = [r for r in chain if is_aa(r, standard=True) and r.resname in protein_letters_3to1]
    return "".join(protein_letters_3to1[r.resname] for r in residues), residues


class ProteinChainSelect(Select):
    def __init__(self, chain_id): self.chain_id = chain_id
    def accept_chain(self, chain): return int(chain.id == self.chain_id)
    def accept_residue(self, residue): return int(is_aa(residue, standard=True) and residue.resname in protein_letters_3to1)
    def accept_atom(self, atom): return int(atom.element not in {"H", "D"})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {"preparation": "Select longest standard-amino-acid chain; save protein atoms only, excluding heterogens, water, hydrogen and deuterium.", "targets": {}}
    for label, (source, kind) in INPUTS.items():
        if not source.exists(): raise SystemExit(f"Missing structure source: {source}; run fetch_target_structures.py first.")
        parser = MMCIFParser(QUIET=True) if kind == "cif" else PDBParser(QUIET=True)
        structure = parser.get_structure(label, str(source))
        model = next(structure.get_models())
        choices = [(len(seq_for_chain(chain)[0]), chain.id, seq_for_chain(chain)[0]) for chain in model]
        length, chain_id, sequence = max(choices, default=(0, "", ""))
        if length < 100: raise ValueError(f"No plausible FtsZ protein chain found in {source.name} (max length {length}).")
        output = OUT / f"{label}_protein_only.pdb"
        io = PDBIO(); io.set_structure(structure); io.save(str(output), ProteinChainSelect(chain_id))
        manifest["targets"][label] = {"source_file": str(source.relative_to(HERE)), "source_sha256": sha256(source),
                                      "source_chain_id": chain_id, "selected_sequence_length": length,
                                      "selected_sequence": sequence, "prepared_file": str(output.relative_to(HERE)),
                                      "prepared_sha256": sha256(output),
                                      "caveat": "Protein-only inspection/receptor input; receptor protonation, side-chain optimization and docking validation are not performed here."}
    (OUT / "prepared_target_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({name: {k: v for k, v in info.items() if k not in {"selected_sequence"}} for name, info in manifest["targets"].items()}, indent=2))


if __name__ == "__main__": main()
