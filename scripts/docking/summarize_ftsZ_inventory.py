"""Summarize RCSB-provided metadata and validation fields for the target PDBs."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "data" / "external" / "ftsZ"


def main() -> None:
    result = {}
    for pdb_id in ("22HQ", "22HR"):
        payload = json.loads((BASE / f"{pdb_id}_rcsb_metadata.json").read_text(encoding="utf-8"))
        entry = payload["entry"]
        polymer = payload["polymer_entity_1"]
        assembly = payload["assembly_1"]["pdbx_struct_assembly"]
        info = entry["rcsb_entry_info"]
        entity_poly = polymer["entity_poly"]
        polymer_entity = polymer["rcsb_polymer_entity"]
        diffraction = entry["pdbx_vrpt_summary_diffraction"][0]
        geometry = entry["pdbx_vrpt_summary_geometry"][0]
        organism = (polymer.get("rcsb_entity_source_organism") or [{}])[0]
        accession = polymer["rcsb_polymer_entity_container_identifiers"]
        result[pdb_id] = {
            "title": entry["struct"]["title"],
            "release_date": entry["rcsb_accession_info"]["initial_release_date"],
            "experimental_method": info["experimental_method"],
            "resolution_angstrom": info["resolution_combined"][0],
            "r_work": diffraction["DCC_R"], "r_free": diffraction["DCC_Rfree"],
            "data_completeness_percent": diffraction["data_completeness"],
            "sequence_length": entity_poly["rcsb_sample_sequence_length"],
            "sequence": entity_poly["pdbx_seq_one_letter_code_can"].replace("\n", "").strip(),
            "chain_ids": polymer["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"],
            "entity_description": polymer_entity["pdbx_description"],
            "uniprot_ids": accession.get("uniprot_ids", []),
            "source_organism": organism.get("ncbi_scientific_name") or organism.get("scientific_name"),
            "mutation_count": entity_poly["rcsb_mutation_count"],
            "modeled_residue_count_all_chains": info["deposited_modeled_polymer_monomer_count"],
            "unmodeled_residue_count_all_chains": info["deposited_unmodeled_polymer_monomer_count"],
            "chain_entity_instance_count": polymer_entity["pdbx_number_of_molecules"],
            "assembly_1_details": assembly,
            "nonpolymer_entity_count": info["nonpolymer_entity_count"],
            "validation": {
                "clashscore": geometry["clashscore"],
                "ramachandran_outliers_percent": geometry["percent_ramachandran_outliers"],
                "rotamer_outliers_percent": geometry["percent_rotamer_outliers"],
                "bond_rmsz": geometry["bonds_RMSZ"], "angle_rmsz": geometry["angles_RMSZ"],
                "data_anisotropy": diffraction["data_anisotropy"],
                "real_space_rsrz_outliers_percent": diffraction["percent_RSRZ_outliers"]
            },
            "coordinate_file": f"{pdb_id}.cif",
            "rcsb_url": f"https://www.rcsb.org/structure/{pdb_id}"
        }
    (BASE / "structure_inventory.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({pdb: {"resolution": v["resolution_angstrom"], "chains": v["chain_ids"],
                            "clashscore": v["validation"]["clashscore"], "r_free": v["r_free"]}
                      for pdb, v in result.items()}, ensure_ascii=True))


if __name__ == "__main__":
    main()
