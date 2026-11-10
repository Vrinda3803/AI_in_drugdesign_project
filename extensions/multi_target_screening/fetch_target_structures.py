"""Acquire reference FtsZ structures for an expanded bacterial target panel.

Includes experimentally determined A. baumannii, E. coli and S. aureus
structures plus precomputed AlphaFold DB monomer models for E. coli/S. aureus.
These are target structures, not de novo predictions of the AMP candidates.
"""
from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "structures"
EXPERIMENTAL = {"A_baumannii_apo_22HR": "22HR", "E_coli_FtsZ_6LL6": "6LL6", "S_aureus_apo_3VPA": "3VPA"}
AFDB_ACCESSIONS = {"E_coli_FtsZ_AFDB": "P0A9A6", "S_aureus_FtsZ_AFDB": "P0A030"}


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "AMP-Multitarget-structure-panel/1.0"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
                if response.status != 200 or not data:
                    raise ValueError(f"Empty/non-200 response at {url}")
                return data
        except (urllib.error.URLError, TimeoutError, ValueError):
            if attempt == 3: raise
            time.sleep(2 ** attempt)
    raise RuntimeError(url)


def save(name: str, url: str, suffix: str) -> dict:
    path = OUT / f"{name}{suffix}"
    data = fetch(url)
    path.write_bytes(data)
    return {"file": str(path.relative_to(ROOT)), "url": url, "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    artifacts = {}
    for name, pdb_id in EXPERIMENTAL.items():
        artifacts[name] = save(name, f"https://files.rcsb.org/download/{pdb_id}.cif", ".cif")
    afdb_meta = {}
    for name, accession in AFDB_ACCESSIONS.items():
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
        payload = json.loads(fetch(api_url))
        if not payload or not payload[0].get("pdbUrl"):
            raise ValueError(f"No AlphaFold DB PDB URL for {accession}")
        entry = payload[0]
        structure = save(name, entry["pdbUrl"], ".pdb")
        artifacts[name] = structure
        afdb_meta[name] = {"accession": accession, "entry_id": entry.get("entryId"),
                           "model_version": entry.get("latestVersion"), "global_metric_value": entry.get("globalMetricValue"),
                           "global_metric_name": entry.get("globalMetricName"), "api_url": api_url}
    manifest = {"retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "sources": {"experimental": "RCSB PDB official coordinate service", "predicted": "AlphaFold Protein Structure Database API; precomputed AlphaFold models"},
                "artifacts": artifacts, "alphafold_db_metadata": afdb_meta,
                "limits": ["AFDB models are precomputed single-chain FtsZ predictions, not team-generated predictions.",
                           "The AFDB models are target structures, not structures of the peptide candidates.",
                           "FtsZ is intracellular; structural models do not establish peptide uptake, target engagement or antibacterial mechanism."]}
    (OUT / "structure_panel_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"structures": artifacts, "alphafold_db_metadata": afdb_meta, "limits": manifest["limits"]}, indent=2))


if __name__ == "__main__": main()
