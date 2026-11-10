"""Fetch the two project-specified A. baumannii FtsZ entries from RCSB PDB."""
from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "external" / "ftsZ"
PDB_IDS = ("22HQ", "22HR")
BASE = "https://files.rcsb.org/download/{pdb_id}.cif"
API = "https://data.rcsb.org/rest/v1/core/{resource}/{pdb_id}{suffix}"


def fetch(url: str, attempts: int = 4) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "AMP-FtsZ-structure-acquisition/1.0"})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
                if response.status != 200 or not data:
                    raise ValueError(f"Empty/unexpected response: {response.status}")
                return data
        except (urllib.error.URLError, TimeoutError, ValueError):
            if attempt == attempts - 1:
                raise
            time.sleep(2**attempt)
    raise RuntimeError(f"Unable to retrieve {url}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    downloaded = []
    for pdb_id in PDB_IDS:
        coord_path = OUT / f"{pdb_id}.cif"
        coordinate = coord_path.read_bytes() if coord_path.exists() else b""
        if not coordinate or not coordinate.decode("utf-8", errors="ignore")[:20].casefold().startswith(f"data_{pdb_id}".casefold()):
            coordinate = fetch(BASE.format(pdb_id=pdb_id))
        text = coordinate.decode("utf-8", errors="strict")
        if not text[:20].casefold().startswith(f"data_{pdb_id}".casefold()):
            raise ValueError(f"Coordinate file identifier did not match {pdb_id}")
        meta_path = OUT / f"{pdb_id}_rcsb_metadata.json"
        try:
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
            if not all(name in metadata and "retrieval_error" not in metadata[name]
                       for name in ("entry", "polymer_entity_1", "assembly_1")):
                raise ValueError("Incomplete cached RCSB metadata")
        except (OSError, json.JSONDecodeError, ValueError):
            metadata = {}
            for name, resource, suffix in (("entry", "entry", ""), ("polymer_entity_1", "polymer_entity", "/1"),
                                           ("assembly_1", "assembly", "/1")):
                url = API.format(resource=resource, pdb_id=pdb_id, suffix=suffix)
                payload = json.loads(fetch(url))
                metadata[name] = payload
            meta_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
        downloaded.append({"pdb_id": pdb_id, "coordinate_file": str(coord_path.relative_to(ROOT)),
                           "coordinate_bytes": len(coordinate), "coordinate_sha256": hashlib.sha256(coordinate).hexdigest(),
                           "metadata_file": str(meta_path.relative_to(ROOT))})
        print(f"{pdb_id}: saved {len(coordinate)} coordinate bytes")
    manifest = {"retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "source": "RCSB PDB official coordinate download and Data API",
                "entries": downloaded}
    (OUT / "structure_download_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
