#!/usr/bin/env python3
"""Prepare chain-A FtsZ inspection inputs and an alignment image with PyMOL."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PYMOL = Path.home() / "AppData/Local/Schrodinger/PyMOL2/Scripts/pymol.exe"
HQ = ROOT / "data/external/ftsZ/22HQ.cif"
HR = ROOT / "data/external/ftsZ/22HR.cif"
OUTDIR = ROOT / "data/external/ftsZ/prepared"
FIGDIR = ROOT / "results/structure_comparison"
PML = ROOT / "scripts/docking/ftsZ_prepare_and_align.pml"
LOG = FIGDIR / "pymol_run.log"
METADATA = OUTDIR / "pymol_preparation_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pymol_binary(requested: str | None) -> Path:
    if requested:
        path = Path(requested).expanduser()
    else:
        found = shutil.which("pymol") or shutil.which("pymol.exe")
        path = Path(found) if found else DEFAULT_PYMOL
    if not path.is_file():
        raise SystemExit(
            "PyMOL executable not found. Pass --pymol-exe with the full path. "
            f"Checked: {path}"
        )
    return path


def write_pml() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    # PyMOL uses relative paths because it is launched with the repository root
    # as its working directory. Keep the saved native inputs separate from the
    # untouched downloaded mmCIF references.
    pml_text = r'''reinitialize
load data/external/ftsZ/22HQ.cif, hq_raw
load data/external/ftsZ/22HR.cif, hr_raw
create hq_chainA_gdp, hq_raw and chain A and (polymer.protein or resn GDP)
create hr_chainA_apo, hr_raw and chain A and polymer.protein
save data/external/ftsZ/prepared/22HQ_chainA_GDPbound.pdb, hq_chainA_gdp
save data/external/ftsZ/prepared/22HR_chainA_apo.pdb, hr_chainA_apo
create hr_aligned, hr_chainA_apo
align hr_aligned and polymer.protein and name CA, hq_chainA_gdp and polymer.protein and name CA, cycles=5
hide everything
show cartoon, hq_chainA_gdp and polymer.protein
show cartoon, hr_aligned and polymer.protein
color marine, hq_chainA_gdp and polymer.protein
color gray70, hr_aligned and polymer.protein
set cartoon_transparency, 0.48, hr_aligned and polymer.protein
quit
'''
    PML.write_text(pml_text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pymol-exe", help="Full path to pymol.exe (auto-detected on this Windows setup).")
    args = parser.parse_args()
    for source in (HQ, HR):
        if not source.is_file():
            raise SystemExit(f"Required downloaded structure is missing: {source}")

    executable = pymol_binary(args.pymol_exe)
    version_run = subprocess.run(
        [str(executable), "-cq", "-d", 'print("PYMOL_VERSION",cmd.get_version()[0])'],
        cwd=ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=60,
    )
    version_match = re.search(r"PYMOL_VERSION\s+([^\s]+)", version_run.stdout)
    if version_run.returncode != 0 or not version_match:
        raise SystemExit("Could not determine the PyMOL version before preparation.")
    write_pml()
    result = subprocess.run(
        [str(executable), "-cq", str(PML)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=300,
    )
    LOG.write_text(result.stdout + "\n--- STDERR ---\n" + result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise SystemExit(f"PyMOL failed with code {result.returncode}; inspect {LOG}")

    hq_pdb = OUTDIR / "22HQ_chainA_GDPbound.pdb"
    hr_pdb = OUTDIR / "22HR_chainA_apo.pdb"
    expected = (hq_pdb, hr_pdb)
    missing = [str(path) for path in expected if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise SystemExit("PyMOL returned without expected output(s): " + ", ".join(missing))

    def pdb_counts(path: Path) -> dict[str, int]:
        atom_lines = [line for line in path.read_text(encoding="ascii", errors="replace").splitlines()
                      if line.startswith(("ATOM  ", "HETATM"))]
        residues = {(line[21:22], line[22:27], line[17:20]) for line in atom_lines}
        return {
            "atom_records": len(atom_lines),
            "residues": len(residues),
            "GDP_residues": sum(1 for _, _, name in residues if name == "GDP"),
        }

    align_lines = [line.strip() for line in (result.stdout + result.stderr).splitlines()
                   if re.search(r"RMSD\s*=|ExecutiveAlign", line, re.IGNORECASE)]
    manifest = {
        "software": "PyMOL",
        "software_version": version_match.group(1),
        "executable": str(executable),
        "source_structures": {
            "22HQ.cif": {"sha256": sha256(HQ), "selection": "chain A protein plus GDP; water/other chains excluded"},
            "22HR.cif": {"sha256": sha256(HR), "selection": "chain A protein; other chains/heterogens excluded"},
        },
        "outputs": {
            hq_pdb.name: {"sha256": sha256(hq_pdb), **pdb_counts(hq_pdb)},
            hr_pdb.name: {"sha256": sha256(hr_pdb), **pdb_counts(hr_pdb)},
        },
        "alignment_log_lines": align_lines,
        "pymol_log": str(LOG.relative_to(ROOT)),
        "reproducible_script": str(PML.relative_to(ROOT)),
        "interpretation": "Prepared inspection/receptor inputs only. This is not peptide docking, binding evidence, or a validated receptor-preparation protocol.",
    }
    METADATA.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
