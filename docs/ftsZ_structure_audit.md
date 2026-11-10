# FtsZ structure audit

## Identifier verification

Both project-specified PDB IDs are valid and correspond to *Acinetobacter baumannii* cell division protein FtsZ, UniProt V5V8B6. RCSB identifies 22HQ as the GDP-bound crystal structure and 22HR as the apo structure. Both are wild type by the deposited mutation count. No substitute structure was needed.

Coordinates and the full RCSB Data API responses are in `data/external/ftsZ/`. Recreate the downloads with `scripts/docking/acquire_ftsZ_structures.py`, then recreate the compact inventory with `scripts/docking/summarize_ftsZ_inventory.py`. SHA-256 checksums and byte counts are in `data/external/ftsZ/structure_download_manifest.json`.

## Structure and validation summary

| PDB | State | Resolution | R-free | FtsZ chains in asymmetric unit | Modeled / unmodeled residues across chains | Clashscore | Ramachandran outliers |
|---|---|---:|---:|---|---:|---:|---:|
| 22HQ | GDP-bound | 2.30 Å | 0.263 | A, B, C, D | 1,191 / 33 | 5.91 | 0.00% |
| 22HR | Apo | 3.00 Å | 0.285 | A, B | 599 / 15 | 10.20 | 0.51% |

Both are experimental X-ray structures, contain the same *A. baumannii* FtsZ sequence family, and are assigned monomeric author-defined biological assemblies in RCSB. The deposited polymer sequences differ by one terminal residue: 22HQ has 306 amino acids and 22HR has 307. 22HQ has better nominal resolution and lower R-free/clashscore, but its validation also reports pseudo-translational symmetry and 16.03% rotamer outliers. 22HR has lower resolution and higher R-free/clashscore, with 18.42% rotamer outliers. The structures should therefore be treated as complementary conformational states rather than interchangeable files.

RCSB reports one nonpolymer ligand entity in 22HQ and none in 22HR; its entry page identifies the 22HQ ligand as GDP. The 22HQ GTPase pocket is consequently occupied by GDP in the deposited coordinates. Receptor preparation must make the ligand treatment explicit and preserve the unmodified downloaded structure as the reference. The author-defined assembly is monomeric although multiple FtsZ copies occur in each asymmetric unit, so later docking should state which chain/assembly representation is used.

## Biological interpretation

The 2026 primary structural study describes the structures as apo and GDP-bound *A. baumannii* FtsZ and reports small-molecule virtual screening at the GTP-binding pocket. That supports the protein identity and provides a known reference site; it does not establish that generated AMPs bind there.

FtsZ is an intracellular division target. Some AMPs can reach intracellular targets, but that is peptide- and condition-dependent. Docking a candidate to FtsZ is therefore a target-binding hypothesis only; candidate entry into bacteria and FtsZ engagement are not shown by the MIC data or by docking. This distinction should remain in the final ranking and project presentation.

## Official records

- [RCSB 22HQ](https://www.rcsb.org/structure/22HQ)
- [RCSB 22HR](https://www.rcsb.org/structure/22HR)
- [Primary study DOI](https://doi.org/10.1016/j.ijbiomac.2026.152295)
- [RCSB entry Data API](https://data.rcsb.org/)

The structural audit is complete. Ten CABS-dock web jobs have since returned representative models and cluster tables against apo 22HR chain A. Systematic pose review, receptor cleanup checks, and interaction interpretation remain future stages. CABS-dock does not report a single affinity score for these jobs; cluster density and RMSD are method-specific outputs and not measured binding energies.

## PyMOL structure inspection and chain extraction

PyMOL 3.1.9 was found on the workstation and used headlessly. The reproducible
driver is `scripts/docking/prepare_ftsZ_with_pymol.py`; its generated PyMOL
command file is `scripts/docking/ftsZ_prepare_and_align.pml`. It exports chain
A from each deposited mmCIF: 22HQ protein plus its GDP, and 22HR protein only.
Other chains and solvent are excluded. The original mmCIF files remain
untouched. These are chain-extracted inspection inputs, not fully prepared
docking receptors; protonation, alternate conformers, missing residues,
terminal chemistry, ligand policy, and force-field preparation have not been
resolved.

PyMOL aligned the chain-A C-alpha atoms with an outlier-refined RMSD of 0.486 Å
for 237 atom pairs. The command log is
`results/structure_comparison/pymol_run.log`; hashes and atom/residue counts
are in `data/external/ftsZ/prepared/pymol_preparation_manifest.json`. PyMOL's
installed evaluation build watermarks rendered figures, so no rendered image
is retained as a research figure. The coordinate exports and alignment log
are unaffected by that display limitation. ChemSketch is installed but is a
chemical drawing application and is not needed for this protein structure
preparation step.
