# CABS-dock pose contact review: published B2RP lead control

## Input verification

Ten downloaded model files for CABS-dock job `3682c266f3d6adc` are stored in
`results/docking/cabsdock/web_models/3682c266f3d6adc/`. Each file contains
chain A with the exact 300-residue sequence in the prepared *A. baumannii*
FtsZ 22HR chain-A structure, and chain D with the exact 21-residue peptide
`GIWKTIKSMGKVFAGAIKQNL`. Per-file SHA-256 checksums and atom counts are in the
JSON and CSV result files.

## First-pass geometric contact summary

The script counts heavy-atom pairs within fixed distance cutoffs between
protein chain A and peptide chain D. “Target residues” and “peptide residues”
are the number of distinct residues with at least one atom contact within 5 Å.
Minimum distance is the closest inter-chain heavy-atom pair.

| CABS model | Minimum distance (Å) | Atom pairs ≤4 Å | Atom pairs ≤5 Å | Target residues ≤5 Å | Peptide residues ≤5 Å |
|---:|---:|---:|---:|---:|---:|
| 1 | 2.653 | 123 | 499 | 39 | 17 |
| 2 | 2.735 | 52 | 241 | 29 | 14 |
| 3 | 2.741 | 122 | 400 | 29 | 17 |
| 4 | 2.808 | 105 | 422 | 32 | 16 |
| 5 | 2.720 | 59 | 258 | 21 | 12 |
| 6 | 2.634 | 95 | 370 | 34 | 17 |
| 7 | 2.625 | 148 | 569 | 35 | 16 |
| 8 | 2.684 | 105 | 416 | 29 | 17 |
| 9 | 2.757 | 90 | 266 | 17 | 15 |
| 10 | 2.740 | 185 | 638 | 40 | 16 |

These are geometric summaries of the returned models, not energies, affinities,
or probabilities of biological binding. CABS model numbers refer to occurrence
in the docking trajectory; the largest contact count does not make model 10 a
better-ranked pose. Chain A is a single apo monomer, so the analysis does not
cover FtsZ filament interfaces, alternate receptor states, peptide uptake, or
whole-cell efficacy. The models have not been validated against a known
peptide–FtsZ complex.

## Reproduce

From the project root run:

```powershell
.\.venv\Scripts\python.exe scripts\docking\analyze_cabsdock_pose_contacts.py
```

Outputs:

- `results/docking/cabsdock/web_model_contact_summary_3682c266f3d6adc.csv`
- `results/docking/cabsdock/web_model_contact_summary_3682c266f3d6adc.json`

The analysis script validates both chain sequences against project inputs before
writing results.
