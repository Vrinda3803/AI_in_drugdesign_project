# FtsZ peptide docking plan

## Method choice

Use CABS-flex standalone 3's CABSdock workflow for exploratory peptide–protein
docking. Its documented workflow samples peptide conformations and positions
globally around a receptor when no binding site is known, with optional contact
restraints, clustering, and all-atom reconstruction. This better matches an
unknown AMP binding mode than applying a small-molecule scoring protocol to a
flexible peptide. The original CABS-dock method was designed for blind flexible
peptide docking ([original server paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4489223/));
the current [standalone 3 documentation](https://cabsflex.lcbio.pl/peptide-protein-docking.html)
describes the modern command-line workflow.

HADDOCK3 can be considered only for a separate, explicitly hypothesis-driven
refinement after evidence supports a site. Its peptide guide recommends an
ensemble of peptide conformations, and its protein–peptide example defines
active site residues in advance ([peptide best-practice guide](https://www.bonvinlab.org/haddock3-user-manual/bpg/peptides.html),
[protein–peptide example](https://www.bonvinlab.org/haddock3-user-manual/docking_scenarios/prot-peptide.html)).
The FtsZ nucleotide pocket is known from the structure, but there is no current
evidence that an AMP binds there; imposing it as the AMP site would bake the
desired answer into the protocol. HADDOCK is therefore not the primary screen.

## Receptor and biological assumptions

- Start with chain A from 22HR as the primary apo receptor. Keep the original
  mmCIF untouched and record every cleanup step and residue mapping in derived
  input metadata.
- Use 22HQ chain A as a receptor-state sensitivity run. Its GDP is part of the
  deposited structure; retain it for a global search and report whether poses
  overlap it. Do not silently strip GDP or claim that the GDP-bound pocket is
  available to a peptide.
- These entries are monomeric author-defined assemblies. Results refer only to
  the chosen monomeric receptor model and do not represent FtsZ filament
  interfaces or dynamics.
- Restrict an initial comparable docking panel to standard-amino-acid linear
  peptides supported by the input protocol. PTMs, disulfide bridges, and
  cyclization need explicit encoding and separate validation before inclusion.
- FtsZ is intracellular. A modeled pose does not show peptide uptake, target
  engagement, or that target inhibition explains antibacterial activity.

## Candidate selection and analysis gates

1. Review resistant-isolate records from
   `data/processed/dbaasp_a_baumannii_resistant_isolate_records.csv` against
   their source papers. Keep exact MIC, unit, isolate, and publication attached;
   resolve candidate citation mappings before using them as independent
   evidence. Do not select a peptide solely because an MIC is favorable in a
   different assay context.
2. Define and publish candidate inclusion criteria before docking. Prioritize
   reproducible, explicitly resistant-isolate observations and report coverage
   of units/isolates. Exclude unreviewed or nonstandard sequences from the first
   comparable panel, while retaining them in the source data.
3. Run global docking with multiple recorded random seeds for each candidate
   against both receptor states. Save exact tool version, command, parameters,
   input checksums, logs, and every output model.
4. Compare cluster recurrence and interface geometry across seeds and receptor
   states. Report the full pose/cluster distribution, not just the lowest
   scoring pose. Any score is a method-specific ranking quantity, not a binding
   free energy or experimental affinity.
5. Interpret target docking alongside MIC evidence and cytotoxicity/hemolysis
   evidence. Do not infer that FtsZ binding causes whole-cell activity.

## Validation and limitations

There is no AMP–FtsZ complex in the project's downloaded structures and no
established AMP positive-control ligand in the current data. Thus, the docking
stage cannot claim validated pose accuracy. First test whether the pipeline
runs reproducibly on a known peptide–protein benchmark with a reference complex;
report benchmark protocol/results separately from the FtsZ candidate screen.
If that benchmark cannot be reproduced, label FtsZ docking exploratory and
avoid using docking rank as a decisive candidate-selection feature.

The current workspace is Windows. `wsl.exe --list --quiet` returned
`Wsl/AccessDenied`, and CABS-flex standalone's standard installation guidance
targets Linux, macOS, or WSL2. To start the peptide-protein screen, ten correct
candidate/control sequences were submitted to the official CABS-dock web
server against apo 22HR chain A at its default 50 simulation cycles. All ten
jobs now return representative models and cluster tables; links and inventory
are in `results/docking/cabsdock/web_jobs.csv`. These one-run apo-state jobs do
not replace the prepared 60-run, three-seed, two-state standalone matrix or a
benchmark-complex validation. Cluster density and RMSD are method outputs, not
binding-affinity scores.

The `model_1` representative pose has been downloaded and checked for all ten
jobs, and a geometric contact comparison is recorded in
`docs/cabsdock_panel_model1_summary.md`. The full ten-pose ensemble was also
downloaded and summarized for the published positive-control lead. These
contact counts are descriptors of the returned structures only; they are not
affinity scores. Systematic interface visualization and benchmark validation
remain future work.

Two additional web jobs were mistakenly submitted with sequences missing Lys
at position 16. They are not panel entries, are excluded from the job inventory
and UI, and were already running without a visible cancellation control. The
exact intended 21-residue sequences were submitted separately and are the only
jobs listed in the project inventory.

PyMOL is installed and has now been used to export chain-A inspection PDBs at
`data/external/ftsZ/prepared/` and record an apo/GDP-bound alignment. These
exports improve readiness of the structural inputs, but do not remove the
Linux docking requirement or count as peptide docking.

## Prepared run matrix

`scripts/docking/prepare_cabsdock_run_matrix.py` creates a reproducible command
matrix and Bash runner for the ten-entry candidate/control FASTA panel. It
plans three recorded random seeds against both the 22HR apo and 22HQ
GDP-bound chain-A receptors (60 runs total), with the documented CABSdock
sampling options `-a 20 -y 50 -s 200`, saved config/contact/simulation files,
and a distinct output directory per run. The script only prepares commands;
it does not install or execute CABSdock. The settings follow the official
[CABS-flex examples](https://cabsflex.lcbio.pl/examples.html) and
[command-line reference](https://cabsflex.lcbio.pl/options.html), but have not
been tuned or benchmarked for FtsZ. Run the matrix only after reviewing compute
availability and the resulting protocol. These remain exploratory runs; any
scores must be interpreted as method-specific rankings and summarized across
seeds and receptor states.

After installing CABS-flex standalone 3 and activating its environment, from
the project root run:

```bash
python3 scripts/docking/prepare_cabsdock_run_matrix.py
bash scripts/docking/run_cabsdock_panel.sh
python3 scripts/docking/summarize_cabsdock_runs.py
```

The generated plan, run matrix, shell script, run inventory, and summary are
written under `results/docking/cabsdock/` and `scripts/docking/`. The collector
requires a success marker in `CABS.log`, at least one ranked `model_N.pdb`, and
a nonempty `output_data/cluster_names.txt` before calling a run complete. It
records artifact checksums and pose-file counts without treating model rank or
cluster density as a binding-affinity score. The official CABSdock output guide
describes ranked medoid models, cluster density, and interaction-energy plots
as separate outputs.

When the UI snapshot is rebuilt, it reflects complete or partial run artifacts
for exact one-letter sequences from the panel. It continues to show a numeric
docking score as unavailable unless a separately validated numeric value is
made available; terminal modifications are not represented by the docking
sequence input.

## Sources

- [CABS-dock blind peptide docking paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4489223/)
- [CABS-flex standalone 3 docking guide](https://cabsflex.lcbio.pl/peptide-protein-docking.html)
- [CABS-flex standalone 3 repository and installation notes](https://github.com/LCBio/CABSflex_standalone)
- [HADDOCK3 peptide best-practice guide](https://www.bonvinlab.org/haddock3-user-manual/bpg/peptides.html)
- [Project FtsZ structure audit](ftsZ_structure_audit.md)
