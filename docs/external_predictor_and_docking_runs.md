# External predictor and docking runs

## Safety predictors

The ten-sequence candidate/control FASTA panel was submitted to the official online prediction forms on 2026-11-09. Results are in the raw tables in `results/candidate_prioritization/external_predictions/` and normalized provenance rows in `results/candidate_prioritization/external_safety_predictions.csv`.

| Tool / endpoint | Job | Model | Threshold | Records |
|---|---:|---|---:|---:|
| ToxinPred3 toxicity classification | 67278 | Hybrid ET+MERCI | 0.38 | 10 |
| HemoPI2 HC50 regression | 78862 | Random Forest Regressor | n/a | 10 |
| HemoPI2 hemolysis classification | 96986 | Hybrid1 ESM2-t6+MERCI | 0.58 | 10 |

These are sequence-only model predictions, not experimental safety measurements or clinical safety decisions. Terminal chemistry is not represented in the submitted one-letter sequences. Scores and endpoint names remain as returned by each server.

## FtsZ peptide docking

Ten correct candidate/control sequences were submitted to the CABS-dock web server against apo *A. baumannii* FtsZ PDB 22HR, chain A, with the default 50 simulation cycles. All ten job pages now return representative models and clustering details. Job links, exact sequences, and top-cluster summaries are in `results/docking/cabsdock/web_jobs.csv`. CABS-dock does not provide a single binding-affinity score for these jobs; its cluster density, RMSD, and model counts are method-specific structural summaries, not binding free energies or experimental affinities.

For the published positive-control lead shown in the UI (`GIWKTIKSMGKVFAGAIKQNL`), cluster 1 density is 26.8454, average RMSD 2.57027, maximum RMSD 32.0078, and 69 members. Full top-cluster summaries for all ten jobs are recorded in the CSV. These are single-run docking outputs, not evidence of FtsZ binding or whole-cell activity. Runs have not been replicated across random seeds or compared against a docking benchmark complex. Representative poses and full cluster tables remain available through each linked CABS-dock page.

The ten downloaded PDB models for positive-control job `3682c266f3d6adc` are
stored under `results/docking/cabsdock/web_models/3682c266f3d6adc/`. Chain A
matches the prepared 22HR target sequence and chain D matches the submitted
21-residue peptide. A reproducible heavy-atom contact summary across the ten
poses is in `docs/cabsdock_pose_contact_review.md` and
`results/docking/cabsdock/web_model_contact_summary_3682c266f3d6adc.csv`.
Contact counts describe geometry only; model occurrence rank and contact counts
are not affinity rankings.

The top representative model (`model_1`) has also been downloaded and
sequence-checked for each of the other nine jobs. The comparable first-pass
contact table for all ten candidates is in
`docs/cabsdock_panel_model1_summary.md` and
`results/docking/cabsdock/web_panel_model1_contact_summary.csv`.

Two additional jobs were submitted in error with 20-residue sequences missing Lys at position 16 (`GIWDTIKSMGKVFAGIRQNL` and `GIWKTIKSMGKVFAGIRQNL`). They are not members of the project FASTA panel and are excluded from the UI and results. The CABS-dock job pages showed them already running; no cancellation control was available. Correct 21-residue panel sequences were submitted in separate jobs and are the only jobs listed in the UI inventory.
