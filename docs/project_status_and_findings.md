# AMP project: current findings and completion status

## Project question and current scope

The initial plan proposes an AI-assisted end-to-end discovery pipeline for
antimicrobial peptides against MDR *Acinetobacter baumannii*, including
activity prediction, sequence generation, safety assessment, and FtsZ docking.
The work completed so far is best framed as a reproducible **data audit and
computational feasibility study with a source-matched published-lead
benchmark**. It does not yet establish a novel AMP or a validated design model.

## Data and preparation

- Acquired the author-exported AbAMPdb tables and all 1,126 DBAASP v4 detail
  records from its *A. baumannii* target search. Source retrieval, licenses,
  citations, and checksums are recorded in `data/raw/dataset_manifest.json`.
- Extracted 2,134 DBAASP activity observations, including 1,834 MIC rows. The
  parsed MIC table preserves ranges, censored values, uncertainty, and
  transition notation rather than converting them to guessed point values.
- AbAMPdb's 250-row table has no MIC unit column, so its values were not pooled
  with DBAASP quantitative MIC data. Exact sequences were checked and sequence
  descriptors were calculated for standard amino acids.
- A note-only resistant-isolate review queue contains 612 MIC records, but
  note-only selection misses paper-supported MDR records with blank or generic
  isolate notes. A separate three-paper source review produced a 101-row
  paper-linked panel; row-level isolate mappings remain to be verified.
- Extracted 1,646 DBAASP hemolysis/cytotoxicity observations across 790 peptide
  IDs. A human erythrocyte 50% hemolysis review queue contains 218 records; 94
  exact numeric records have standard-amino-acid sequences, split across µM
  and µg/mL. No universal safe/unsafe threshold was applied.
- Joined activity and hemolysis evidence by peptide ID. In this snapshot, 208
  peptide IDs have both evidence types and 73 have at least one exact MIC and
  exact H50 result with standard sequence text. Raw records remain distinct.

Full details: `docs/data_audit.md`, `docs/resistant_paper_review.md`, and
`docs/toxicity_data_audit.md`.

## Activity prediction feasibility

An exploratory baseline compared a median predictor, Ridge, RBF-SVR, and
Random Forest on separate µM and µg/mL exact-point MIC cohorts. Splits grouped
similar sequences and, in a sensitivity analysis, candidate publication IDs.
Under the strictest current grouping, Random Forest results were:

| Unit | MAE (log2 MIC) | R² | Interpretation |
|---|---:|---:|---|
| µM | 1.40 | 0.16 | Modest; not a reliable selected predictor |
| µg/mL | 2.26 | -0.44 | Worse than the fold test-set mean baseline in R² |

Performance was more favorable under exact-sequence splits, indicating that
ordinary splits can overstate generalization. The model excludes many
range/censored observations and lacks consistent strain/assay descriptors. It
does not establish MDR-specific prediction. Full protocol and scores:
`docs/activity_baseline.md` and `results/activity_predictions/`.

A second classification endpoint was implemented as a **relative potency
feasibility analysis** because the downloaded sources contain no validated
inactive labels. It classifies lower-than-cohort-median versus at-or-above-
median MIC, separately by unit. With five repeated grouped splits combining
80% edit-identity clusters and shared candidate publications, the best tested
µM result was RBF-SVM ROC-AUC 0.645 (SD 0.093); the best µg/mL result was
Random Forest ROC-AUC 0.569 (SD 0.077). These modest exploratory results do
not support an active/inactive call or candidate-level prediction. See
`docs/activity_classification.md` and `results/activity_classification/`.

## Published lead benchmark

The strongest integrated source-matched evidence currently available is a
four-variant B2RP panel linked to PMID 19793185. Among the named NM isolates,
the [D4K][K16A][L18K] variant has MIC 3 µM across five isolates and a reported
human-RBC 50% hemolysis concentration above 200 µM. Its conservative
same-study concentration margin is therefore `>66.67`. This makes it a useful
published positive-control lead to reproduce, not a newly designed candidate
or a validated therapeutic index. See `docs/candidate_prioritization.md` and
the [primary paper](https://pubmed.ncbi.nlm.nih.gov/19793185/).

A source-informed analogue hypothesis generator now proposes five Arg
substitutions at the experimentally varied B2RP positions 4 and 18; the four
published peptides remain explicit controls. All five proposed analogues were
absent from the local 981-sequence exact-match panel. This is not global
novelty and carries no activity/safety prediction; details are in
`docs/sequence_generation.md`.

## FtsZ target and docking

The specified FtsZ structures 22HQ (GDP-bound, 2.30 Å) and 22HR (apo, 3.00 Å)
were verified and downloaded. Structure state, validation concerns, and
intracellular-target caveats are documented in `docs/ftsZ_structure_audit.md`.
A workstation PyMOL run extracted chain-A inspection PDBs and aligned the
structures at a 0.486 Å outlier-refined C-alpha RMSD across 237 atom pairs.
A peptide-compatible global docking method and controls are documented in
`docs/docking_workflow.md`. The local CABS-flex standalone tool is unavailable
because WSL returned `Wsl/AccessDenied` and Docker's Linux engine is not
running. Ten correct candidate/control sequences have instead been submitted
to the CABS-dock web server against apo 22HR chain A using its default 50-cycle
protocol. All ten now return representative models and cluster tables. These
one-run apo-state outputs remain exploratory: cluster density and RMSD are not
affinity scores, and target binding and bacterial uptake remain hypotheses.
Representative model 1 was downloaded and sequence-checked for each candidate;
heavy-atom contact summaries were prepared across the panel. The ten-model
ensemble for the published lead control also has a contact summary. Full visual
interface review, multi-seed replication, and benchmark-complex validation
remain.

Official sequence-only predictions have now been run for the ten-sequence
panel: ToxinPred3 Hybrid ET+MERCI (threshold 0.38), HemoPI2 HC50 regression
(Random Forest Regressor), and HemoPI2 classification (Hybrid1 ESM2-t6+MERCI,
threshold 0.58). Raw outputs, job IDs, and imported rows are recorded in
`docs/external_predictor_and_docking_runs.md` and
`results/candidate_prioritization/external_predictions/`. These are model
predictions, not experimental measurements or safety clearance. Hemolytik2 is a curated
experimental-data resource, not the hemolysis predictor; HemoPI2 is the
separate classifier/regressor for hemolysis. The project therefore retains
source-reported assays separately and does not fabricate prediction scores.
An additional local H50 regression feasibility run found negative grouped
holdout R² and no improvement over the median baseline; exact feature matches
include 54 µM sequences and 16 µg/mL sequences. This is not a useful candidate
safety model. See `docs/toxicity_data_audit.md`,
`docs/h50_regression_feasibility.md`, and `docs/remaining_tools_setup.md`.

## What is and is not complete

Completed: dataset acquisition and preparation, descriptive safety-data
extraction, exploratory MIC regression, relative-potency classification, and
an assay-specific H50 regression feasibility analysis; three online safety
prediction runs for the labeled candidate panel; focused source-paper
cross-checks; structure verification and chain extraction; a source-matched
published-lead benchmark; an auditable source-informed sequence hypothesis
panel; and an offline peptide lookup UI.

Not completed: a validated MDR-specific predictor; an ML-trained sequence
design model; systematic review of all returned CABS-dock poses; a benchmarked
FtsZ docking experiment; and experimental confirmation. The
current data do not
support claiming that a novel peptide has been designed or that an existing
candidate's efficacy was experimentally reproduced here.

## Reproduction

Follow the ordered commands in the root `README.md`. The data extraction,
cohort-building, and benchmark scripts are deterministic against the saved
source snapshots. The exploratory model has fixed seeds and is described in
its own methods note. Raw files are checksummed in the dataset manifest.

## Recommended research claim

The supportable claim at this stage is that public peptide activity and
toxicity data for *A. baumannii* are heterogeneous and sparse for
publication-independent MDR prediction; similarity/publication-aware
validation weakens sequence-only model performance; and source-linked review
can recover promising published AMP analogues for reproducible positive
controls. The work is a dataset and feasibility study, not an AI-designed
therapeutic discovery result.
