# AI-assisted AMP prioritization for *A. baumannii*

This repository is a reproducible data and computational feasibility pipeline. It covers source-data acquisition, activity and toxicity extraction, MIC parsing, sequence descriptors, exploratory unit-specific MIC regression and relative-potency classification, resistant-isolate paper review, source-informed analogue hypotheses, FtsZ structure preparation, external sequence-only toxicity predictions, CABS-dock pose clustering, and an offline peptide workbench. Neither activity model is a reliable validated predictor; the analogue generator is not a trained sequence model. Ten CABS-dock jobs have returned models and cluster summaries; these are exploratory outputs, not affinity scores or experimental validation.

The completed data-driven feasibility study, its supported conclusions, and its remaining limitations are summarized in `docs/project_status_and_findings.md`.
The item-by-item comparison with the original workflow is in
`docs/planned_work_completion_matrix.md`; the lab team handoff is in
`docs/experimental_validation_handoff.md`.
The updated mid-review presentation, with speaker notes and the latest extension
results, is in `deliverables/AMP_Project_Update_and_Viva_Presentation_Reviewed.pptx`.
The mid-review completion map for the expanded AMP-generation plan, including
separate AMP/non-AMP and LSTM extensions, multi-species evidence, target
structures, and remaining gaps, is in
[`extensions/PROJECT_COMPLETION_STATUS.md`](extensions/PROJECT_COMPLETION_STATUS.md).

## Current data

- AbAMPdb author exports, including the published 250-row table.
- DBAASP API v4 index and all available detail records from its *A. baumannii* search (1,126 records).
- 2,134 *A. baumannii* activity observations; 1,834 MIC rows across 1,066 DBAASP peptide IDs in the downloaded API snapshot.
- Source, license/citation notes, retrieval date, and SHA-256 checksums: `data/raw/dataset_manifest.json`.
- Detailed audit, limitations, and planned modeling decisions: `docs/data_audit.md`.
- Explicitly labeled resistant-isolate DBAASP subcohort (still source-paper unreviewed): `data/processed/dbaasp_a_baumannii_resistant_isolate_records.csv` and `docs/data_audit.md`.
- Focused source-paper cross-check and 101-row paper-linked review panel: `docs/resistant_paper_review.md` and `data/processed/dbaasp_a_baumannii_reviewed_paper_records.csv`.
- Regression cohort and grouped-CV results/limitations: `docs/activity_baseline.md`.
- Exploratory relative MIC classification and label caveat: `docs/activity_classification.md` and `results/activity_classification/`.
- Exploratory unit-specific human-RBC H50 regression and limitations: `docs/h50_regression_feasibility.md` and `results/safety_predictions/`.
- PDB ID verification, structural quality, and target caveats: `docs/ftsZ_structure_audit.md`.
- PyMOL chain-A extraction and C-alpha alignment of apo/GDP-bound FtsZ: `data/external/ftsZ/prepared/` and `results/structure_comparison/pymol_run.log`.
- Peptide-compatible FtsZ docking rationale and planned controls: `docs/docking_workflow.md`.
- Prepared 60-run, multi-seed CABSdock command matrix (not executed): `scripts/docking/prepare_cabsdock_run_matrix.py`.
- DBAASP hemolysis/cytotoxicity extraction, human erythrocyte H50 candidate cohort, and endpoint limitations: `docs/toxicity_data_audit.md`.
- Per-peptide MIC/hemolysis evidence inventory, retaining observation-level provenance: `data/processed/dbaasp_activity_toxicity_evidence_panel.csv`.
- Reproducible within-paper B2RP positive-control lead prioritization: `docs/candidate_prioritization.md`.
- Source-informed B2RP sequence hypotheses and local novelty audit: `docs/sequence_generation.md`.
- Candidate FASTA panel for official safety predictors and docking, with exact setup instructions: `results/candidate_prioritization/b2rp_candidate_panel.fasta` and `docs/remaining_tools_setup.md`.
- Official ToxinPred3 and HemoPI2 outputs for all ten panel sequences, with settings and raw server tables: `results/candidate_prioritization/external_predictions/` and `docs/external_predictor_and_docking_runs.md`.
- Ten completed CABS-dock web jobs against apo FtsZ 22HR chain A, including returned representative models and cluster tables: `results/docking/cabsdock/web_jobs.csv`. Top-cluster density and RMSD summaries are recorded for all ten; these method-specific outputs are not binding affinities.
- Ten downloaded models for the published-control job and a sequence-checked geometric contact summary: `results/docking/cabsdock/web_models/3682c266f3d6adc/` and `docs/cabsdock_pose_contact_review.md`.
- Model-1 pose comparison and geometric contacts for all ten CABS-dock candidates: `docs/cabsdock_panel_model1_summary.md` and `results/docking/cabsdock/web_panel_model1_contact_summary.csv`.
- Offline AMP Design Workbench for bacterium selection, source-informed sequence generation, MIC evidence, toxicity predictions, and live docking job status: `ui/README.md`.
- Separate LSTM hypothesis set (50 sequences) and its training/validation report: `extensions/lstm_generation/output/`.
- LSTM candidate descriptor and nearest local sequence audit (not an activity/safety model): `extensions/lstm_generation/output/lstm_candidate_audit.csv`.
- Weak-label AMP-vs-putative-non-AMP grouped-CV prototype: `extensions/amp_binary_classifier/output/`.
- Cross-species MIC evidence and grouped baseline results (five species): `extensions/multi_target_screening/output/`.

The raw source values are retained. MIC units are not combined; ranges, censoring, and reported uncertainty stay distinguishable. AbAMPdb MIC values are not pooled because the export does not state units. Missing activity is not treated as inactivity.

## Data acknowledgements

DBAASP records in this repository were obtained from the [DBAASP](https://dbaasp.org), an open-access AMP data resource. DBAASP's [current data use policy](https://www.dbaasp.org/terms-and-conditions) permits redistribution with public acknowledgement. Please cite Pirtskhalava et al., *Nucleic Acids Research* 49(D1), D288–D297 (2021), [doi:10.1093/nar/gkaa991](https://doi.org/10.1093/nar/gkaa991), and review the policy before reusing the records. AbAMPdb source files are accompanied by their MIT license in `data/raw/AbAMPdb/LICENSE`.

## Recreate the current preparation

Use Python 3.11 or newer. The acquisition and preprocessing scripts use the standard library; the regression stage uses pinned scikit-learn dependencies.
Acquisition requires internet access. The CABS-flex standalone application cannot run in the current Windows environment because WSL is denied and Docker's Linux engine is unavailable. CABS-dock web jobs are being used for the flexible peptide-protein screen instead.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe scripts\data_collection\download_dbaasp_ab_details.py
.\.venv\Scripts\python.exe scripts\data_collection\prepare_amp_datasets.py
.\.venv\Scripts\python.exe scripts\preprocessing\parse_mic_observations.py
.\.venv\Scripts\python.exe scripts\features\engineer_sequence_features.py
.\.venv\Scripts\python.exe scripts\preprocessing\extract_resistant_subcohort.py
.\.venv\Scripts\python.exe scripts\preprocessing\extract_paper_review_panel.py
.\.venv\Scripts\python.exe scripts\preprocessing\extract_toxicity_records.py
.\.venv\Scripts\python.exe scripts\preprocessing\prepare_human_rbc_h50_candidates.py
.\.venv\Scripts\python.exe scripts\preprocessing\build_activity_toxicity_evidence_panel.py
.\.venv\Scripts\python.exe scripts\activity_prediction\build_b2rp_lead_benchmark.py
.\.venv\Scripts\python.exe scripts\sequence_generation\generate_b2rp_hypotheses.py
.\.venv\Scripts\python.exe scripts\sequence_generation\export_candidate_fasta.py
.\.venv\Scripts\python.exe scripts\activity_prediction\mic_regression_baselines.py
.\.venv\Scripts\python.exe scripts\activity_prediction\mic_relative_classifier.py
.\.venv\Scripts\python.exe scripts\activity_prediction\h50_regression_baselines.py
.\.venv\Scripts\python.exe scripts\docking\acquire_ftsZ_structures.py
.\.venv\Scripts\python.exe scripts\docking\summarize_ftsZ_inventory.py
.\.venv\Scripts\python.exe scripts\docking\prepare_ftsZ_with_pymol.py
.\.venv\Scripts\python.exe scripts\docking\prepare_cabsdock_run_matrix.py
```

The downloader is resumable and reuses verified detail JSON files. It needs internet access to `dbaasp.org`. The PyMOL preparation command requires PyMOL to be installed; pass `--pymol-exe` if its executable is not auto-detected.

For a one-command rebuild of the local analysis artifacts from the cached
source snapshot, run:

```powershell
.\scripts\run_reproducible_pipeline.ps1
```

The pipeline stops on the first failed stage, refreshes any missing DBAASP
details only with `-RefreshDbaasp`, and rebuilds the UI snapshot only with
`-BuildUi`. Use `-SkipPyMol` if PyMOL is unavailable. It does not submit
sequences externally or start docking; imported predictor results and web-job
links remain recorded separately. Each run writes its stages and the current
external-result inventory to `results/pipeline_run_manifest.json`.

## Open the demo UI

```powershell
.\.venv\Scripts\python.exe ui\build_ui_data.py
```

Open `ui/index.html` in a browser. The interface embeds its local snapshot and requires no Streamlit installation or server. If local-file scripts are blocked, run `python -m http.server 8765 --directory ui` and open <http://localhost:8765/>.

## Research boundary

Current outputs are data preparation and exploratory modeling artifacts, not evidence of efficacy. Before treating an activity predictor as reliable, validate the reference mapping, account for assay/strain context, use an independent publication-level test split, and review the sequence-clustering method. Toxicity, hemolysis, FtsZ docking, and ranking are later computational assessments and require separate validation; none should be represented as experimental activity or clinical effectiveness.
