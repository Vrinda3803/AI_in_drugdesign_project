# AMP Design Workbench

An offline-friendly browser UI for selecting the supported bacterium and
generating a source-informed candidate from the project's five B2RP analogue
hypotheses. It also supports searching records and shows the candidate sequence,
source-reported *A. baumannii* MIC evidence, FtsZ docking status, and toxicity
evidence.

## Run for the presentation

First rebuild the local lookup snapshot from the project root:

```powershell
.\.venv\Scripts\python.exe ui\build_ui_data.py
```

Then open `ui/index.html` in a browser. Choose *Acinetobacter baumannii* and
select **Generate AMP sequence** to cycle through the five source-informed
hypotheses. This selects candidates from the reproducible local design panel;
it does not imply that a trained generative model or validated activity model
produced them. The UI is self-contained for offline use; no Streamlit or
external packages are required. Use search to browse published records and
controls. If the browser blocks local files, run
`python -m http.server 8765 --directory ui` from the project root and open
<http://localhost:8765/>.

## What the fields mean

- **AMP sequence:** sequence recorded in DBAASP, AbAMPdb, or the project's
  source-informed hypothesis panel. Terminal modifications are shown where
  available but are not encoded in the one-letter sequence.
- **FtsZ docking:** ten exact candidate/control sequences were submitted to
  the CABS-dock web server against apo FtsZ 22HR chain A. All ten jobs now
  return representative models and cluster tables. The UI shows a cluster
  top-cluster density from the result inventory; this is a method-specific
  clustering metric, not an affinity score. Every record links to the full
  server results. A pose does not establish target binding or antibacterial
  activity.
- **Generated sequence:** the five hypotheses are source-informed residue
  substitutions at positions studied in a published B2RP series. They have no
  measured MIC results. External sequence-only toxicity predictions are shown
  separately from experimental records. Docking results are exploratory and
  do not substitute for MIC or toxicity measurements.
- **Antibacterial activity:** source-reported DBAASP MIC records for *A.
  baumannii* show the isolate/target, raw concentration and unit, medium,
  inoculum where supplied, and candidate publication reference. Values are not
  converted or pooled across assay conditions. Missing evidence is unknown.
- **Toxicity evidence:** source-reported hemolysis/cytotoxicity rows and their
  raw values, units, target cell, and candidate publication reference. These
  are not a universal safe/unsafe label. No toxicity row means unknown.
- **Model predictions:** official ToxinPred3 and HemoPI2 results have been
  imported. Their endpoint, prediction, score, model, threshold, and run date
  appear in a separate table above experimental observations. HC50 values are
  model estimates in µM; they are not experimental results or safety clearance.

The data snapshot is built from downloaded local files by
`build_ui_data.py`. If source datasets change, rerun that script before starting
the server. Source sequence casing is preserved because it can encode
stereochemistry or source-specific residue notation; sequence matching is
case-insensitive only for user search, not for source aliasing or plan status.

To import official candidate predictor outputs, use
`scripts/activity_prediction/import_external_safety_predictions.py` as
described in `docs/remaining_tools_setup.md`. It validates the FASTA candidate
ID and exact sequence, preserves model/threshold/date, then add its normalized
rows to `results/candidate_prioritization/external_safety_predictions.csv`.
Rebuild the snapshot after import.
