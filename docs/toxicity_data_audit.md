# DBAASP hemolysis and cytotoxicity data audit

## Extraction

`scripts/preprocessing/extract_toxicity_records.py` extracts every
`hemoliticCytotoxicActivities` observation from the 1,126 cached DBAASP v4
detail records in the *A. baumannii* target search. The API field's spelling is
retained in the script; the resulting file uses the clearer term “toxicity.”
The extraction preserves peptide identity, termini/modification fields,
target-cell text, the raw lysis/cytotoxicity group and value, concentration,
unit, assay conditions, notes, reference index, and candidate article mapping.

The snapshot contains **1,646 observations across 790 peptide IDs**. Human
erythrocytes are the most frequent target-cell type (829 observations), followed
by HEK293T cells (111) and mouse erythrocytes (99). The measurement labels
include percentages of hemolysis, cytotoxicity, and killing as well as LC50,
IC50, and other endpoints. At least 159 observations lack a measurement-group
label; units include 835 µM, 652 µg/mL, and 159 missing units. See the full
machine-generated counts in `data/processed/toxicity_data_report.json`.

## Interpretation limits

These are not interchangeable safety measurements. Cell type, species,
concentration, unit, endpoint, exposure time (often absent in the API record),
and source protocol all matter. A reported percentage at one concentration is
not an LC50, and erythrocyte hemolysis is not general mammalian cytotoxicity.
The values are preserved verbatim; no binary “safe” label, unit conversion,
or pooled toxicity score is generated.

Reference-to-publication assignment is retained as a candidate 1-based article
index mapping and still needs record/source-paper validation. Do not train a
single toxicity predictor on the complete extracted table. A defensible next
step is to select one endpoint and compatible cell type (for example, human
erythrocyte hemolysis with a clearly interpretable concentration/response),
review individual source papers, and only then define a modeling cohort and
publication/sequence-grouped validation. Missing toxicity data is unknown,
not nontoxicity.

## Human erythrocyte 50% hemolysis candidate endpoint

As a more focused source-review queue, `scripts/preprocessing/prepare_human_rbc_h50_candidates.py`
selects records where the target is human erythrocytes and the reported value
is exactly “50% Hemolysis.” This yields 218 observations across 165 sequences;
143 are standard-amino-acid sequences. Of these, 118 concentrations are
syntactically exact numeric values, 80 are censored, and 20 include reported
uncertainty notation. Among exact values with standard sequences, 94 rows span
both µM (63 rows) and µg/mL (31 rows), so they remain separate. This is a
candidate assay-specific H50 endpoint, not a standardized or validated
cytotoxicity model; publication and protocol review is still required.

An exploratory unit-specific H50 regression was subsequently run on exact
sequence-feature matches. It found negative grouped holdout R² and no
improvement over the median dummy. Case-sensitive sequences without feature
rows were retained in the source target table and excluded from model fitting;
they were not uppercased because case may encode stereochemistry. See
`docs/h50_regression_feasibility.md` and `results/safety_predictions/`.

The filtered raw-plus-parsed records are in
`data/processed/dbaasp_human_rbc_h50_candidate_records.csv`; counts are in
`data/processed/human_rbc_h50_candidate_report.json`. Exact/censored values
remain distinct and there is no safe/unsafe threshold.

## Activity and hemolysis overlap

`scripts/preprocessing/build_activity_toxicity_evidence_panel.py` joins the
DBAASP MIC and human-RBC H50 candidate observations by DBAASP peptide ID and
stores all raw observations as JSON lists in
`data/processed/dbaasp_activity_toxicity_evidence_panel.csv`. In the current
snapshot, 208 peptide IDs have both evidence types; 73 have at least one exact
MIC and one exact H50 observation with standard-amino-acid sequence text. The
report also counts same-unit MIC/H50 observation pairs, but those are
cross-products of observations, not independent peptide sample sizes, and do
not guarantee the assays share conditions. No therapeutic index is calculated.

## Reproduce

```powershell
.\.venv\Scripts\python.exe scripts\preprocessing\extract_toxicity_records.py
.\.venv\Scripts\python.exe scripts\preprocessing\prepare_human_rbc_h50_candidates.py
.\.venv\Scripts\python.exe scripts\preprocessing\build_activity_toxicity_evidence_panel.py
```
