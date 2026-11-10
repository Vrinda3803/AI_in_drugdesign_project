# Exploratory human-RBC H50 regression feasibility

## Endpoint and cohort

This local safety-modeling experiment uses exact numeric concentrations from
DBAASP's “50% Hemolysis” endpoint for human erythrocytes. Concentrations are
modelled as log2(H50), with µM and µg/mL kept separate. This is a narrow,
assay-specific hemolysis endpoint; it is not general toxicity, a calibrated
safety classification, or a clinical cutoff.

The source-review cohort has 63 exact µM observations over 61 sequences and 31
exact µg/mL observations over 31 sequences. The feature table has exact,
case-sensitive matches for 54 µM sequences (55 observations) and 16 µg/mL
sequences (16 observations). Seven µM and 15 µg/mL sequences were excluded
because no matching feature row exists. The feature pipeline accepts only
uppercase standard-AA strings. Since sequence case may encode stereochemistry
or another sequence convention, the analysis does not uppercase these strings
automatically.

## Evaluation and result

The median baseline, Ridge, RBF-SVR, and Random Forest were evaluated over five
fixed-seed 80/20 `GroupShuffleSplit` holdouts, grouped by connected 80% global
edit-identity components and candidate PubMed IDs. There are 22 grouped
components in the µM cohort and 12 in the µg/mL cohort. Holdouts overlap, the
publication-index mapping is not comprehensively paper-verified, and there is
no independent external test set.

| Unit | Model | Sequences after feature join | MAE (log2 H50) | R² |
|---|---|---:|---:|---:|
| µM | Median dummy | 54 | 2.03 ± 0.49 | -0.27 ± 0.21 |
| µM | Random forest | 54 | 2.04 ± 0.42 | -0.31 ± 0.29 |
| µg/mL | Median dummy | 16 | 2.20 ± 0.75 | -1.32 ± 0.75 |
| µg/mL | Random forest | 16 | 2.27 ± 0.78 | -1.54 ± 0.91 |

The learned models do not improve over the simple baseline. Their holdout R²
values are negative and errors correspond to roughly four-to-six-fold
concentration error. This shows the local H50 cohort and descriptors do not
support a useful standalone safety regressor. Do not use these models to score
generated hypotheses. Censored, range, and uncertain values remain preserved
in the source data but are excluded from this exact-point baseline. Terminal
modifications are also absent from the sequence-only features.

## Reproduction

```powershell
.\.venv\Scripts\python.exe scripts\activity_prediction\h50_regression_baselines.py
```

Fold metrics, predictions, sequence-level targets, and the summary JSON are in
`results/safety_predictions/`. This is a reproducible feasibility analysis,
not a selected predictor.
