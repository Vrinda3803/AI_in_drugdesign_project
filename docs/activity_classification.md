# Relative MIC classification feasibility

## Endpoint and cohort

This is an exploratory two-class task over exact, parseable, species-level
*Acinetobacter baumannii* MIC observations. MIC values are grouped by exact
sequence and summarized by their median, with µM and µg/mL kept separate. The
label is lower-than-cohort-median concentration versus at-or-above the median
(4 µM for the 255-sequence µM cohort; 12.5 µg/mL for the 256-sequence mass
concentration cohort). **These are relative potency bins, not active/inactive
labels.** The source set has no validated inactive controls.

## Evaluation

Four models (majority dummy, class-weighted logistic regression, RBF SVM, and
class-weighted random forest) were evaluated over five fixed-seed 80/20
GroupShuffleSplit holdouts. Connected components of the explicit 80% global
edit-identity screen were joined when sequences shared a candidate PubMed ID.
There are 70 grouped components in the µM cohort and 69 in the µg/mL cohort.
Mean and standard deviation are across repeated, overlapping holdouts; this is
not an independent external benchmark.

| Unit | Model | Accuracy | F1 | ROC-AUC | PR-AUC |
|---|---|---:|---:|---:|---:|
| µM | RBF SVM | 0.618 ± 0.089 | 0.654 ± 0.134 | 0.645 ± 0.093 | 0.668 ± 0.124 |
| µg/mL | Random forest | 0.580 ± 0.104 | 0.569 ± 0.034 | 0.569 ± 0.077 | 0.580 ± 0.109 |

Performance is modest and variable. The unit-specific cohort median is dataset-
dependent; changing the source mix changes both labels and apparent scores.
The classifier is not used to assign a clinical or biological activity call,
score generated hypotheses, or claim MDR-specific generalization. Full fold
metrics and machine-readable summary are in `results/activity_classification/`.

## Reproduction

From the project root, run:

```powershell
.\.venv\Scripts\python.exe scripts\activity_prediction\mic_relative_classifier.py
```

The script reuses the MIC target construction, sequence descriptors, and
similarity/publication grouping functions from
`scripts/activity_prediction/mic_regression_baselines.py`.
