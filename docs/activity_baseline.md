# Exploratory MIC regression baseline

## Endpoint

This is an initial regression experiment, not a final potency model. It uses only DBAASP MIC observations where the target name is exactly `Acinetobacter baumannii` (no strain name is reported), the sequence contains only the 20 standard amino acids, and the concentration can be parsed as an exact point. Named strains are excluded here because assay/strain context is not yet represented consistently. The endpoint is the median log2(MIC) for each exact sequence within each source unit; µM and µg/ml are trained and predicted separately. Exact observations for the same peptide are median-aggregated, while their DBAASP record IDs and candidate publication IDs stay in the target table.

This cohort does **not** establish MDR-specific prediction: the target field does not identify a resistant isolate. It also excludes range, censored, uncertainty, and transition observations from this simple point regressor; those rows remain in the parsed source data.

## Cohort sizes

| Unit | Source observations | Unique valid sequences | Sequences with repeat observations |
|---|---:|---:|---:|
| µM | 271 | 255 | 14 |
| µg/ml | 292 | 256 | 11 |

## Evaluation

Four small-data baselines were compared: median dummy, Ridge, RBF-SVR, and Random Forest. Scores are means and standard deviations over five overlapping 80/20 group holdouts (seed 42). Targets and predictions are log2 concentrations in the stated source unit. The error factor is `2^MAE`, a multiplicative scale interpretation of log2-MAE.

Three grouping protocols are reported:

1. Exact sequence grouping (a reference point; each sequence has already been reduced to one row).
2. Connected components with global edit identity of at least 80%, where identity is `(max_length - Levenshtein distance) / max_length`.
3. The same sequence components joined when records share a *candidate* PubMed ID. DBAASP activity references are numeric indices into each peptide response's ordered article list in most records; this mapping is retained as a candidate and needs source-level confirmation. It is a leakage sensitivity analysis, not treated as verified bibliographic truth.

Full fold-level scores and predictions are in `results/activity_predictions/`. Summary results below compare the median dummy and Random Forest on the stricter two grouping protocols.

| Unit | Grouping | Median dummy MAE (log2) / R² | Random Forest MAE (log2) / R² |
|---|---|---:|---:|
| µM | 80% identity | 1.51 / -0.04 | 1.26 / 0.25 |
| µM | 80% identity + candidate publication | 1.60 / -0.12 | 1.40 / 0.16 |
| µg/ml | 80% identity | 2.43 / -0.28 | 2.28 / -0.06 |
| µg/ml | 80% identity + candidate publication | 2.25 / -0.21 | 2.26 / -0.44 |

## Interpretation and limits

The exact-sequence split looked more favorable than the similarity-controlled split, so ordinary random peptide splits would overstate generalization. Adding candidate publication grouping gives weaker and variable performance. The Random Forest is only a provisional comparator; it is not selected as the project model. In µg/ml, cross-validated R² remained below zero in both stricter protocols. The µM result is modest and still needs publication-level reference verification, an independent test set, a better-supported sequence clustering method, and assay-context review.

The present models use sequence descriptors only, do not encode assay medium or bacterial strain, and are trained on species-level records with no reported strain. They should not be applied as evidence that a candidate is active against MDR *A. baumannii*. No binary activity labels were created.

## Reproduce

Install the pinned dependencies into the project virtual environment and run:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe scripts\activity_prediction\mic_regression_baselines.py
```

The script saves unit-specific full-cohort model objects for future research use, along with their preprocessing pipelines. Those fitted artifacts are not independently validated models.
