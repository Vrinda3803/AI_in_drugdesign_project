# Data audit: antimicrobial peptides against *A. baumannii*

## Scope and source files

The workspace initially contained no project documents, datasets, or code. In particular, the project-summary document, review slides, and named preprocessed MIC CSV were not present. The current audit therefore starts from author-published AbAMPdb exports and the public DBAASP v4 API. Source locations, retrieval date, API version, licenses/citation notes, and SHA-256 checksums are in `data/raw/dataset_manifest.json`.

Downloaded full DBAASP peptide detail records are cached by numeric record ID under `data/raw/DBAASP/a_baumannii_details/`. Rerunning `scripts/data_collection/download_dbaasp_ab_details.py` reuses validated cached responses and retries missing records. The script uses only Python's standard library.

## Findings

| Source / view | Observed result |
|---|---:|
| AbAMPdb `ampdb-Final.csv` | 250 source rows; 248 sequences contain only standard amino acids; 237 unique exact valid sequences |
| DBAASP *A. baumannii* API search | 1,126 peptide records; 1,126 detail responses downloaded |
| DBAASP activity rows against *A. baumannii* (including named strains) | 2,134 measurements across assay types |
| DBAASP MIC rows | 1,834 measurement rows for 1,066 unique DBAASP peptide IDs |
| DBAASP unique nonempty sequences among MIC rows | 926 raw sequences; 803 distinct sequences restricted to the 20 standard amino acids |
| Exact unique standard-amino-acid sequences across AbAMPdb and DBAASP | 981 feature rows |

The 2,134 activity count and 1,066 MIC-associated peptide IDs match the prior project notes. The current API returns 1,834 MIC records, five more than the previously reported 1,829; all retrieved records are retained pending comparison with the missing prior processed file. Record IDs and raw concentration strings remain available to resolve this discrepancy.

The source `ampdb-Final.csv` has an `MIC value` field without a unit column. It is not currently suitable for quantitative MIC training or direct pooling with DBAASP. `ampdb-Final3.csv` has no header row and `ampdb-Final2.csv` is a one-record fragment; neither is being interpreted as an independent table. `data.json` was excluded because the public repository fixture contains unrelated account credential hashes and is unnecessary for peptide analysis.

## Outputs created

- `data/processed/abampdb_records_qc.csv`: AbAMPdb source rows with sequence validity flags; original source columns retained.
- `data/processed/dbaasp_a_baumannii_activity_records.csv`: all extracted target activity rows, with raw values, units, strain/isolate notes, medium, conditions, and reference fields.
- `data/processed/dbaasp_a_baumannii_mic_measurements.csv`: MIC-only view of the same raw assay rows. Range/censored concentration strings are not converted or collapsed.
- `data/processed/sequence_provenance_panel.csv`: exact-sequence grouping across sources, retaining record IDs.
- `data/processed/dbaasp_a_baumannii_mic_parsed.csv`: syntax-level concentration parsing while retaining source strings. Current parse outcomes: 1,274 exact decimals, 1 exact comma-decimal, 259 ranges, 267 censored values, 31 values with reported uncertainty, and 2 reported transitions. The uncertainty and transition formats are deliberately not given interval/point semantics.
- `data/processed/sequence_features.csv`: amino-acid composition and sequence-derived physicochemical descriptors for 981 unique valid sequences.
- `data/processed/feature_engineering_metadata.json`: equations, scales, and terminal/modification assumptions for the descriptor table.
- `data/processed/dataset_qc_report.json`: reproducible counts, sequence-QC results, and current limitations.
- `data/processed/dbaasp_a_baumannii_resistant_isolate_records.csv`: records selected using explicit resistance language in the isolate note; this is a review panel, not a harmonized MDR benchmark.
- `data/processed/resistant_isolate_subcohort_report.json`: cohort selection rule and counts.
- `docs/activity_baseline.md`: exploratory, unit-separated regression cohort and similarity/publication-grouped evaluation.

Run the extraction again with:

```powershell
python scripts\data_collection\prepare_amp_datasets.py
python scripts\preprocessing\parse_mic_observations.py
python scripts\features\engineer_sequence_features.py
```

The feature table uses mean Kyte-Doolittle hydrophobicity, an Eisenberg-scale 100-degree alpha-helical hydrophobic moment, and Bjellqvist pKa charge/pI calculations. Molecular weights and charge/pI are theoretical for a linear peptide with free, unmodified termini. These are not measured properties; modification-aware structures are outside this initial table.

## Remaining research decisions

The explicit resistant-isolate extraction is a first curation aid. Review its
source records and papers, harmonize resistance phenotype and assay context,
and keep publication and sequence groups intact before evaluating any model.

1. Parse concentrations only when a numeric value and unit can be interpreted without guessing. Keep exact, range, and censored observations distinguishable; report rows that cannot be normalized.
2. Review assay conditions, organism strain names, and explicit resistance annotations. Do not infer MDR from species identity or from an unlabeled strain.
3. Resolve repeated MIC observations per peptide/strain/assay condition; retain source rows and document any aggregation rule.
4. The first regression baseline groups close sequence neighbors using a global edit identity screen. Compare with a standard protein alignment/cluster tool and improve the validation before interpreting performance.
5. Activity reference numbers were mapped to article-array entries as candidate PubMed IDs. Verify this mapping against DBAASP's semantics before treating publication-grouped CV as definitive.
6. The initial continuous endpoint uses exact point observations only. Add interval/censoring-aware regression methods to incorporate ranges and bounds without converting them to guessed point values.
7. No binary activity threshold has been set. Define one only with defensible biological criteria; do not treat missing activity as inactivity.
8. The exploratory regressors do not model assay medium or strain and are not MDR-specific. Improve assay-context handling and obtain an independent publication-level evaluation before selecting a predictor.

The data-audit stage produced no activity classes or docking claims. Exploratory MIC model metrics are reported separately in `docs/activity_baseline.md` and must not be presented as validated predictions.
# Resistant-isolate subcohort

`scripts/preprocessing/extract_resistant_subcohort.py` extracts rows only when
the DBAASP isolate/strain note itself explicitly contains MDR, XDR, PDR, CRAB,
carbapenem-resistant, multidrug-resistant, or extensively drug-resistant
language. Publication titles and peptide names are not used as labels. The
derived table is `data/processed/dbaasp_a_baumannii_resistant_isolate_records.csv`.

The downloaded snapshot yielded 612 explicitly labeled rows, including 406
exact-point MIC rows, 452 sequences with any sequence text, and 382 unique
standard-amino-acid sequences. The data still mixes units, strains, protocols,
and labels more than MDR alone; terms in source notes are not a harmonized
resistance phenotype. This is a traceable curation panel for follow-up review,
not a ready-to-train MDR classifier or independent test set. The CSV retains
the original notes and provenance fields so records can be traced back to
DBAASP and checked against their papers.

Recreate it with:

```powershell
.\.venv\Scripts\python.exe scripts\preprocessing\extract_resistant_subcohort.py
```
