# AMP versus non-AMP classifier extension

This is an isolated extension for Problem 2. It does not change the main
pipeline. It adds the missing instability-index descriptor and a reproducible
binary-classification workflow, but it deliberately requires a curated
negative-sequence CSV before fitting a model.

## Data needed

The first candidate set can be recreated from UniProt using:

```powershell
.\.venv\Scripts\python.exe extensions\amp_binary_classifier\fetch_uniprot_putative_negatives.py
```

This fetches reviewed *E. coli* proteins, excludes obvious AMP/antibiotic/
defense/toxin annotations, then samples random contiguous fragments with the
same length distribution as local positive sequences. Each fragment retains
its UniProt accession and parent coordinates. These are **putative negatives**:
the parent annotation does not prove that an isolated fragment has no
antimicrobial activity. Do not call these experimentally validated negatives.
UniProt describes reviewed entries as manually annotated, while annotation
absence does not prove biological absence of a function
([UniProtKB overview](https://www.uniprot.org/help/uniprotkb)).

Alternatively, provide `data/non_amp_sequences.csv` with columns:

```csv
sequence,label,source,source_id
ACDEFGHIKLMNPQRSTVWY,0,UniProt,ACCESSION
```

The builder rejects nonstandard residues and sequences outside 5–100 amino
acids. It excludes exact positive/negative sequence collisions and records
them in the QC report. The classifier keeps fragments from the same UniProt
parent in one split group, in addition to grouping close sequence neighbors.

Positive sequences are sourced from the local AbAMPdb export and DBAASP cache.
Because the dataset is assembled from AMP databases rather than an independent
curated benchmark, the output is exploratory and should not be presented as a
validated AMP detector.

## Run

Install the extension-only feature dependency if needed, then create the
negative CSV and run the two scripts from the project root:

```powershell
.\.venv\Scripts\python.exe -m pip install -r extensions\amp_binary_classifier\requirements.txt
.\.venv\Scripts\python.exe extensions\amp_binary_classifier\build_dataset.py
.\.venv\Scripts\python.exe extensions\amp_binary_classifier\train_classifier.py
```

Outputs go only to `extensions/amp_binary_classifier/output/`: the assembled
feature table and QC report, grouped cross-validation metrics, final fitted
pipeline, and model metadata. Evaluation groups sequences into connected
components at 80% normalized Levenshtein identity before splitting to reduce
close-homologue leakage; fragments from one UniProt parent accession are also
kept in one group. Results still require an independent external test set.
Because AMP-positive records and negative fragments come from different
sources, unusually high scores may reflect source or composition artifacts.
The main project's `requirements.txt` is unchanged.

## Initial run

The current snapshot has 974 AMP-positive sequences and 974 putative-negative
fragments from 365 eligible parent proteins. Five grouped holdouts produced
ROC-AUC values around 0.986–0.991 across the three baseline models. Treat this
as a **dataset-separability result only**, not evidence that a model can detect
true AMPs: the labels are weak, the sources differ, and high scores may
reflect source/composition signatures. An independently curated benchmark and
experimental negatives are needed before making biological performance claims.

## Added descriptor

The instability index uses the Guruprasad dipeptide instability weights.
Values below 40 are conventionally described as stable and values at/above 40
as unstable in that heuristic. This is a sequence-based estimate, not a
measured stability assay. Scores are not calculated for sequences shorter than
two residues.
