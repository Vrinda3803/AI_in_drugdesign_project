# Raw dataset acquisition and provenance

Retrieved 2026-11-09 from the primary/public sources listed in `dataset_manifest.json`.

- **AbAMPdb**: author repository exports. `ampdb-Final.csv` has 250 rows; 248 sequences use only the 20 standard amino acids and there are 237 unique exact sequences. The exported `MIC value` field does not specify units, so it is retained only as source data and is not mixed with DBAASP values. The exported `ampdb-Final3.csv` is headerless and `ampdb-Final2.csv` is a one-row fragment; these are not treated as separate validated data tables. Cite Anwer et al. (2024), DOI: 10.1093/database/baae096. Repository license is MIT; the database page describes access as free.
- **DBAASP**: full 25,542-peptide public API index plus all 1,126 detail responses returned by the *A. baumannii* search. The detail responses contain 2,134 target activity records, including 1,834 MIC measurement rows associated with 1,066 unique DBAASP peptide IDs. Those counts closely match the earlier project notes; the current API returned five more MIC rows than the cited 1,829. Units, concentration strings, strain/isolate notes, assay media, and references are kept as reported.

`data/processed/` contains source-level tables and a QC report. No MIC ranges/censoring were parsed into point values, no activity threshold was applied, and no missing assay result was labeled inactive. SHA-256 values for downloaded source files are recorded in the manifest. The prompt-referenced cleaned MIC CSV was not in the original workspace; this dataset was freshly retrieved from DBAASP.
