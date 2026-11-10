# Focused source review: resistant-isolate MIC records

## Review scope

I checked the DBAASP candidate PMID mapping against the article array in the
cached peptide detail records and reviewed the corresponding PubMed records.
The following three primary studies are a focused manual audit, not yet a
complete review of all 109 candidate PMIDs seen in the resistant-note panel.
Record counts below are from the local downloaded DBAASP MIC snapshot and count
measurements, not independent experiments.

## Reviewed studies

| PMID | Evidence in the paper | DBAASP rows mapped to paper | Interpretation |
|---|---|---:|---|
| 19793185 | B2RP and analogues inhibited nosocomial MDR *A. baumannii* isolates; the abstract reports B2RP MIC 3–6 µM and [Lys4]B2RP 1.5–3 µM. | 28 rows, 8 peptide records, all exact µM values | Direct MDR evidence. The paper title/abstract identify MDR, but most local isolate-note values are blank or only say “Clinical isolate.” Thus the explicit-note cohort alone has false negatives for paper-supported resistant evidence. |
| 24071034 | LL-37, KS-30, KR-20, and KR-12 were evaluated against MDR clinical isolates. The abstract reports MIC ranges, including 16–32 µg/mL for LL-37 and 8–16 µg/mL for KS-30. | 8 rows, 4 peptide records, µg/mL; four rows explicitly name MDR clinical isolates, four lack that note | Direct MDR evidence. Keep ranges as ranges and do not average them into point MICs. |
| 40056800 | A 2025 study synthesized cathelicidin-derived dodecapeptides for CRAB. The abstract reports low MICs of 2.5–20 µg/mL for its leading candidates and describes membrane binding/penetration; d24 is hydrocarbon-stapled. | 65 rows, 31 peptide records; 53 exact point values and remaining values include non-point formats, all µg/mL | Direct CRAB evidence. Preserve modified peptide chemistry: a linear one-letter sequence is not a complete representation of the stapled d24 molecule. |

The article identities and claims above are supported by the respective
[PubMed records for B2RP](https://pubmed.ncbi.nlm.nih.gov/19793185/),
[LL-37 fragments](https://pubmed.ncbi.nlm.nih.gov/24071034/), and
[CRAB dodecapeptides](https://pubmed.ncbi.nlm.nih.gov/40056800/). In the cached
DBAASP detail JSON, each activity's numeric `reference` matches the
corresponding 1-based article position for these records. That verifies the
index interpretation for these examples only; it does not automatically
validate the mapping for every peptide response.

## Consequences for analysis

1. Keep two provenance strata: explicit resistance wording in the activity's
   strain/isolate field, and resistance evidence established in a reviewed
   primary paper. Do not label every record from an article as resistant if
   the article includes susceptible or unspecified isolates; inspect the paper's
   methods and its mapping of peptide/strain/assay records first.
2. Treat the 612-row explicit-note table as a high-specificity review queue,
   not a complete resistant activity set. Paper-linked records such as B2RP
   show that note-only selection misses relevant studies.
3. Keep µM and mass concentration units separate unless peptide molecular form
   and molecular mass are confirmed. Keep ranges, bounds, and point values
   separate.
4. Do not regard 31 records from one paper as 31 independent publications or
   as independent external validation. Group validation by source publication
   and sequence family.
5. For modified peptides (including stapled analogues), retain the exact
   chemical modification record. Exclude them from an unmodified linear
   sequence model unless a chemistry-aware encoding is added.

## Reproduction and next curation step

The 101-row linked review panel is
`data/processed/dbaasp_a_baumannii_reviewed_paper_records.csv`; its extraction
script is `scripts/preprocessing/extract_paper_review_panel.py`. It adds the
reviewed PMID and a paper-level resistance-evidence note, while explicitly
marking each row-level isolate assignment as not yet reviewed. Recreate it with:

```powershell
.\.venv\Scripts\python.exe scripts\preprocessing\extract_paper_review_panel.py
```

The row counts can also be regenerated from
`data/processed/dbaasp_a_baumannii_mic_parsed.csv` by matching each
`candidate_pubmed_ids_json` field to the PMIDs above. Next review should verify
paper-specific strain susceptibility, peptide identities/modifications,
assay conditions, and whether each MIC is associated with resistant isolates.
Do that at the source-table level before adding labels or making a resistant
benchmark. No additional activity labels are created by this document.
