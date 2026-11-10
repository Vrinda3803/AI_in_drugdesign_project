# Evidence-based pilot candidate prioritization

## Scope

The strongest current integrated activity/safety evidence is a source-matched
set of B2RP analogues reported in one paper. This is a **retrospective lead
benchmark**, not a novel generated sequence and not an independent validation
set. The five named NM isolates in DBAASP are linked to PubMed PMID 19793185;
the paper's abstract describes activity against nosocomial MDR *A. baumannii*
isolates and reports MICs and hemolytic activity. See the
[primary paper record](https://pubmed.ncbi.nlm.nih.gov/19793185/).

The pipeline matched MIC rows and 50% human erythrocyte hemolysis observations
by DBAASP peptide ID and candidate PMID. It considered MIC values against the
five named NM isolates and used the worst (highest) MIC among those isolates.
For hemolysis, it used the reported 50% hemolysis concentration. Ratios are
shown only where both are in µM, and a censored H50 such as `>200` yields a
lower bound. These ratios are descriptive concentration margins from one
paper, not validated therapeutic indices; assays can differ in protocol.

## Results

| Within-paper rank | DBAASP ID | Variant | MIC across five NM isolates (µM) | Highest MIC | Human erythrocyte 50% hemolysis (µM) | Apparent H50 / highest MIC |
|---:|---:|---|---|---:|---:|---:|
| 1 | 705 | [D4K][K16A][L18K] B2RP | 3, 3, 3, 3, 3 | 3 | >200 | >66.67 |
| 2 | 703 | [D4K][L18K] B2RP | 6, 6, 3, 6, 6 | 6 | >200 | >33.33 |
| 3 | 699 | [D4K] B2RP | 3, 3, 3, 1.5, 1.5 | 3 | 95 | 31.67 |
| 4 | 14 | B2RP parent | 6, 6, 3, 6, 6 | 6 | 90 | 15.00 |

The first-ranked variant is a reasonable **published positive-control lead to
reproduce**, given its consistent reported MIC across the five named isolates
and the strongest conservative hemolysis-to-MIC margin in this small panel.
The D4K variant is a useful potency comparator; the parent peptide is a
reference. The ranking only sorts these four variants under this transparent
rule and does not indicate that the first-ranked sequence is novel or
experimentally confirmed in the present project.

## Reproduction and limitations

Run:

```powershell
.\.venv\Scripts\python.exe scripts\activity_prediction\build_b2rp_lead_benchmark.py
```

Outputs are `results/candidate_prioritization/b2rp_published_lead_benchmark.csv`
and `.json`. The calculation uses source-level reported measurements already
present in DBAASP; no classifier or activity predictor is involved. Before
using the result for prospective design, inspect the full methods and tables,
confirm strain identities/resistance phenotypes, assay conditions, and the
chemical form of each analogue. These four sequences already exist in the
source publication and cannot be presented as newly designed AMPs.
