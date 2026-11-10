# Source-informed B2RP analogue hypotheses

## Purpose and design rule

The current MIC regressor did not generalize reliably under sequence- and
publication-aware evaluation, so it is not used to generate or rank sequences.
Instead, this script creates a compact, transparent hypothesis panel around
the B2RP parent using the source paper's experimentally varied positions:
Asp4 and Leu18. The source paper reports D4K as an activity-improving change;
L18K reduces potency but hemolysis; and the combined D4K/L18K variant retains
activity with lower hemolytic activity. These known variants are controls.

Five untested Arg substitutions explore conservative basic-residue alternatives
at these positions. This is a hand-specified SAR hypothesis, not an AI model
prediction. All hypotheses are unvalidated. K16A was not used as a design
direction because the source paper reports increased hemolytic activity for
the relevant K16A/K16L changes. See the [primary B2RP paper](https://pubmed.ncbi.nlm.nih.gov/19793185/).

## Generated panel

| Variant label | Sequence | Status in the local dataset |
|---|---|---|
| D4R | `GIWRTIKSMGKVFAGKILQNL` | Not found in local exact-sequence panel |
| L18R | `GIWDTIKSMGKVFAGKIRQNL` | Not found in local exact-sequence panel |
| D4K/L18R | `GIWKTIKSMGKVFAGKIRQNL` | Not found in local exact-sequence panel |
| D4R/L18K | `GIWRTIKSMGKVFAGKIKQNL` | Not found in local exact-sequence panel |
| D4R/L18R | `GIWRTIKSMGKVFAGKIRQNL` | Not found in local exact-sequence panel |

The local panel contains 981 exact standard-amino-acid sequences from the
downloaded AbAMPdb/DBAASP data. “Not found” is only a local exact-match result;
it is not a literature or patent novelty search. The output includes a nearest
local sequence and simple edit identity for audit, not as a biological
similarity score.

## Preliminary public exact-sequence search

On 2026-11-10, exact-phrase web searches were run for each of the five
hypothesis sequences above. The search returned no indexed exact-sequence
matches. This is a quick preliminary screen only; search indexing is incomplete
and this does not establish global novelty, patent freedom, or absence from
non-indexed databases. A formal literature and patent search remains necessary
before making a novelty claim.

## Chemical assumptions and limits

The source reports a C-terminally amidated B2RP parent. Candidate sequences
inherit that intended terminal modification as an assumption; confirm it before
synthesis. No predicted MIC, hemolysis, cytotoxicity, bacterial uptake,
FtsZ-binding, or efficacy values are produced. Do not describe these sequences
as active, safe, validated, globally novel, or AI-predicted. The next gate is
source/method review followed by synthesis and matched MIC/hemolysis testing.

## Reproduce

```powershell
.\.venv\Scripts\python.exe scripts\sequence_generation\generate_b2rp_hypotheses.py
```

Outputs:

- `results/sequence_generation/b2rp_analogue_hypotheses.csv`
- `results/sequence_generation/b2rp_analogue_generation_report.json`
