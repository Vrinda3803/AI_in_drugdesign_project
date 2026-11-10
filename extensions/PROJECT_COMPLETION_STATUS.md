# Project completion status against the mid-review plan

This status sheet distinguishes the original plan from analyses completed in
the main project and its isolated extensions. “Implemented” means an analysis
or artifact exists; it does not mean biological validation.

## Scope and objective

The project studies whether public AMP evidence can support reproducible
computational prioritization for MDR *Acinetobacter baumannii*. FtsZ is a
structural target hypothesis. The current evidence does not show that FtsZ
binding explains antibacterial activity or that peptides enter bacterial
cells.

## Plan compared with delivered work

| Mid-review stage | Delivered work | Status and safe interpretation |
|---|---|---|
| Define MDR pathogen and FtsZ target | Project question and FtsZ rationale documented; 612-row resistant-note queue and 101-row paper-linked panel | Partial: isolate/resistance mappings need source-paper verification. |
| Collect AMP data | AbAMPdb export and 1,126 DBAASP detail records cached with provenance | Complete for the saved snapshot. AbAMPdb MIC units are absent, so those values are not pooled with DBAASP. |
| Build multi-bacterium evidence | Separate extension extracts 7,744 MIC rows for five species from the A. baumannii-selected DBAASP detail records | Partial: useful cross-target evidence, but not a comprehensive or assay-harmonized multi-species database. |
| Preprocess sequences and MICs | Source QC, exact deduplication, unit separation, and range/censor preservation | Complete for the documented source snapshot; uncertain values remain outside exact-point models. |
| Engineer features | Length, amino-acid composition, MW, charge, pI, hydrophobicity, hydrophobic moment; extension adds instability index | Implemented; features are sequence-derived estimates and omit many peptide modifications. |
| Train activity prediction | Unit-separated MIC regression and relative-potency classification; multi-target extension adds five-target grouped regression | Feasibility analysis only. Core RF strict-group R² is 0.16 (µM) and -0.44 (µg/mL). Multi-target RF R² is generally near or below zero. No reliable MDR-specific candidate MIC predictor is established. |
| AMP versus non-AMP classification | Separate binary classifier trained on 974 local AMP sequences and 974 length-matched UniProt-derived fragments | Prototype only. Negatives are weak labels; unusually high grouped AUC (about 0.986–0.991) may reflect source/composition differences and is not validated detector performance. |
| Generate new peptides | Five auditable B2RP SAR hypotheses plus a character-LSTM baseline that generated 50 exact-local-match-novel sequences from 960 unique training peptides; a local descriptor and nearest-sequence audit was added for all 50 | Implemented as sequence hypotheses. LSTM validation loss remained higher than training loss; neither set has predicted or measured activity evidence. The new audit is not an activity predictor or global novelty search. No conditional generation or deep RL. |
| Safety screen | ToxinPred3 and HemoPI2 predictions for the original ten-sequence panel; source-reported H50 records kept separately | Complete for that panel only. An authorized attempt to submit the 50 LSTM sequences to both public predictors did not return a prediction: browser submission stalled and direct requests returned HTTP 400. Keep LSTM safety as not screened. Existing predictions are not safety clearance. |
| AlphaFold structure prediction | Precomputed AlphaFold DB FtsZ models retrieved for *E. coli* and *S. aureus*; experimental FtsZ structures are available for *A. baumannii*, *E. coli*, and *S. aureus* | Partial: these are target structures, not de novo predictions of peptide candidates. No peptide AlphaFold prediction was used. AF3 Server outputs are not used in docking; its terms prohibit using those outputs in docking/screening. |
| Multiple targets and docking | Five-species MIC evidence panel; protein-only FtsZ receptor inputs prepared for three species; ten CABS-dock jobs already completed against *A. baumannii* apo 22HR | Partial: no multi-species docking jobs were run. Existing 22HR poses are exploratory and have no affinity score. Standalone CABS-flex was unavailable in the Windows environment described in the project notes. |
| Final prioritization | Transparent same-paper B2RP lead benchmark and measured-evidence Pareto tiering; generated sequences remain unscored | Candidate set ready for review, not an evidence-backed ranking of untested hypotheses. |
| Experimental validation | No synthesis, MIC, hemolysis, or target-engagement experiment in the repository | Not done; required before claiming an active or safe candidate. |

## Current candidate groups

1. **Published positive control:** [D4K][K16A][L18K] B2RP, reported in the source study; it is neither newly discovered nor experimentally reproduced by this project.
2. **Five source-informed hypotheses:** Arg substitutions at B2RP positions 4 and/or 18; no measured MIC or toxicity result.
3. **Fifty LSTM samples:** zero exact matches against the local 960-sequence training panel; validation loss indicates overfitting, and no activity or safety prediction was used to rank them. Sequence descriptors and nearest local sequences are in `extensions/lstm_generation/output/lstm_candidate_audit.csv`.

Do not combine docking geometry, source MIC evidence, and weak model scores into one candidate score. Keep measured controls, hypotheses, and generated samples clearly distinguished.

## Reproduction entry points

- Main pipeline: `scripts/run_reproducible_pipeline.ps1`
- AMP/non-AMP prototype: `extensions/amp_binary_classifier/README.md`
- LSTM baseline: `extensions/lstm_generation/README.md`
- Multi-target MIC and structure panel: `extensions/multi_target_screening/README.md`

## Completion boundary

The computational project now contains a binary-classification prototype,
multi-target activity evidence, a baseline LSTM generator, AlphaFold Database
target structures, and an expanded receptor panel. It remains a computational
feasibility project. Independent validation, comprehensive multi-target
docking, peptide-specific structure assessment, expert resistance mapping,
and wet-lab testing are still needed for a biological discovery claim.
