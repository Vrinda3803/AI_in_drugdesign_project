# Experimental and prospective validation handoff

## What the computation establishes

- Five B2RP source-informed sequences are hypotheses, not tested antimicrobial
  agents. Each currently inherits an intended C-terminal amide as an
  assumption; confirm exact terminal chemistry and purity before ordering or
  interpreting results.
- The MIC models have only modest grouped cross-validation performance and
  do not support reliable prospective predictions for these five sequences.
  A candidate with no source MIC record must be reported as unmeasured, not
  inactive.
- ToxinPred3 and HemoPI2 outputs are sequence-only model predictions. They
  cannot establish low toxicity or safety.
- CABS-dock returned pose clusters for ten sequences against one apo FtsZ
  monomer state. This is exploratory docking, not experimental target
  engagement or proof that FtsZ mediates antibacterial activity.

## Recommended next work with the wet-lab team

1. Confirm the exact sequence and terminal modifications for each hypothesis
   and reference control against the source paper and intended synthesis
   specification. Record purity and lot information for every synthesized
   peptide.
2. Pre-register a blinded, replicated comparison panel: the five hypotheses,
   B2RP parent, relevant published substitution controls, the published
   positive-control lead, and an assay-appropriate negative control.
3. Measure antibacterial activity against authenticated *A. baumannii*
   reference and clinically relevant resistant isolates using the lab's
   validated antimicrobial susceptibility workflow. Report raw replicate
   values, MIC units, medium, inoculum, isolate identity/resistance profile,
   and all controls. Do not pool results across different assay conditions.
4. Measure hemolysis and a relevant mammalian-cell cytotoxicity endpoint for
   the same peptide lots and concentration range. Report raw dose-response
   data, replicate variability, cell type, and assay controls. Treat model
   predictions as a comparison only.
5. Repeat unexpected or borderline results, then report activity and toxicity
   side by side. Advance a candidate only using pre-agreed criteria and
   independent review of the raw data.

## Recommended next work with the computational team

- Review the returned CABS-dock poses and contact maps for receptor artifacts,
  peptide placement, and consistency across clusters. Do not rank candidates
  by cluster density as if it were a binding free energy.
- Validate the docking setup against a known peptide-protein complex, and
  repeat candidate docking across recorded random seeds and receptor states
  before using pose recurrence for prioritization.
- Rebuild activity models only after verifying paper, isolate, assay, and
  sequence mappings. Use publication-aware external validation and compare
  against simple baselines. The current results are not suitable for candidate
  MIC prediction.
- When new MIC and toxicity measurements arrive, keep them as a prospectively
  collected validation set. Do not use them to tune the model before first
  reporting its blinded predictions.

## Data package to return to the project

For each peptide, return the verified sequence and chemistry, supplier/lot and
purity, isolate identifiers, assay protocol reference, raw replicate data,
controls, units, dates, and exclusions. The analysis script should ingest raw
values and produce a traceable summary without converting missing results to
inactive or safe labels.
