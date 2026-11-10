# Planned workflow completion matrix

This matrix compares the project plan with the current implementation and
states what the results support. “Done” means the code or source-data output
exists; it does not imply biological or clinical validation.

## End-to-end workflow

| Planned stage | Status | Current output and remaining work |
|---|---|---|
| Define MDR *A. baumannii* and FtsZ question | Done | Project scope, FtsZ structures, and intracellular-target limitations are documented. |
| Collect public AMP data | Done | AbAMPdb export and 1,126 DBAASP detail records are cached with source provenance. |
| Clean and QC sequences/data | Done with limits | Sequence QC, descriptors, and duplicate/provenance panels are available. Nonstandard chemistry is preserved or excluded where a model requires standard residues. |
| Build *A. baumannii* activity/MIC dataset | Done | 2,134 activity observations and 1,834 parsed MIC records; units and censored/range values are retained. |
| Annotate MDR/XDR evidence | Partial | A 612-row note-based queue and 101-row paper-linked panel exist. Isolate/resistance mappings still need expert source-paper verification. |
| Deduplicate and prevent train/test leakage | Partial | Exact sequences are collapsed for modeling; grouped validation joins ≥80% edit-identity components and candidate publication IDs. Some article mappings are unresolved, and there is no independent external test set. |
| Engineer sequence/physicochemical features | Done for standard amino acids | Length, amino-acid composition, charge, molecular weight, pI, hydrophobicity, and hydrophobic moment are calculated. Learned embeddings were not added. |
| Train and compare activity models | Done as feasibility analysis | Regression: median, Ridge, RBF-SVR, Random Forest. Relative-potency classification: majority dummy, logistic regression, RBF-SVM, Random Forest. Standard requested metrics are saved. No validated active/inactive labels exist; results are not reliable candidate MIC predictions. |
| Generate novel AMP sequences | Partial by design | Five source-informed B2RP substitutions are generated. They are hand-specified hypotheses, not LSTM/AI-generated sequences; local exact-match checks do not prove global novelty. |
| Filter novelty and diversity | Partial | Exact local matches/nearest-sequence distances are reported, and preliminary exact-phrase web searches found no indexed matches for the five hypotheses. No comprehensive literature/patent search or broad library diversity optimization has been completed. |
| Screen toxicity/hemolysis | Done computationally for the panel | ToxinPred3 and HemoPI2 predictions exist for ten sequences; source assay records are also displayed. Predictions do not establish safety. Hemolytik2 is used as a source-data resource, not as a prediction model. |
| Evaluate peptides against FtsZ | Initial screen done | Ten CABS-dock jobs against apo 22HR chain A returned representative poses and clusters; top-cluster summaries are recorded for all ten. Model 1 was downloaded and sequence-checked for each job, and all ten have a first-pass contact comparison. The ten-model ensemble for the published lead control also has contact summaries. Full visual review, repeated runs, and benchmark validation remain. There is no binding-affinity score. |
| Multi-objective prioritization | Partial | Matched published MIC/H50 controls are Pareto-tiered; generated hypotheses remain unscored without measured activity. The UI presents separate evidence rather than a fabricated combined score. |
| Produce final experimental candidate set | Candidate shortlist ready for review | Five untested hypotheses and published controls are in the panel. Synthesis chemistry and lab selection still require the project team. |
| Test generated peptides experimentally | Not done | Requires peptide synthesis and validated microbiology/hemolysis/cytotoxicity assays. No experimental data are present for the generated hypotheses. See `docs/experimental_validation_handoff.md`. |

## Planned methods and software

| Plan item | Status |
|---|---|
| Logistic regression | Implemented for relative-potency classification. |
| Random Forest | Implemented for MIC regression, classification, and H50 feasibility analysis. |
| SVM/SVR | Implemented for relative-potency classification and MIC regression. |
| Deep-learning activity model | Not implemented. Grouped classical-model performance is modest, cohorts are about 255 sequences per MIC unit, and this dataset does not yet support a meaningful deep-learning claim. |
| Classification metrics | Accuracy, precision, recall, F1, ROC-AUC, and PR-AUC are saved by fold. Labels are relative MIC bins, not active/inactive. |
| Regression metrics | RMSE, MAE, and R² are saved by fold on log2 MIC. |
| LSTM sequence generator | Not implemented. Current candidate generator is a reproducible five-sequence SAR hypothesis panel; adding an LSTM on this small, heterogeneous dataset would create unvalidated sequences without a sound activity objective. |
| Conditional generation / deep RL | Not implemented; no validated objective or prospective activity signal is available to optimize. |
| PyMOL / structure preparation | Done for chain-A structure inspection and alignment. Full receptor preparation and docking validation remain. |
| SwissDock | Not used for AMP docking because its documented purpose is protein–small-molecule docking. Flexible peptide–protein jobs used CABS-dock. |
| Python / scikit-learn / NumPy | Used by the reproducible pipeline. Data are retained as CSV/JSON; SQLite, Pandas, Colab, PyTorch, and TensorFlow are not required by the current implementation. |
| User interface | Done as an offline-friendly browser workbench with candidate generation, source evidence, toxicity predictions, and direct pose/cluster links. |

## Research completion boundary

The reproducible software and the first computational screen are working. The
planned project is not complete as a biological discovery study until the
resistance mappings and poses are reviewed and the generated peptides are
synthesized and tested. Activity/deep-learning claims should wait for an
independent, publication-aware validation set or prospective MIC measurements.
