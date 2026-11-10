# Multi-bacterium activity evidence extension

This extension extracts MIC evidence for several bacterial species from the
same cached DBAASP peptide detail records used by the core project. It does
not merge assay units, convert MIC values, infer resistance from a species
name, or claim that all records form a harmonized multi-organism benchmark.

Run from the project root:

```powershell
.\.venv\Scripts\python.exe extensions\multi_target_screening\extract_multitarget_mic_panel.py
.\.venv\Scripts\python.exe extensions\multi_target_screening\multitarget_mic_baselines.py
.\.venv\Scripts\python.exe extensions\multi_target_screening\fetch_target_structures.py
.\.venv\Scripts\python.exe extensions\multi_target_screening\prepare_target_receptors.py
```

Outputs are isolated in `extensions/multi_target_screening/output/`. Each row
retains peptide ID, sequence, exact target string, reported concentration,
unit, assay context, and candidate publication provenance. The panel supports
an expanded descriptive analysis and target selection; target-specific model
training still needs assay and isolate harmonization.

The current evidence panel has 7,744 source MIC rows across *A. baumannii*,
*E. coli*, *P. aeruginosa*, *K. pneumoniae* and *S. aureus*. Target-specific
grouped baselines are in `output/multitarget_mic_cv_summary.json`; results are
weak, with Random Forest R-squared values near or below zero. The structure
panel adds experimental FtsZ structures for three species and precomputed
AlphaFold Database FtsZ models for *E. coli* and *S. aureus*. These are target
structures, not AlphaFold predictions of the peptide candidates. Multi-target
docking has not been executed.
