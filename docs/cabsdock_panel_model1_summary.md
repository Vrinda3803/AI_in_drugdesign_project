# Representative CABS-dock pose comparison across the panel

The representative `model_1` for each of the ten jobs was downloaded, checked
for exact chain-A target and chain-D peptide sequence, and analyzed for simple
heavy-atom contacts. The model files are in
`results/docking/cabsdock/web_models/<job_id>/`. The full table is
`results/docking/cabsdock/web_panel_model1_contact_summary.csv`.

| Candidate | Minimum inter-chain distance (Å) | Atom pairs ≤4 Å | Atom pairs ≤5 Å | Target residues ≤5 Å | Peptide residues ≤5 Å | CABS top-cluster density |
|---|---:|---:|---:|---:|---:|---:|
| D4R hypothesis | 2.735 | 199 | 687 | 42 | 20 | 41.009 |
| D4K control | 2.841 | 104 | 400 | 30 | 15 | 27.0812 |
| L18K control | 2.741 | 174 | 574 | 37 | 18 | 42.569 |
| D4K/L18K control | 2.701 | 156 | 577 | 37 | 15 | 42.0051 |
| L18R hypothesis | 2.935 | 84 | 302 | 21 | 12 | 35.8457 |
| D4K/L18R hypothesis | 2.732 | 218 | 748 | 45 | 19 | 32.9871 |
| D4R/L18K hypothesis | 2.679 | 143 | 521 | 33 | 19 | 34.4454 |
| D4R/L18R hypothesis | 2.707 | 117 | 459 | 24 | 17 | 38.9495 |
| B2RP parent | 2.741 | 129 | 477 | 27 | 17 | 31.9966 |
| Published positive-control lead | 2.653 | 123 | 499 | 39 | 17 | 26.8454 |

This is a first-pass comparison of single-run representative poses from the
same server protocol and the same apo monomer. It is not an affinity ranking.
Contact counts depend on pose geometry; cluster density is a separate
method-specific clustering statistic. Larger values for either field do not
establish stronger binding. The five generated hypotheses still have no
measured MIC or experimental safety results. Multi-seed docking, an appropriate
benchmark, visual inspection of interfaces, and experimental validation remain
necessary before using these outputs to prioritize biological claims.

Reproduce the full-panel summary from the project root:

```powershell
.\.venv\Scripts\python.exe scripts\docking\summarize_cabsdock_panel_model1_contacts.py
```
