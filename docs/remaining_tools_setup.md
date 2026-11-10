# Run the remaining peptide predictors and docking workflow

## Current setup status

| Component | Present in project | What is still needed |
|---|---|---|
| Activity data and exploratory models | Yes | No validated active/inactive labels or useful prospective predictor are available from these data. |
| Peptide controls and hypotheses | Yes | Generated five B2RP analogues; all still need prospective evidence. |
| FtsZ structures | Yes | 22HR apo and 22HQ GDP-bound structures are prepared for inspection. |
| ToxinPred3 | Online run complete | Ten sequences; Hybrid ET+MERCI, threshold 0.38; output imported. |
| Hemolysis predictor | Online runs complete | HemoPI2 HC50 regression and Hybrid1 classification on ten sequences; outputs imported. |
| Peptide–FtsZ docking | Web jobs complete | Ten correct panel sequences returned representative models and cluster tables from CABS-dock against apo 22HR chain A. No single binding-affinity score is reported; systematic pose review and benchmark validation remain. |
| Lookup UI | Yes | Shows candidate generation, imported toxicity predictions, experimental evidence, and direct docking-job links. |

Job IDs and raw predictor tables are recorded in
`docs/external_predictor_and_docking_runs.md`,
`results/candidate_prioritization/external_predictions/`, and
`results/docking/cabsdock/web_jobs.csv`.

The prepared sequence panel is
`results/candidate_prioritization/b2rp_candidate_panel.fasta` (one published
parent, three published controls, the top-ranked published lead as a positive
control, and five generated hypotheses). Its headers
preserve candidate IDs. The source table
`results/sequence_generation/b2rp_analogue_hypotheses.csv` records which
sequences are hypotheses and the assumed terminal modification. FASTA input
contains only the one-letter sequence; it does not encode the assumed
C-terminal amidation. The predictors and initial docking protocol therefore
evaluate unmodified sequence strings unless their separate inputs explicitly
support modifications.

## Safety prediction without installing packages

These websites receive the submitted peptide sequences, so use the batch file
only if sharing this candidate panel with the predictor providers is acceptable.
It contains public controls plus unpublished, unvalidated sequence hypotheses.

1. Open [ToxinPred3 Prediction](https://webs.iiitd.edu.in/raghava/toxinpred3/prediction.php).
2. Upload `b2rp_candidate_panel.fasta` or paste its contents. Select the hybrid
   model (ET+MERCI or DL+MERCI), record which model and threshold you choose,
   submit, then download/save the result table.
3. Open [HemoPI2 Prediction](https://webs.iiitd.edu.in/raghava/hemopi2/prediction.html)
   and submit the same FASTA once in **Regression** mode and once in
   **Classification** mode. Save the predicted HC50 table and the class/score
   table separately; record the selected model and threshold.
4. Save the downloaded result CSVs. Keep FASTA IDs unchanged so results can be
   joined to the panel. Use the project importer commands below; it verifies
   both ID and returned sequence before importing anything.

ToxinPred3 estimates general peptide toxicity; HemoPI2 estimates erythrocyte
hemolysis, including a modelled HC50. These are different endpoints and are
model predictions, not measurements or clinical safety judgments. HemoPI2's
paper defines its training classification threshold at HC50 ≤100 µM, but do
not treat that as a universal assay or therapeutic cutoff. The project does
not send these sequences to either website automatically.

**Hemolytik2 is a curated experimental-data resource, not the HemoPI2
predictor.** It can be used to search for matching experimental records, but a
missing match is not a negative result. Prefer the same sequence IDs and exact
assay conditions when adding source measurements to the UI.

## Local predictor installation (optional)

Keep these tools in separate environments rather than replacing the project's
pinned scikit-learn version. The official standalone packages list additional
dependencies; use their own instructions and model downloads.

```powershell
py -3 -m venv .tools\toxinpred3
.tools\toxinpred3\Scripts\python.exe -m pip install toxinpred3
.tools\toxinpred3\Scripts\toxinpred3.exe -h

py -3 -m venv .tools\hemopi2
.tools\hemopi2\Scripts\python.exe -m pip install hemopi2
.tools\hemopi2\Scripts\hemopi2_classification.exe -h
.tools\hemopi2\Scripts\hemopi2_regression.exe -h
```

If a package does not provide those Windows launchers or needs unavailable
model files, use the official websites or follow the respective standalone
installation guide. Do not install either package into the project's `.venv`;
ToxinPred3/HemoPI2 dependencies may conflict with the project's scikit-learn
version. The online route is simpler for this ten-sequence panel.

## Peptide–FtsZ docking installation

The prepared apo receptor is
`data/external/ftsZ/prepared/22HR_chainA_apo.pdb`. CABS-flex standalone 3's
official installer supports Linux, macOS, and Windows through WSL2. This
workspace currently returns `Wsl/AccessDenied`; an administrator may need to
enable/install WSL2, or use a Linux/macOS machine. Do not disable security
controls to work around that restriction. An Ubuntu app is installed in the
Windows app inventory, but it did not expose a targetable terminal from the
project session; if it opens normally when you launch it yourself, you can use
that Ubuntu environment for the commands below.

After an authorized Ubuntu/WSL2 installation, install from the official source:

```bash
sudo apt update
sudo apt install -y git curl tar
git clone https://github.com/LCBio/CABSflex_standalone.git
cd CABSflex_standalone
bash install.sh
micromamba activate cabs
CABSdock --help
```

Use the full CABSdock peptide–protein docking workflow documented at
[CABS-flex standalone 3](https://cabsflex.lcbio.pl/peptide-protein-docking.html).
Do not treat the maintainer's very short quick-start simulation as production
evidence. Start with the apo 22HR chain-A input; record the full command,
software/version, parameters, random seeds, input hashes, run logs, pose
clusters, and method-specific scores. Repeat across seeds and perform a
separate 22HQ receptor-state sensitivity run. Scores are docking ranks, not
binding free energies or experimental affinity. The project has no validated
AMP–FtsZ positive-control complex, so any candidate docking remains
exploratory.

The prepared batch runner covers all ten panel sequences, both receptor states,
and three random seeds. Review its 60-run plan and resource demand first; it is
not a quick demo command. After CABS-flex is installed, generate and launch it
from the project root with the two commands in `docs/docking_workflow.md`.

## Import predictor outputs into the interface

Run the commands from the project root in PowerShell, replacing the input file
and model/threshold with the values used for that run:

```powershell
.\.venv\Scripts\python.exe scripts\activity_prediction\import_external_safety_predictions.py --tool toxinpred3 --endpoint toxicity --input "C:\path\to\toxinpred3_results.csv" --model "Hybrid ET+MERCI" --threshold 0.38
.\.venv\Scripts\python.exe scripts\activity_prediction\import_external_safety_predictions.py --tool hemopi2 --endpoint hemolysis-hc50 --input "C:\path\to\hemopi2_regression.csv" --model "Random Forest Regressor"
.\.venv\Scripts\python.exe scripts\activity_prediction\import_external_safety_predictions.py --tool hemopi2 --endpoint hemolysis-classification --input "C:\path\to\hemopi2_classification.csv" --model "Hybrid1 ESM2-t6+MERCI" --threshold "<selected value>"
.\.venv\Scripts\python.exe ui\build_ui_data.py
```

The importer writes normalized provenance rows to
`results/candidate_prioritization/external_safety_predictions.csv`. It stops
without writing if an output row's ID is unknown, the sequence is missing or
does not exactly match the labeled FASTA, or no recognized result value exists.
The UI shows model-based toxicity/hemolysis results separately from
source-reported experimental measurements. FtsZ docking remains “Not run” until
a reproducible CABSdock result exists.

## Official tool references

- [ToxinPred3 prediction form](https://webs.iiitd.edu.in/raghava/toxinpred3/prediction.php)
- [ToxinPred3 standalone source and documentation](https://github.com/raghavagps/toxinpred3)
- [HemoPI2 prediction form](https://webs.iiitd.edu.in/raghava/hemopi2/prediction.html)
- [HemoPI2 standalone source and documentation](https://github.com/raghavagps/HemoPI2)
- [Hemolytik2 database and API guide](https://webs.iiitd.edu.in/raghava/hemolytik2/usr_guide.html)
- [CABS-flex/CABSdock installation guide](https://cabsflex.lcbio.pl/installation.html)
- [CABS-flex peptide–protein docking workflow](https://cabsflex.lcbio.pl/peptide-protein-docking.html)
