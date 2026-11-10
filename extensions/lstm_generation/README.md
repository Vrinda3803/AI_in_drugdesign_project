# Character-level LSTM generation baseline

This optional extension implements the LSTM baseline proposed at mid-review.
It trains on the unique standard-amino-acid sequences in the local AMP panel
and samples new sequence hypotheses. It does not train an activity predictor,
conditional generator, or reinforcement-learning agent. Generated sequences
are not evidence of activity, safety, novelty, or target engagement.

Install the extension-only dependency and run from the project root:

```powershell
.\.venv\Scripts\python.exe -m pip install -r extensions\lstm_generation\requirements.txt
.\.venv\Scripts\python.exe extensions\lstm_generation\train_lstm_generator.py
.\.venv\Scripts\python.exe extensions\lstm_generation\audit_generated_candidates.py
```

The run records the training/validation split, seed, hyperparameters, loss,
and generation filters under `extensions/lstm_generation/output/`. Exact local
matches are removed. This is a small-data proof of implementation; report
candidate sequences as hypotheses requiring synthesis and laboratory tests.
The current run generated 50 locally exact-match-novel sequences from 960
unique training peptides. Validation loss remained above training loss, which
is consistent with overfitting; the model has no activity objective and these
samples have not been experimentally tested. `output/lstm_candidate_audit.csv`
adds sequence-derived descriptors and nearest-sequence comparison against the
same 960-sequence local training panel. This audit is not an activity or safety
predictor and is not an external novelty search. Toxicity/hemolysis predictions
for these 50 candidates are not available: submissions to the official public
forms did not return results (browser automation stalled; direct requests
returned HTTP 400). The attempt and its status are recorded in
`output/external_prediction_attempts.json`.
