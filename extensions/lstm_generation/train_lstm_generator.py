"""Small-data character LSTM baseline for generating AMP sequence hypotheses."""
from __future__ import annotations

import csv
import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PANEL = ROOT / "data/processed/sequence_provenance_panel.csv"
OUT = HERE / "output"
AA = "ACDEFGHIKLMNPQRSTVWY"
PAD, BOS, EOS = "<PAD>", "<BOS>", "<EOS>"
SEED, EPOCHS, HIDDEN, EMBED, BATCH, N_GENERATE = 2026, 45, 96, 32, 64, 50


def read_sequences():
    sequences = set()
    with PANEL.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            seq = (row.get("sequence") or "").strip().upper()
            if 8 <= len(seq) <= 80 and set(seq) <= set(AA):
                sequences.add(seq)
    return sorted(sequences)


def main():
    try:
        import torch
        from torch import nn
        from torch.utils.data import DataLoader, Dataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required. Install extensions/lstm_generation/requirements.txt in .venv.") from exc

    random.seed(SEED); torch.manual_seed(SEED)
    torch.set_num_threads(max(1, min(4, torch.get_num_threads())))
    sequences = read_sequences()
    if len(sequences) < 100:
        raise SystemExit(f"Too few AMP sequences for the proposed prototype: {len(sequences)}")
    tokens = [PAD, BOS, EOS, *AA]
    to_id = {token: i for i, token in enumerate(tokens)}
    encoded = [[to_id[BOS], *[to_id[a] for a in s], to_id[EOS]] for s in sequences]
    shuffled = list(range(len(encoded))); random.Random(SEED).shuffle(shuffled)
    n_val = max(1, int(len(shuffled) * 0.2))
    val_ids, train_ids = shuffled[:n_val], shuffled[n_val:]

    class SeqDataset(Dataset):
        def __init__(self, ids): self.rows = [encoded[i] for i in ids]
        def __len__(self): return len(self.rows)
        def __getitem__(self, i): return self.rows[i]

    def collate(batch):
        width = max(map(len, batch)) - 1
        x = torch.zeros((len(batch), width), dtype=torch.long)
        y = torch.full((len(batch), width), -100, dtype=torch.long)
        for i, row in enumerate(batch):
            x[i, :len(row)-1] = torch.tensor(row[:-1])
            y[i, :len(row)-1] = torch.tensor(row[1:])
        return x, y

    class CharLSTM(nn.Module):
        def __init__(self):
            super().__init__()
            self.embedding = nn.Embedding(len(tokens), EMBED, padding_idx=0)
            self.lstm = nn.LSTM(EMBED, HIDDEN, batch_first=True)
            self.output = nn.Linear(HIDDEN, len(tokens))
        def forward(self, x):
            z, _ = self.lstm(self.embedding(x))
            return self.output(z)

    model = CharLSTM()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.003)
    loss_fn = nn.CrossEntropyLoss(ignore_index=-100)
    train_loader = DataLoader(SeqDataset(train_ids), batch_size=BATCH, shuffle=True, collate_fn=collate)
    val_loader = DataLoader(SeqDataset(val_ids), batch_size=BATCH, shuffle=False, collate_fn=collate)
    history = []
    for epoch in range(1, EPOCHS + 1):
        model.train(); train_losses = []
        for x, y in train_loader:
            optimizer.zero_grad(); loss = loss_fn(model(x).reshape(-1, len(tokens)), y.reshape(-1))
            loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
            train_losses.append(float(loss.detach()))
        model.eval(); val_losses = []
        with torch.no_grad():
            for x, y in val_loader:
                val_losses.append(float(loss_fn(model(x).reshape(-1, len(tokens)), y.reshape(-1))))
        history.append({"epoch": epoch, "train_loss": sum(train_losses)/len(train_losses),
                        "validation_loss": sum(val_losses)/len(val_losses)})

    known = set(sequences); generated = []
    model.eval()
    # Generate a larger temporary pool and deduplicate exact local matches.
    for attempt in range(N_GENERATE * 200):
        current = torch.tensor([[to_id[BOS]]], dtype=torch.long)
        result = []
        with torch.no_grad():
            for _ in range(80):
                logits = model(current)[0, -1]
                logits[to_id[PAD]] = -1e9; logits[to_id[BOS]] = -1e9
                probs = torch.softmax(logits / 0.85, dim=-1)
                nxt = int(torch.multinomial(probs, 1))
                if nxt == to_id[EOS]: break
                result.append(tokens[nxt]); current = torch.cat([current, torch.tensor([[nxt]])], dim=1)
        seq = "".join(result)
        if 10 <= len(seq) <= 50 and seq not in known and seq not in generated and set(seq) <= set(AA):
            generated.append(seq)
        if len(generated) >= N_GENERATE: break

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "lstm_hypotheses.csv").open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["candidate_id", "sequence", "length", "local_exact_match", "activity_status"])
        w.writeheader()
        for i, seq in enumerate(generated, 1):
            w.writerow({"candidate_id": f"LSTM_{i:03d}", "sequence": seq, "length": len(seq),
                        "local_exact_match": "no", "activity_status": "unmeasured; hypothesis only"})
    with (OUT / "lstm_hypotheses.fasta").open("w", encoding="ascii") as f:
        for i, seq in enumerate(generated, 1): f.write(f">LSTM_{i:03d}|unvalidated_hypothesis\n{seq}\n")
    torch.save(model.state_dict(), OUT / "amp_char_lstm_state.pt")
    report = {"model": "single-layer character LSTM", "seed": SEED, "n_unique_training_panel_sequences": len(sequences),
              "n_train_sequences": len(train_ids), "n_validation_sequences": len(val_ids),
              "epochs": EPOCHS, "embedding_dim": EMBED, "hidden_dim": HIDDEN,
              "final_train_loss": history[-1]["train_loss"], "final_validation_loss": history[-1]["validation_loss"],
              "training_history": history, "hypotheses_generated": len(generated),
              "generation_filter": "Standard amino acids; length 10-50; remove exact local sequence matches and duplicates; temperature 0.85.",
              "limitations": ["Small, heterogeneous AMP-only training data; random sequence split is not a novelty benchmark.",
                              "The generator optimizes next-residue likelihood, not activity, safety, novelty, or target binding.",
                              "No generated sequence has experimental activity or toxicity evidence."]}
    (OUT / "lstm_generation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("n_unique_training_panel_sequences", "n_train_sequences", "n_validation_sequences", "final_train_loss", "final_validation_loss", "hypotheses_generated", "limitations")}, indent=2))


if __name__ == "__main__": main()
