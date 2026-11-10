"""Train a cautious, grouped baseline AMP classifier on extension data."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, average_precision_score, balanced_accuracy_score,
                             f1_score, precision_score, recall_score, roc_auc_score)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

HERE = Path(__file__).resolve().parent
DATA = HERE / "output" / "amp_binary_dataset.csv"
OUT = HERE / "output"
FEATURES = ["length", "molecular_weight_g_mol", "net_charge_pH7", "isoelectric_point",
            "hydrophobicity_gravy", "hydrophobic_moment", "instability_index"]
FEATURES += [f"AAC_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"]


def identity_groups(sequences: list[str], sources: list[str], source_ids: list[str], threshold: float = 0.8) -> list[int]:
    parent = list(range(len(sequences)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    def union(i, j):
        a, b = find(i), find(j)
        if a != b: parent[max(a, b)] = min(a, b)
    def edit(a, b):
        prev = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            cur = [i]
            for j, cb in enumerate(b, 1):
                cur.append(min(cur[-1] + 1, prev[j] + 1, prev[j-1] + (ca != cb)))
            prev = cur
        return prev[-1]
    for i, a in enumerate(sequences):
        for j in range(i + 1, len(sequences)):
            b = sequences[j]
            n = max(len(a), len(b))
            if n and abs(len(a)-len(b))/n <= 1-threshold and (n-edit(a,b))/n >= threshold:
                union(i,j)
    first_by_parent = {}
    for i, (source, accession) in enumerate(zip(sources, source_ids)):
        if source.startswith("UniProtKB/") and accession:
            if accession in first_by_parent:
                union(i, first_by_parent[accession])
            else:
                first_by_parent[accession] = i
    roots = [find(i) for i in range(len(sequences))]
    return [roots.index(root) for root in roots]


def model_set():
    return {
        "logistic_regression": make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=3000, random_state=42)),
        "rbf_svm": make_pipeline(StandardScaler(), SVC(C=1.0, kernel="rbf", class_weight="balanced", probability=True, random_state=42)),
        "random_forest": RandomForestClassifier(n_estimators=400, min_samples_leaf=3, class_weight="balanced_subsample", random_state=42, n_jobs=-1),
    }


def main():
    if not DATA.exists():
        raise SystemExit("Run build_dataset.py first; a curated negative CSV is required.")
    with DATA.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    y = np.array([int(r["label"]) for r in rows])
    if len(set(y)) != 2 or min(np.bincount(y)) < 5:
        raise SystemExit("Need both classes and at least 5 examples per class; inspect dataset_qc.json.")
    X = np.array([[float(r[k]) for k in FEATURES] for r in rows])
    seqs = [r["sequence"] for r in rows]
    groups = np.array(identity_groups(seqs, [r["source"] for r in rows], [r["source_id"] for r in rows]))
    splitter = GroupShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
    scores, final_models = [], {}
    for name, model in model_set().items():
        fold_rows = []
        for fold, (tr, te) in enumerate(splitter.split(X, y, groups), 1):
            if len(set(y[tr])) < 2 or len(set(y[te])) < 2:
                continue
            model.fit(X[tr], y[tr])
            prob = model.predict_proba(X[te])[:, 1]
            pred = (prob >= 0.5).astype(int)
            fold_rows.append({"fold": fold, "n_train": len(tr), "n_test": len(te),
                              "accuracy": accuracy_score(y[te], pred),
                              "balanced_accuracy": balanced_accuracy_score(y[te], pred),
                              "precision": precision_score(y[te], pred, zero_division=0),
                              "recall": recall_score(y[te], pred, zero_division=0),
                              "f1": f1_score(y[te], pred, zero_division=0),
                              "roc_auc": roc_auc_score(y[te], prob),
                              "pr_auc": average_precision_score(y[te], prob)})
        if not fold_rows:
            raise SystemExit("Grouped splits could not provide both classes. Add more balanced, independently sourced examples.")
        scores.extend({"model": name, **r} for r in fold_rows)
        final_models[name] = {k: float(np.mean([r[k] for r in fold_rows])) for k in fold_rows[0] if k not in {"fold", "n_train", "n_test"}}
    summary = {"evaluation": "5 repeated 80/20 GroupShuffleSplit holdouts; groups join >=80% normalized edit identity components and fragments from the same UniProt parent accession",
               "n_sequences": len(rows), "class_counts": {str(c): int(sum(y == c)) for c in (0,1)},
               "group_count": len(set(groups)), "mean_metrics_by_model": final_models,
               "limitations": ["Putative negatives may contain undiscovered AMPs.", "AMP positives and UniProt-derived negative fragments come from different biological sources; source or composition differences may inflate scores.", "No independent external test set; exploratory only. Do not interpret high CV scores as validated biological performance."]}
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "cv_metrics.csv").open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(scores[0])); w.writeheader(); w.writerows(scores)
    (OUT / "cv_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    for name, model in model_set().items(): model.fit(X, y); joblib.dump(model, OUT / f"{name}.joblib")
    (OUT / "model_metadata.json").write_text(json.dumps({"features": FEATURES, **summary}, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__": main()
