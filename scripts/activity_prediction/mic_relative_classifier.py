"""Exploratory relative MIC classification with peptide/publication grouped splits.

Labels mean below vs. at/above the cohort median MIC concentration within each
unit. They are not experimentally validated active/inactive labels.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import sklearn
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, average_precision_score, f1_score,
                             precision_score, recall_score, roc_auc_score)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

import mic_regression_baselines as base

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "activity_classification"
FEATURES = base.FEATURES
SPLITS, TEST_SIZE, SEED = 5, 0.2, 42


def models():
    return {
        "majority_dummy": DummyClassifier(strategy="most_frequent"),
        "logistic_regression": make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=3000, random_state=SEED)),
        "svm_rbf": make_pipeline(StandardScaler(), SVC(C=1.0, kernel="rbf", class_weight="balanced", random_state=SEED)),
        "random_forest": RandomForestClassifier(n_estimators=400, min_samples_leaf=3, class_weight="balanced_subsample", random_state=SEED, n_jobs=-1),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mic = base.read_csv(base.MIC_FILE)
    feature_map = {r["sequence"]: r for r in base.read_csv(base.FEATURE_FILE)}
    cohorts = base.build_targets(mic)
    metrics, summary = [], {
        "endpoint": "relative potency class from exact species-level MIC point measurements",
        "label_definition": "lower_mic_than_unit_cohort_median vs at_or_above_unit_cohort_median, per unit",
        "not_an_activity_call": "No validated inactive labels are present; these labels do not mean active/inactive.",
        "evaluation": {"split": "GroupShuffleSplit", "test_size": TEST_SIZE, "repeats": SPLITS, "seed": SEED,
                       "groups": "connected 80% global edit-identity components joined across shared candidate PubMed IDs"},
        "python_version": sys.version.split()[0], "scikit_learn_version": sklearn.__version__, "units": {}
    }
    for unit, rows in cohorts.items():
        rows = [r for r in rows if r["sequence"] in feature_map]
        if len(rows) < 30:
            summary["units"][base.UNITS[unit]] = {"status": "insufficient cohort", "n_sequences": len(rows)}
            continue
        cut = sorted(float(r["median_log2_mic"]) for r in rows)[len(rows)//2]
        X = [[float(feature_map[r["sequence"]][f]) for f in FEATURES] for r in rows]
        y = [int(float(r["median_log2_mic"]) >= cut) for r in rows]
        seqs = [r["sequence"] for r in rows]
        sim = base.sequence_identity_clusters(seqs, 0.8)
        pubs = [json.loads(r["candidate_pubmed_ids_json"]) for r in rows]
        groups = base.combine_similarity_and_publication_groups(seqs, sim, pubs)
        by_model = {name: [] for name in models()}
        split_iter = GroupShuffleSplit(n_splits=SPLITS, test_size=TEST_SIZE, random_state=SEED).split(X, y, groups)
        for fold, (train, test) in enumerate(split_iter, 1):
            for name, estimator in models().items():
                estimator.fit([X[i] for i in train], [y[i] for i in train])
                pred = estimator.predict([X[i] for i in test])
                X_test = [X[i] for i in test]
                scores = (estimator.decision_function(X_test) if hasattr(estimator, "decision_function")
                          else estimator.predict_proba(X_test)[:, 1])
                yt = [y[i] for i in test]
                row = {"unit_slug": base.UNITS[unit], "fold": fold, "model": name,
                       "n_train": len(train), "n_test": len(test), "accuracy": accuracy_score(yt, pred),
                       "precision": precision_score(yt, pred, zero_division=0),
                       "recall": recall_score(yt, pred, zero_division=0), "f1": f1_score(yt, pred, zero_division=0),
                       "roc_auc": roc_auc_score(yt, scores) if len(set(yt)) == 2 else None,
                       "pr_auc": average_precision_score(yt, scores) if len(set(yt)) == 2 else None}
                metrics.append(row); by_model[name].append(row)
        summary["units"][base.UNITS[unit]] = {
            "n_sequences": len(rows), "n_positive_lower_mic": sum(v == 0 for v in y),
            "n_at_or_above_median": sum(y), "median_cutoff_mic_source_unit": 2 ** cut,
            "n_group_components": len(set(groups)),
            "metrics_mean_sd": {name: {metric: {"mean": sum(r[metric] for r in values)/len(values),
                                                   "sd": (sum((r[metric]-sum(z[metric] for z in values)/len(values))**2 for r in values)/len(values))**0.5}
                                         for metric in ("accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc")}
                              for name, values in by_model.items()}}
    base.write_csv(OUT/"mic_relative_classification_metrics_by_fold.csv", metrics,
                   ["unit_slug", "fold", "model", "n_train", "n_test", "accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc"])
    (OUT/"mic_relative_classification_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary["units"], indent=2))


if __name__ == "__main__":
    main()
