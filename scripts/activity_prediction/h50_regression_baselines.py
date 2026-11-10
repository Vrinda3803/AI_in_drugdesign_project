"""Exploratory, unit-specific H50 regression feasibility analysis.

Uses only exact reported human-RBC 50%-hemolysis concentrations and standard
amino-acid sequences. It does not define safe/unsafe labels or predict general
toxicity. This small, assay-heterogeneous cohort is not a validated model.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

import sklearn
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mic_regression_baselines as shared

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data/processed/dbaasp_human_rbc_h50_candidate_records.csv"
FEATURE_FILE = ROOT / "data/processed/sequence_features.csv"
OUT = ROOT / "results/safety_predictions"
SEED, SPLITS, TEST_SIZE = 42, 5, 0.2
UNITS = {"µM": "micromolar", "µg/ml": "ug_per_ml"}
FEATURES = shared.FEATURES


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def estimators():
    return {
        "median_dummy": DummyRegressor(strategy="median"),
        "ridge": make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
        "svr_rbf": make_pipeline(StandardScaler(), SVR(C=1.0, epsilon=0.1, kernel="rbf")),
        "random_forest": RandomForestRegressor(n_estimators=400, min_samples_leaf=3,
                                                max_features=0.75, random_state=SEED, n_jobs=-1),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = read_csv(INPUT)
    feature_map = {row["sequence"]: row for row in read_csv(FEATURE_FILE)}
    buckets: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if (row["sequence_qc"] != "valid_standard_aa" or
                row["concentration_parse_status"] not in {"exact", "exact_comma_decimal"} or
                row["unit_raw"] not in UNITS):
            continue
        value = float(row["concentration_value"])
        if value > 0 and math.isfinite(value):
            buckets[(row["unit_raw"], row["sequence"])].append(row)

    targets = []
    metrics = []
    predictions = []
    summary = {
        "endpoint": "log2 concentration reported for exact 50% human erythrocyte hemolysis",
        "scope_limit": "This is an exploratory assay-specific hemolysis endpoint, not general toxicity or a safe/unsafe classification.",
        "selection": "Human erythrocytes; exact 50% Hemolysis endpoint; exact numeric concentration; standard amino acids; units separate.",
        "censoring": "Ranges, censored values, and uncertainty notation are excluded from this baseline and retained in the source cohort.",
        "sequence_scope": "One-letter sequence descriptors; terminal and other chemical modifications are not represented in model features.",
        "publication_grouping": "Candidate PubMed IDs are mapped through the source's article-index field and have not been comprehensively paper-verified.",
        "evaluation": {"split": "GroupShuffleSplit", "test_size": TEST_SIZE, "repeats": SPLITS, "seed": SEED,
                       "groups": "connected 80% global edit-identity components joined across shared candidate PubMed IDs"},
        "python_version": sys.version.split()[0], "scikit_learn_version": sklearn.__version__, "units": {}
    }
    targets_by_unit: dict[str, list[dict]] = defaultdict(list)
    for (unit, sequence), observations in buckets.items():
        values = [float(row["concentration_value"]) for row in observations]
        pmids = sorted({row["candidate_pubmed_id"].strip() for row in observations
                        if row.get("candidate_pubmed_id", "").strip()})
        target = {"unit_raw": unit, "unit_slug": UNITS[unit], "sequence": sequence,
                  "n_observations": len(observations), "median_h50": median(values),
                  "median_log2_h50": median(math.log2(value) for value in values),
                  "dbaasp_ids_json": json.dumps(sorted({row["dbaasp_id"] for row in observations})),
                  "candidate_pubmed_ids_json": json.dumps(pmids)}
        targets.append(target)
        if sequence in feature_map:
            targets_by_unit[unit].append(target)

    for unit in UNITS:
        unit_targets = targets_by_unit[unit]
        all_unit_targets = [row for row in targets if row["unit_raw"] == unit]
        if len(unit_targets) < 10:
            summary["units"][UNITS[unit]] = {
                "n_candidate_sequences_before_feature_join": len(all_unit_targets),
                "n_candidate_observations_before_feature_join": sum(int(row["n_observations"]) for row in all_unit_targets),
                "n_sequences_excluded_for_missing_exact_feature_match": len(all_unit_targets) - len(unit_targets),
                "n_sequences": len(unit_targets), "status": "too few for planned grouped holdouts"}
            continue
        sequences = [row["sequence"] for row in unit_targets]
        X = [[float(feature_map[row["sequence"]][feature]) for feature in FEATURES] for row in unit_targets]
        y = [float(row["median_log2_h50"]) for row in unit_targets]
        identity = shared.sequence_identity_clusters(sequences, 0.8)
        pubids = [json.loads(row["candidate_pubmed_ids_json"]) for row in unit_targets]
        groups = shared.combine_similarity_and_publication_groups(sequences, identity, pubids)
        if len(set(groups)) < 2:
            summary["units"][UNITS[unit]] = {"n_sequences": len(unit_targets), "n_group_components": len(set(groups)),
                                             "status": "too few independent groups for grouped holdouts"}
            continue

        model_rows: dict[str, list[dict]] = {name: [] for name in estimators()}
        split_iter = GroupShuffleSplit(n_splits=SPLITS, test_size=TEST_SIZE, random_state=SEED).split(X, y, groups)
        for fold, (train_idx, test_idx) in enumerate(split_iter, 1):
            for name, estimator in estimators().items():
                estimator.fit([X[i] for i in train_idx], [y[i] for i in train_idx])
                predicted = estimator.predict([X[i] for i in test_idx])
                actual = [y[i] for i in test_idx]
                row = {"unit_slug": UNITS[unit], "fold": fold, "model": name,
                       "n_train_sequences": len(train_idx), "n_test_sequences": len(test_idx),
                       "MAE_log2_H50": mean_absolute_error(actual, predicted),
                       "RMSE_log2_H50": math.sqrt(mean_squared_error(actual, predicted)),
                       "R2_log2_H50": r2_score(actual, predicted) if len(actual) > 1 else None}
                row["MAE_fold_change"] = 2 ** row["MAE_log2_H50"]
                metrics.append(row)
                model_rows[name].append(row)
                for idx, estimate in zip(test_idx, predicted):
                    predictions.append({"unit_slug": UNITS[unit], "split_protocol": "identity_0.8_plus_candidate_publication",
                                        "model": name, "fold": fold, "sequence": sequences[idx],
                                        "y_true_log2_H50": y[idx], "y_pred_log2_H50": float(estimate),
                                        "y_pred_H50": 2 ** float(estimate), "H50_unit": unit,
                                        "n_source_observations": unit_targets[idx]["n_observations"]})

        summary["units"][UNITS[unit]] = {
            "n_candidate_sequences_before_feature_join": len(all_unit_targets),
            "n_candidate_observations_before_feature_join": sum(int(row["n_observations"]) for row in all_unit_targets),
            "n_sequences_excluded_for_missing_exact_feature_match": len(all_unit_targets) - len(unit_targets),
            "n_unique_sequences": len(unit_targets),
            "n_source_observations": sum(int(row["n_observations"]) for row in unit_targets),
            "n_group_components": len(set(groups)),
            "models": {name: {metric: {"mean": sum(row[metric] for row in fold_rows if row[metric] is not None) /
                                                              sum(row[metric] is not None for row in fold_rows),
                                        "sd": (sum((row[metric] - sum(other[metric] for other in fold_rows if other[metric] is not None) /
                                                                  sum(other[metric] is not None for other in fold_rows)) ** 2
                                                   for row in fold_rows if row[metric] is not None) /
                                                sum(row[metric] is not None for row in fold_rows)) ** 0.5}
                              for metric in ("MAE_log2_H50", "RMSE_log2_H50", "R2_log2_H50", "MAE_fold_change")}
                       for name, fold_rows in model_rows.items()}}

    write_csv(OUT / "human_rbc_h50_sequence_targets.csv", targets,
              ["unit_raw", "unit_slug", "sequence", "n_observations", "median_h50", "median_log2_h50",
               "dbaasp_ids_json", "candidate_pubmed_ids_json"])
    write_csv(OUT / "human_rbc_h50_grouped_metrics_by_fold.csv", metrics,
              ["unit_slug", "fold", "model", "n_train_sequences", "n_test_sequences", "MAE_log2_H50",
               "RMSE_log2_H50", "R2_log2_H50", "MAE_fold_change"])
    write_csv(OUT / "human_rbc_h50_grouped_predictions.csv", predictions,
              ["unit_slug", "split_protocol", "model", "fold", "sequence", "y_true_log2_H50",
               "y_pred_log2_H50", "y_pred_H50", "H50_unit", "n_source_observations"])
    (OUT / "human_rbc_h50_baseline_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary["units"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
