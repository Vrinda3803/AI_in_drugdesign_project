"""Exploratory unit-specific MIC regressors with peptide-grouped evaluation."""
from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

import joblib
import sklearn
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.linear_model import Ridge

ROOT = Path(__file__).resolve().parents[2]
MIC_FILE = ROOT / "data" / "processed" / "dbaasp_a_baumannii_mic_parsed.csv"
FEATURE_FILE = ROOT / "data" / "processed" / "sequence_features.csv"
OUT = ROOT / "results" / "activity_predictions"
MODELS = ROOT / "models" / "activity_regression"
SEED = 42
SPLITS = 5
TEST_SIZE = 0.2
TARGET_SPECIES = "Acinetobacter baumannii"  # species-level records; named isolates excluded for this first endpoint
UNITS = {"µg/ml": "ug_per_ml", "µM": "micromolar"}
FEATURES = ["length", "molecular_weight_average_unmodified_g_mol", "net_charge_unmodified_pH7",
            "isoelectric_point_unmodified", "gravy_kyte_doolittle",
            "hydrophobic_moment_eisenberg_100deg_full_sequence",
            "hydrophobic_moment_eisenberg_100deg_max_11mer"] + [f"AAC_{aa}" for aa in "ACDEFGHIKLMNPQRSTVWY"]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_targets(mic_rows: list[dict[str, str]]) -> dict[str, list[dict]]:
    # Use exact, parseable point values only for this simple baseline. Keep all
    # censored/range/uncertainty rows in the source dataset for later methods.
    bucket: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in mic_rows:
        unit = row["unit_raw"]
        if (row["target_species_raw"] != TARGET_SPECIES or unit not in UNITS or
                row["sequence_qc"] != "valid_standard_aa" or
                row["concentration_parse_status"] not in {"exact", "exact_comma_decimal"}):
            continue
        value = float(row["concentration_value"])
        if not math.isfinite(value) or value <= 0:
            continue
        key = (unit, row["sequence"])
        try:
            pmids = json.loads(row.get("candidate_pubmed_ids_json") or "[]")
        except json.JSONDecodeError:
            pmids = []
        bucket[key].append({**row, "candidate_pmids": pmids,
                            "concentration_numeric": value, "log2_mic": math.log2(value)})
    result: dict[str, list[dict]] = defaultdict(list)
    for (unit, sequence), measurements in bucket.items():
        logs = [x["log2_mic"] for x in measurements]
        result[unit].append({"unit_raw": unit, "unit_slug": UNITS[unit], "sequence": sequence,
                             "n_source_mic_observations": len(measurements),
                             "median_mic_in_source_unit": median(x["concentration_numeric"] for x in measurements),
                             "median_log2_mic": median(logs),
                             "dbaasp_ids_json": json.dumps(sorted({x["dbaasp_id"] for x in measurements})),
                             "activity_record_ids_json": json.dumps(sorted({x["activity_record_id"] for x in measurements})),
                             "candidate_pubmed_ids_json": json.dumps(sorted({pmid for x in measurements for pmid in x["candidate_pmids"]})),
                             "target_species_raw": TARGET_SPECIES})
    return result


def make_models() -> dict:
    return {
        "median_dummy": DummyRegressor(strategy="median"),
        "ridge": make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
        "svr_rbf": make_pipeline(StandardScaler(), SVR(C=1.0, epsilon=0.1, kernel="rbf")),
        "random_forest": RandomForestRegressor(n_estimators=500, min_samples_leaf=3, max_features=0.75,
                                               random_state=SEED, n_jobs=-1),
    }


def sequence_identity_clusters(sequences: list[str], threshold: float = 0.8) -> list[str]:
    """Connected components at global edit identity >= threshold.

    Identity is (max length - Levenshtein edit distance) / max length. This is
    an explicitly defined screening similarity, not a substitution-matrix
    protein alignment. Exact duplicates are already aggregated by sequence.
    """
    parent = list(range(len(sequences)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    def edit_distance(a: str, b: str) -> int:
        if len(a) < len(b):
            a, b = b, a
        previous = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            current = [i]
            for j, cb in enumerate(b, 1):
                current.append(min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (ca != cb)))
            previous = current
        return previous[-1]

    order = sorted(range(len(sequences)), key=lambda i: (len(sequences[i]), sequences[i]))
    for oi, i in enumerate(order):
        a = sequences[i]
        for j in order[oi + 1:]:
            b = sequences[j]
            maximum = max(len(a), len(b))
            if maximum and abs(len(a) - len(b)) / maximum > 1 - threshold:
                continue
            if (maximum - edit_distance(a, b)) / maximum >= threshold:
                union(i, j)
    roots = [find(i) for i in range(len(sequences))]
    root_labels: dict[int, str] = {}
    for i, root in enumerate(roots):
        root_labels[root] = min(sequences[i], root_labels.get(root, sequences[i]))
    return [root_labels[root] for root in roots]


def combine_similarity_and_publication_groups(sequences: list[str], similarity_groups: list[str],
                                              publication_ids: list[list[str]]) -> list[str]:
    """Join sequence-similarity components that share a candidate PubMed ID."""
    parent = list(range(len(sequences)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    first_by_similarity: dict[str, int] = {}
    first_by_pmid: dict[str, int] = {}
    for i, (similarity_group, pmids) in enumerate(zip(similarity_groups, publication_ids)):
        if similarity_group in first_by_similarity:
            union(i, first_by_similarity[similarity_group])
        else:
            first_by_similarity[similarity_group] = i
        for pmid in pmids:
            if pmid in first_by_pmid:
                union(i, first_by_pmid[pmid])
            else:
                first_by_pmid[pmid] = i
    roots = [find(i) for i in range(len(sequences))]
    labels: dict[int, str] = {}
    for i, root in enumerate(roots):
        labels[root] = min(sequences[i], labels.get(root, sequences[i]))
    return [labels[root] for root in roots]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    MODELS.mkdir(parents=True, exist_ok=True)
    mic_rows = read_csv(MIC_FILE)
    features = {row["sequence"]: row for row in read_csv(FEATURE_FILE)}
    cohorts = build_targets(mic_rows)
    target_rows = [row for unit_rows in cohorts.values() for row in unit_rows]
    target_fields = ["unit_raw", "unit_slug", "sequence", "n_source_mic_observations", "median_mic_in_source_unit",
                     "median_log2_mic", "dbaasp_ids_json", "activity_record_ids_json", "candidate_pubmed_ids_json", "target_species_raw"]
    write_csv(OUT / "mic_sequence_level_targets.csv", target_rows, target_fields)

    metric_rows: list[dict] = []
    prediction_rows: list[dict] = []
    run_summary: dict = {"endpoint": "median log2(MIC concentration) per exact peptide sequence, species-level target only",
                         "python_version": sys.version.split()[0], "scikit_learn_version": sklearn.__version__,
                         "target_species_filter": TARGET_SPECIES, "excluded_named_strains": True,
                         "unit_specific_models": True, "censoring_and_ranges": "excluded from this simple point-regression cohort but retained in parsed source table",
                         "sequence_qc": "only sequences containing the 20 standard amino acids",
                         "aggregation": "median over exact parseable MIC observations per exact sequence and unit; source rows/IDs retained in target table",
                         "splits": {"method": "GroupShuffleSplit", "test_size": TEST_SIZE, "repeats": SPLITS, "seed": SEED,
                                    "protocols": {"exact_sequence_group": "groups identical sequence records; one median target per sequence",
                                                  "global_edit_identity_0.8": "connected components at (max_length-Levenshtein_distance)/max_length >= 0.80; conservative sequence similarity screen, not a substitution-matrix alignment",
                                                  "identity_0.8_plus_candidate_publication": "combine edit-identity components and sequences sharing candidate PubMed IDs mapped from DBAASP reference indices; candidate reference mapping requires audit"}},
                         "feature_columns": FEATURES, "unit_results": {}, "metrics_interpretation": "descriptive mean and SD across overlapping repeated holdouts; not a final unbiased benchmark"}

    for unit, data in cohorts.items():
        data = [d for d in data if d["sequence"] in features]
        if len(data) < 30:
            run_summary["unit_results"][UNITS[unit]] = {"n_sequences": len(data), "status": "too few for planned holdouts"}
            continue
        X = [[float(features[d["sequence"]][feature]) for feature in FEATURES] for d in data]
        y = [float(d["median_log2_mic"]) for d in data]
        groups = [d["sequence"] for d in data]
        model_defs = make_models()
        homology_groups = sequence_identity_clusters(groups, threshold=0.8)
        publication_ids = [json.loads(d["candidate_pubmed_ids_json"]) for d in data]
        combined_groups = combine_similarity_and_publication_groups(groups, homology_groups, publication_ids)
        protocol_metrics: dict[str, dict[str, list[dict]]] = {}
        protocol_descriptions = {"exact_sequence_group": groups, "global_edit_identity_0.8": homology_groups,
                                 "identity_0.8_plus_candidate_publication": combined_groups}
        protocol_cluster_counts = {}
        for protocol, split_groups in protocol_descriptions.items():
            protocol_cluster_counts[protocol] = len(set(split_groups))
            protocol_metrics[protocol] = {name: [] for name in model_defs}
            split_iter = GroupShuffleSplit(n_splits=SPLITS, test_size=TEST_SIZE, random_state=SEED).split(X, y, split_groups)
            for fold, (train_idx, test_idx) in enumerate(split_iter, 1):
                for model_name, model in model_defs.items():
                    model.fit([X[i] for i in train_idx], [y[i] for i in train_idx])
                    pred = model.predict([X[i] for i in test_idx])
                    mae = mean_absolute_error([y[i] for i in test_idx], pred)
                    rmse = math.sqrt(mean_squared_error([y[i] for i in test_idx], pred))
                    r2 = r2_score([y[i] for i in test_idx], pred) if len(test_idx) > 1 else float("nan")
                    fold_row = {"unit_slug": UNITS[unit], "split_protocol": protocol, "model": model_name, "fold": fold,
                                "n_train_sequences": len(train_idx), "n_test_sequences": len(test_idx),
                                "MAE_log2_MIC": mae, "RMSE_log2_MIC": rmse, "R2_log2_MIC": r2,
                                "MAE_fold_change": 2 ** mae}
                    metric_rows.append(fold_row)
                    protocol_metrics[protocol][model_name].append(fold_row)
                    for idx, estimate in zip(test_idx, pred):
                        prediction_rows.append({"unit_slug": UNITS[unit], "split_protocol": protocol, "model": model_name, "fold": fold,
                                                "sequence": data[idx]["sequence"], "y_true_median_log2_mic": y[idx],
                                                "y_pred_log2_mic": float(estimate),
                                                "y_pred_mic_in_source_unit": 2 ** float(estimate),
                                                "n_source_mic_observations": data[idx]["n_source_mic_observations"]})
        # Store full-cohort fits for downstream candidate scoring; holdout
        # scores remain the evaluation and must not be conflated with these fits.
        for model_name, model in model_defs.items():
            model.fit(X, y)
            joblib.dump(model, MODELS / f"{UNITS[unit]}_{model_name}.joblib")
        run_summary["unit_results"][UNITS[unit]] = {
            "n_unique_sequences": len(data),
            "n_source_MIC_observations_before_sequence_aggregation": sum(int(d["n_source_mic_observations"]) for d in data),
            "n_sequences_with_repeated_MIC": sum(int(d["n_source_mic_observations"]) > 1 for d in data),
            "n_sequences_with_candidate_publication_ids": sum(bool(x) for x in publication_ids),
            "sequence_cluster_counts": protocol_cluster_counts,
            "models_by_split_protocol": {protocol: {
                name: {metric: {"mean": sum(float(r[metric]) for r in rows) / len(rows),
                                "sd": (sum((float(r[metric]) - sum(float(z[metric]) for z in rows)/len(rows))**2 for r in rows)/len(rows))**0.5}
                      for metric in ("MAE_log2_MIC", "RMSE_log2_MIC", "R2_log2_MIC", "MAE_fold_change")}
                for name, rows in model_metrics.items()}
                for protocol, model_metrics in protocol_metrics.items()}}

    write_csv(OUT / "mic_regression_cv_metrics_by_fold.csv", metric_rows,
              ["unit_slug", "split_protocol", "model", "fold", "n_train_sequences", "n_test_sequences", "MAE_log2_MIC", "RMSE_log2_MIC", "R2_log2_MIC", "MAE_fold_change"])
    write_csv(OUT / "mic_regression_cv_predictions.csv", prediction_rows,
              ["unit_slug", "split_protocol", "model", "fold", "sequence", "y_true_median_log2_mic", "y_pred_log2_mic", "y_pred_mic_in_source_unit", "n_source_mic_observations"])
    (OUT / "mic_regression_baseline_summary.json").write_text(json.dumps(run_summary, indent=2, ensure_ascii=False), encoding="utf-8")
    (MODELS / "metadata.json").write_text(json.dumps(run_summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({slug: info for slug, info in run_summary["unit_results"].items()}, ensure_ascii=True))


if __name__ == "__main__":
    main()
