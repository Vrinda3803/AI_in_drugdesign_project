"""Exploratory target- and unit-specific MIC regressions for 5 bacteria."""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PANEL = HERE / "output/multitarget_mic_evidence.csv"
FEATURES_PATH = ROOT / "data/processed/sequence_features.csv"
OUT = HERE / "output"
NUMERIC = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*$")
UNITS = {"µM": "micromolar", "µg/ml": "ug_per_ml", "µg/mL": "ug_per_ml"}
FEATURES = ["length", "molecular_weight_average_unmodified_g_mol", "net_charge_unmodified_pH7",
            "isoelectric_point_unmodified", "gravy_kyte_doolittle",
            "hydrophobic_moment_eisenberg_100deg_full_sequence",
            "hydrophobic_moment_eisenberg_100deg_max_11mer"] + [f"AAC_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"]


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def main():
    feature_map = {r["sequence"].upper(): r for r in read(FEATURES_PATH)}
    buckets = defaultdict(list)
    source_rows = read(PANEL)
    for r in source_rows:
        unit = UNITS.get(r["unit_raw"].strip())
        m = NUMERIC.fullmatch(r["mic_raw"] or "")
        if not unit or not m or r["sequence"].upper() not in feature_map:
            continue
        value = float(m.group(1))
        if value <= 0:
            continue
        buckets[(r["target_group"], unit, r["sequence"].upper())].append((math.log2(value), r))
    targets = defaultdict(list)
    for (target, unit, seq), observations in buckets.items():
        target_pmids = sorted({p for _, r in observations for p in json.loads(r["candidate_pubmed_ids_json"] or "[]")})
        targets[(target, unit)].append({"sequence": seq, "log2_mic": median(x[0] for x in observations),
                                        "n_observations": len(observations), "pmids": target_pmids})

    try:
        sys.path.insert(0, str(ROOT))
        from scripts.activity_prediction.mic_regression_baselines import sequence_identity_clusters, combine_similarity_and_publication_groups
    except ImportError as exc:
        raise SystemExit(f"Could not import existing grouped split utilities: {exc}")

    metric_rows, cohort_summary = [], {}
    for (target, unit), rows in sorted(targets.items()):
        if len(rows) < 30:
            cohort_summary[f"{target}|{unit}"] = {"status": "too few sequences for exploratory CV", "n_sequences": len(rows)}
            continue
        sequences = [r["sequence"] for r in rows]
        X = [[float(feature_map[s][k]) for k in FEATURES] for s in sequences]
        y = [r["log2_mic"] for r in rows]
        similarity = sequence_identity_clusters(sequences, 0.8)
        groups = combine_similarity_and_publication_groups(sequences, similarity, [r["pmids"] for r in rows])
        split = GroupShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
        model_metrics = defaultdict(list)
        for fold, (tr, te) in enumerate(split.split(X, y, groups), 1):
            for model_name, model in (("median_dummy", DummyRegressor(strategy="median")),
                                      ("random_forest", RandomForestRegressor(n_estimators=300, min_samples_leaf=3,
                                                                               max_features=0.75, random_state=42, n_jobs=-1))):
                model.fit([X[i] for i in tr], [y[i] for i in tr])
                pred = model.predict([X[i] for i in te])
                scores = {"mae_log2_mic": mean_absolute_error([y[i] for i in te], pred),
                          "rmse_log2_mic": mean_squared_error([y[i] for i in te], pred) ** 0.5,
                          "r2_log2_mic": r2_score([y[i] for i in te], pred)}
                metric_rows.append({"target_group": target, "unit": unit, "model": model_name,
                                    "fold": fold, "n_sequences": len(rows), "n_train": len(tr), "n_test": len(te), **scores})
                model_metrics[model_name].append(scores)
        cohort_summary[f"{target}|{unit}"] = {
            "status": "exploratory grouped CV", "n_sequences": len(rows),
            "group_components": len(set(groups)), "median_MIC_source_unit": 2 ** median(y),
            "mean_metrics": {name: {k: sum(v[k] for v in values) / len(values) for k in values[0]}
                            for name, values in model_metrics.items()},
        }
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "multitarget_mic_cv_metrics.csv").open("w", encoding="utf-8-sig", newline="") as f:
        fields = list(metric_rows[0]) if metric_rows else ["target_group", "unit", "model", "fold"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(metric_rows)
    report = {"scope": "Species-grouped, unit-separated exact-point MIC feasibility baselines from the cached A. baumannii-selected DBAASP peptide panel.",
              "endpoint": "Median log2 MIC by exact peptide sequence and broad target group; raw conditions remain in source evidence CSV.",
              "evaluation": "Five repeated grouped 80/20 holdouts; >=80% edit-identity components joined across candidate PubMed IDs.",
              "cohorts": cohort_summary,
              "limits": ["Other-organism records are only from the A. baumannii query-selected peptide subset.",
                         "Species groups combine strain/assay contexts and are not resistance-specific.",
                         "Only exact numeric MIC values are modeled; ranges and censored records are excluded.",
                         "Candidate PubMed mappings need source verification; there is no external validation set."]}
    (OUT / "multitarget_mic_cv_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"evaluated_cohorts": len(metric_rows), "cohorts": cohort_summary}, indent=2))


if __name__ == "__main__": main()
