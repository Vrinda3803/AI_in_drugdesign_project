"""Build provenance-preserving analysis tables from the acquired raw sources."""
from __future__ import annotations

import csv
import json
import re
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
DETAILS = RAW / "DBAASP" / "a_baumannii_details"
AA = set("ACDEFGHIKLMNPQRSTVWY")
ABAMP_CSV = RAW / "AbAMPdb" / "ampdb-Final.csv"
TARGET_RE = re.compile(r"^Acinetobacter\s+baumannii\b", re.I)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="cp1252", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def seq_qc(raw: str) -> tuple[str, str, str]:
    sequence = (raw or "").strip().upper()
    invalid = "".join(sorted(set(sequence) - AA))
    status = "valid_standard_aa" if sequence and not invalid else ("missing" if not sequence else "nonstandard_or_invalid")
    return sequence, status, invalid


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    abamp_raw = read_csv(ABAMP_CSV)
    abamp_rows: list[dict] = []
    for row in abamp_raw:
        sequence, status, invalid = seq_qc(row.get("Sequence", ""))
        abamp_rows.append({**row, "Sequence": sequence, "sequence_qc": status, "invalid_characters": invalid,
                           "source": "AbAMPdb", "source_record_id": row.get("AMP", "")})

    peptide_rows: list[dict] = []
    activity_rows: list[dict] = []
    failures: list[dict] = []
    details_found = 0
    species_counts: Counter[str] = Counter()
    assay_counts: Counter[str] = Counter()
    for path in sorted(DETAILS.glob("*.json"), key=lambda p: int(p.stem)):
        try:
            peptide = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            failures.append({"file": path.name, "error": str(exc)})
            continue
        details_found += 1
        seq, status, invalid = seq_qc(peptide.get("sequence", ""))
        monomers = peptide.get("monomers") or []
        articles = peptide.get("articles") or []
        mono = [seq_qc(m.get("sequence", ""))[0] for m in monomers]
        peptide_rows.append({"source": "DBAASP", "dbaasp_id": peptide.get("id"), "dbaasp_accession": peptide.get("dbaaspId"),
                             "name": peptide.get("name"), "sequence": seq, "sequence_qc": status, "invalid_characters": invalid,
                             "sequence_length_reported": peptide.get("sequenceLength"), "complexity": peptide.get("complexity", {}).get("name") if isinstance(peptide.get("complexity"), dict) else peptide.get("complexity"),
                             "monomer_sequences_json": json.dumps(mono, ensure_ascii=False),
                             "n_terminus": (peptide.get("nTerminus") or peptide.get("NTerminus") or {}).get("name") if isinstance(peptide.get("nTerminus") or peptide.get("NTerminus"), dict) else peptide.get("nTerminus") or peptide.get("NTerminus"),
                             "c_terminus": (peptide.get("cTerminus") or peptide.get("CTerminus") or {}).get("name") if isinstance(peptide.get("cTerminus") or peptide.get("CTerminus"), dict) else peptide.get("cTerminus") or peptide.get("CTerminus")})
        for activity in peptide.get("targetActivities") or []:
            species_obj = activity.get("targetSpecies") or {}
            species = species_obj.get("name", "") if isinstance(species_obj, dict) else str(species_obj)
            if not TARGET_RE.match(species):
                continue
            measure_obj = activity.get("activityMeasureGroup") or {}
            measure = measure_obj.get("name", "") if isinstance(measure_obj, dict) else str(measure_obj)
            unit_obj = activity.get("unit") or {}
            unit = unit_obj.get("name", "") if isinstance(unit_obj, dict) else str(unit_obj)
            medium_obj = activity.get("medium") or {}
            medium = medium_obj.get("name", "") if isinstance(medium_obj, dict) else str(medium_obj)
            reference_raw = str(activity.get("reference") or "").strip()
            reference_indices = [int(x) for x in re.findall(r"\d+", reference_raw)]
            article_matches = [articles[i - 1] for i in reference_indices if 1 <= i <= len(articles)]
            if reference_indices and len(article_matches) == len(reference_indices):
                reference_status = "candidate_1_based_article_index_mapping"
            elif article_matches:
                reference_status = "partially_mapped"
            else:
                reference_status = "unresolved"
            pubmed_ids = [str((article.get("pubmed") or {}).get("pubmedId") or "") for article in article_matches]
            article_titles = [str(article.get("title") or "") for article in article_matches]
            activity_rows.append({"source": "DBAASP", "dbaasp_id": peptide.get("id"), "dbaasp_accession": peptide.get("dbaaspId"),
                                  "name": peptide.get("name"), "sequence": seq, "sequence_qc": status,
                                  "activity_record_id": activity.get("id"), "target_species_raw": species,
                                  "activity_measure": measure, "activity_measure_value_raw": activity.get("activityMeasureValue"),
                                  "concentration_raw": activity.get("concentration"), "unit_raw": unit,
                                  "dbaasp_activity_numeric_raw": activity.get("activity"), "strain_or_isolate_note_raw": activity.get("note"),
                                  "medium_raw": medium, "ph_raw": activity.get("ph"), "ionic_strength_raw": activity.get("ionicStrength"),
                                  "salt_type_raw": activity.get("saltType"), "cfu_raw": activity.get("cfu"),
                                  "reference_raw": activity.get("reference"),
                                  "reference_resolution_status": reference_status,
                                  "reference_article_indices_json": json.dumps(reference_indices),
                                  "candidate_pubmed_ids_json": json.dumps([x for x in pubmed_ids if x]),
                                  "candidate_reference_titles_json": json.dumps(article_titles, ensure_ascii=False)})
            species_counts[species] += 1
            assay_counts[measure] += 1

    ab_fields = list(abamp_rows[0]) if abamp_rows else []
    write_csv(OUT / "abampdb_records_qc.csv", abamp_rows, ab_fields)
    activity_fields = ["source", "dbaasp_id", "dbaasp_accession", "name", "sequence", "sequence_qc", "activity_record_id", "target_species_raw", "activity_measure", "activity_measure_value_raw", "concentration_raw", "unit_raw", "dbaasp_activity_numeric_raw", "strain_or_isolate_note_raw", "medium_raw", "ph_raw", "ionic_strength_raw", "salt_type_raw", "cfu_raw", "reference_raw", "reference_resolution_status", "reference_article_indices_json", "candidate_pubmed_ids_json", "candidate_reference_titles_json"]
    write_csv(OUT / "dbaasp_a_baumannii_activity_records.csv", activity_rows, activity_fields)
    mic_rows = [r for r in activity_rows if str(r["activity_measure"]).strip().casefold() == "mic"]
    write_csv(OUT / "dbaasp_a_baumannii_mic_measurements.csv", mic_rows, activity_fields)

    seq_groups: dict[str, list[str]] = defaultdict(list)
    for row in abamp_rows:
        if row["Sequence"]:
            seq_groups[row["Sequence"]].append(f"AbAMPdb:{row['AMP']}")
    for row in peptide_rows:
        if row["sequence"]:
            seq_groups[row["sequence"]].append(f"DBAASP:{row['dbaasp_id']}")
    sequence_panel = [{"sequence": seq, "length": len(seq), "sources_and_ids_json": json.dumps(ids),
                       "source_record_count": len(ids), "cross_source_or_within_source_duplicate": len(ids) > 1,
                       "sequence_qc": seq_qc(seq)[1], "invalid_characters": seq_qc(seq)[2]}
                      for seq, ids in sorted(seq_groups.items())]
    write_csv(OUT / "sequence_provenance_panel.csv", sequence_panel,
              ["sequence", "length", "sources_and_ids_json", "source_record_count", "cross_source_or_within_source_duplicate", "sequence_qc", "invalid_characters"])

    valid_ab = [r for r in abamp_rows if r["sequence_qc"] == "valid_standard_aa"]
    valid_dbaasp = [r for r in peptide_rows if r["sequence_qc"] == "valid_standard_aa"]
    mic_ids = {str(r["dbaasp_id"]) for r in mic_rows}
    mic_sequences = {r["sequence"] for r in mic_rows if r["sequence"]}
    mic_valid_sequences = {r["sequence"] for r in mic_rows if r["sequence_qc"] == "valid_standard_aa"}
    exact_ab_dups = len(valid_ab) - len({r["Sequence"] for r in valid_ab})
    exact_dbaasp_dups = len(valid_dbaasp) - len({r["sequence"] for r in valid_dbaasp})
    report = {
        "scope": "Source-level audit and raw activity extraction; no activity thresholding or MIC unit conversion applied.",
        "abampdb": {"source_rows": len(abamp_rows), "valid_standard_aa_sequences": len(valid_ab),
                    "missing_or_invalid_sequences": len(abamp_rows) - len(valid_ab), "exact_duplicate_rows_among_valid": exact_ab_dups,
                    "MIC_value_units_known_from_export": False},
        "dbaasp_a_baumannii": {"search_index_records": 1126, "detail_responses_found": details_found,
                                "activity_records_species_prefix_match": len(activity_rows), "MIC_measurement_rows": len(mic_rows),
                                "unique_peptide_ids_with_MIC": len(mic_ids),
                                "unique_nonempty_raw_sequences_among_MIC_records": len(mic_sequences),
                                "unique_standard_AA_sequences_among_MIC_records": len(mic_valid_sequences),
                                "peptide_detail_records_with_valid_parent_sequence": len(valid_dbaasp),
                                "exact_duplicate_rows_among_valid_parent_sequences": exact_dbaasp_dups,
                                "activity_record_counts_by_measure": dict(assay_counts),
                                "activity_reference_resolution_counts": dict(Counter(r["reference_resolution_status"] for r in activity_rows)),
                                "activity_record_counts_by_raw_species_name": dict(species_counts),
                                "missing_or_invalid_parent_sequences": len(peptide_rows) - len(valid_dbaasp),
                                "detail_failures_while_parsing": failures},
        "combined_sequence_panel": {"unique_nonempty_parent_sequences": len(sequence_panel),
                                     "panel_records_with_duplicate_sequence": sum(x["cross_source_or_within_source_duplicate"] for x in sequence_panel)},
        "limitations": ["Sequence deduplication is exact only at the parent sequence level; multimers/modified peptides require separate treatment.",
                        "AbAMPdb exported MIC units and source citations are not explicit in this CSV and are not merged into DBAASP MIC values.",
                        "DBAASP activity values remain in source strings and units; censored/range parsing is deferred, not guessed."]
    }
    (OUT / "dataset_qc_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    # Keep console output compatible with the default Windows code page.
    print(json.dumps(report, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
