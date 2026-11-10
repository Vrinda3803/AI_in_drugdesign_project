"""Extract multi-bacterium MIC evidence from cached DBAASP detail records."""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DETAILS = ROOT / "data/raw/DBAASP/a_baumannii_details"
OUT = Path(__file__).resolve().parent / "output"
AA = set("ACDEFGHIKLMNPQRSTVWY")
TARGETS = {
    "Acinetobacter baumannii": "A. baumannii",
    "Escherichia coli": "E. coli",
    "Pseudomonas aeruginosa": "P. aeruginosa",
    "Klebsiella pneumoniae": "K. pneumoniae",
    "Staphylococcus aureus": "S. aureus",
}


def obj_name(value):
    return value.get("name", "") if isinstance(value, dict) else str(value or "")


def main():
    rows, counts, units = [], Counter(), Counter()
    for path in sorted(DETAILS.glob("*.json"), key=lambda x: int(x.stem)):
        try:
            peptide = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        seq = (peptide.get("sequence") or "").strip()
        if not seq or not set(seq.upper()) <= AA:
            continue
        articles = peptide.get("articles") or []
        for assay in peptide.get("targetActivities") or []:
            measure = obj_name(assay.get("activityMeasureGroup")).strip()
            if measure.casefold() != "mic":
                continue
            target = obj_name(assay.get("targetSpecies")).strip()
            target_group = next((group for prefix, group in TARGETS.items() if target.casefold().startswith(prefix.casefold())), None)
            if target_group is None:
                continue
            ref = str(assay.get("reference") or "").strip()
            ref_ids = []
            for token in ref.replace(";", ",").split(","):
                try:
                    idx = int(token.strip()) - 1
                    if 0 <= idx < len(articles):
                        pmid = (articles[idx].get("pubmed") or {}).get("pubmedId")
                        if pmid: ref_ids.append(str(pmid))
                except ValueError:
                    continue
            unit = obj_name(assay.get("unit"))
            row = {
                "dbaasp_id": peptide.get("id", path.stem),
                "name": peptide.get("name", ""), "sequence": seq,
                "target_group": target_group, "target_species_raw": target,
                "mic_raw": str(assay.get("concentration") or ""), "unit_raw": unit,
                "medium_raw": obj_name(assay.get("medium")), "strain_or_isolate_note_raw": assay.get("note", ""),
                "ph_raw": assay.get("ph", ""), "ionic_strength_raw": assay.get("ionicStrength", ""),
                "cfu_raw": assay.get("cfu", ""), "reference_raw": ref,
                "candidate_pubmed_ids_json": json.dumps(ref_ids),
                "label_caveat": "Source-reported MIC; not harmonized across strains or protocols; resistance status only if supported by source evidence.",
            }
            rows.append(row); counts[target_group] += 1; units[(target_group, unit)] += 1
    OUT.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with (OUT / "multitarget_mic_evidence.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    report = {
        "records": len(rows), "unique_peptide_ids": len({r["dbaasp_id"] for r in rows}),
        "counts_by_target_group": dict(counts),
        "counts_by_target_and_unit": {f"{k[0]} | {k[1]}": v for k, v in units.items()},
        "source_scope": "Cached DBAASP detail records downloaded for the A. baumannii query; other-target counts are a selected peptide subset, not a comprehensive species database.",
        "limits": ["MICs remain raw and unit-separated.", "Species labels do not prove MDR status.", "Same-peptide cross-target comparisons need strain and assay harmonization.", "Candidate publication mappings are not fully source-paper verified."],
    }
    (OUT / "multitarget_mic_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__": main()
