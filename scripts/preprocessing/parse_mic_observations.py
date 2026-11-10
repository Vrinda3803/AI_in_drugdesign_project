"""Parse MIC concentration syntax without collapsing ranges or censoring."""
from __future__ import annotations

import csv
import json
import re
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data" / "processed" / "dbaasp_a_baumannii_mic_measurements.csv"
OUTPUT = ROOT / "data" / "processed" / "dbaasp_a_baumannii_mic_parsed.csv"
NUMBER = r"(?:\d+(?:\.\d*)?|\.\d+)"
EXACT_RE = re.compile(rf"^\s*({NUMBER})\s*$")
COMMA_DECIMAL_RE = re.compile(r"^\s*(\d+),(\d+)\s*$")
CENSOR_RE = re.compile(rf"^\s*(<=|>=|<|>|≤|≥)\s*({NUMBER})\s*$")
RANGE_RE = re.compile(rf"^\s*({NUMBER})\s*[-–—]\s*({NUMBER})\s*$")
UNCERTAINTY_RE = re.compile(rf"^\s*({NUMBER})\s*(?:±|\+/-)\s*({NUMBER})\s*$")
TRANSITION_RE = re.compile(rf"^\s*({NUMBER})\s*->\s*({NUMBER})\s*$")


def parse(raw: str) -> dict[str, str]:
    text = (raw or "").strip().replace("−", "-")
    result = {"concentration_parse_status": "unparsed", "concentration_operator": "",
              "concentration_value": "", "concentration_lower": "", "concentration_upper": "",
              "reported_uncertainty": "", "reported_transition_start": "", "reported_transition_end": ""}
    if not text:
        result["concentration_parse_status"] = "missing"
        return result
    try:
        match = EXACT_RE.fullmatch(text)
        if match:
            value = Decimal(match.group(1))
            result.update(concentration_parse_status="exact", concentration_operator="=",
                          concentration_value=str(value), concentration_lower=str(value), concentration_upper=str(value))
            return result
        match = COMMA_DECIMAL_RE.fullmatch(text)
        if match:
            value = Decimal(".".join(match.groups()))
            result.update(concentration_parse_status="exact_comma_decimal", concentration_operator="=",
                          concentration_value=str(value), concentration_lower=str(value), concentration_upper=str(value))
            return result
        match = CENSOR_RE.fullmatch(text)
        if match:
            operator, number = match.groups()
            canonical = {"≤": "<=", "≥": ">="}.get(operator, operator)
            value = str(Decimal(number))
            result.update(concentration_parse_status="censored", concentration_operator=canonical,
                          concentration_value=value)
            if canonical in ("<", "<="):
                result["concentration_upper"] = value
            else:
                result["concentration_lower"] = value
            return result
        match = RANGE_RE.fullmatch(text)
        if match:
            low, high = map(Decimal, match.groups())
            if low <= high:
                result.update(concentration_parse_status="range", concentration_operator="range",
                              concentration_lower=str(low), concentration_upper=str(high))
                return result
        match = TRANSITION_RE.fullmatch(text)
        if match:
            start, end = map(Decimal, match.groups())
            result.update(concentration_parse_status="reported_transition_not_interpreted",
                          reported_transition_start=str(start), reported_transition_end=str(end))
            return result
        match = UNCERTAINTY_RE.fullmatch(text)
        if match:
            center, uncertainty = map(Decimal, match.groups())
            result.update(concentration_parse_status="reported_uncertainty_not_interpreted",
                          concentration_value=str(center), reported_uncertainty=str(uncertainty))
            return result
    except InvalidOperation:
        pass
    return result


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    counts: Counter[str] = Counter()
    output: list[dict[str, str]] = []
    for row in rows:
        parsed = parse(row.get("concentration_raw", ""))
        counts[parsed["concentration_parse_status"]] += 1
        output.append({**row, **parsed})
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output[0]) if output else []
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output)
    report_path = ROOT / "data" / "processed" / "dataset_qc_report.json"
    report = json.loads(report_path.read_text(encoding="utf-8-sig"))
    report["dbaasp_a_baumannii"]["MIC_concentration_parse_status_counts"] = dict(counts)
    report["dbaasp_a_baumannii"]["MIC_concentration_units_raw_counts"] = dict(Counter(r.get("unit_raw", "") for r in rows))
    report["dbaasp_a_baumannii"]["MIC_parsing_rule"] = "Decimal syntax only; ranges and censoring are represented explicitly; uncertainty notation is preserved without interpreting its statistical meaning; no unit conversion or activity threshold applied."
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"rows": len(output), "parse_status_counts": dict(counts)}, ensure_ascii=True))


if __name__ == "__main__":
    main()
