"""Download full DBAASP details for peptides in the A. baumannii search result.

Raw detail responses are cached by DBAASP numeric peptide ID. Existing files are
reused so interrupted downloads can be resumed without repeating requests.
"""
from __future__ import annotations

import concurrent.futures
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "DBAASP"
DETAILS = RAW / "a_baumannii_details"
PAGES = [RAW / "a_baumannii_search_page_0000.json", RAW / "a_baumannii_search_page_1000.json"]
API = "https://dbaasp.org/peptides/{}"


def load_ids() -> list[int]:
    ids: set[int] = set()
    for page in PAGES:
        payload = json.loads(page.read_text(encoding="utf-8-sig"))
        ids.update(int(item["id"]) for item in payload["data"])
    return sorted(ids)


def download_one(peptide_id: int) -> tuple[int, bool, str]:
    target = DETAILS / f"{peptide_id}.json"
    if target.exists() and target.stat().st_size:
        try:
            json.loads(target.read_text(encoding="utf-8"))
            return peptide_id, True, "cached"
        except (json.JSONDecodeError, UnicodeDecodeError):
            target.unlink()
    request = urllib.request.Request(API.format(peptide_id), headers={"User-Agent": "AMP-research-data-acquisition/1.0"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                raw = response.read()
            # Validate before writing so interrupted/invalid responses are not cached.
            payload = json.loads(raw)
            if int(payload.get("id", -1)) != peptide_id:
                raise ValueError(f"Unexpected response ID for {peptide_id}")
            target.write_bytes(raw)
            return peptide_id, True, "downloaded"
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, ValueError) as exc:
            if attempt == 4:
                return peptide_id, False, str(exc)
            time.sleep(2**attempt)
    return peptide_id, False, "unreachable"


def main() -> None:
    DETAILS.mkdir(parents=True, exist_ok=True)
    ids = load_ids()
    print(f"Unique DBAASP IDs in A. baumannii search: {len(ids)}")
    failures: list[tuple[int, str]] = []
    completed = 0
    # Low concurrency and retry/backoff to respect the public API.
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for peptide_id, success, message in pool.map(download_one, ids):
            completed += 1
            if not success:
                failures.append((peptide_id, message))
            if completed % 100 == 0 or completed == len(ids):
                print(f"{completed}/{len(ids)} details handled; failures={len(failures)}")
    if failures:
        (RAW / "a_baumannii_detail_failures.json").write_text(
            json.dumps([{"id": peptide_id, "error": message} for peptide_id, message in failures], indent=2),
            encoding="utf-8",
        )
        raise SystemExit(f"Incomplete detail download: {len(failures)} failures; rerun to resume")
    print(f"All {len(ids)} detail responses are present and valid in {DETAILS}")


if __name__ == "__main__":
    main()
