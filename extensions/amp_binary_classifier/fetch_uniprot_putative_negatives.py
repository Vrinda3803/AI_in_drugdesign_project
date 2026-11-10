"""Fetch length-matched, provenance-tracked putative-negative fragments.

These are not experimentally confirmed non-AMPs. Each is a deterministic
random window from a reviewed E. coli protein without obvious AMP/defense
annotation; isolated fragments may still have activity.
"""
from __future__ import annotations

import csv
import json
import random
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "data" / "non_amp_sequences.csv"
MANIFEST = HERE / "data" / "uniprot_negative_dataset_manifest.json"
QUERY = "reviewed:true AND organism_id:562 AND length:[101 TO 1000]"
FIELDS = "accession,protein_name,organism_name,length,sequence,keyword,cc_function"
BLOCK = re.compile(r"antimicrob|antibiotic|bacteriocin|defensin|cathelicidin|host.?defen[cs]e|innate.?immun|toxin|hemolysin|antibacterial", re.I)
AA = set("ACDEFGHIKLMNPQRSTVWY")


def fetch():
    params = urllib.parse.urlencode({"query": QUERY, "format": "tsv", "fields": FIELDS, "size": 500})
    url = "https://rest.uniprot.org/uniprotkb/search?" + params
    request = urllib.request.Request(url, headers={"User-Agent": "AMPProjectNegativeSet/1.0", "Accept": "text/tab-separated-values"})
    with urllib.request.urlopen(request, timeout=45) as response:
        return response.geturl(), list(csv.DictReader(response.read().decode("utf-8").splitlines(), delimiter="\t"))


def main():
    source_url, records = fetch()
    eligible, excluded = [], []
    for row in records:
        seq = (row.get("Sequence") or "").strip().upper()
        accession = row.get("Entry") or ""
        name, keywords, function = row.get("Protein names", ""), row.get("Keywords", ""), row.get("Function [CC]", "")
        if not (101 <= len(seq) <= 1000 and seq and set(seq) <= AA):
            excluded.append({"accession": accession, "reason": "parent sequence invalid or outside 101-1000 aa"}); continue
        if BLOCK.search(" ".join((name, keywords, function))):
            excluded.append({"accession": accession, "reason": "AMP/defense/toxin-related annotation or name"}); continue
        eligible.append({"sequence": seq, "accession": accession, "name": name,
                         "organism": row.get("Organism") or "Escherichia coli (taxon 562)", "keywords": keywords})

    panel_path = ROOT / "data" / "processed" / "sequence_provenance_panel.csv"
    with panel_path.open(encoding="utf-8-sig", newline="") as f:
        positives = [r["sequence"].strip().upper() for r in csv.DictReader(f)
                     if 5 <= len(r.get("sequence", "")) <= 100 and set(r["sequence"].strip().upper()) <= AA]
    if not eligible or not positives:
        raise SystemExit("UniProt returned no eligible parent proteins or local positive sequences.")

    rng, seen, accepted = random.Random(42), set(positives), []
    # One requested negative fragment for every local positive sequence, with
    # exactly the same length distribution to prevent length-only separation.
    for length in [len(s) for s in positives]:
        for _ in range(500):
            parent = rng.choice(eligible)
            protein = parent["sequence"]
            start = rng.randrange(0, len(protein) - length + 1)
            fragment = protein[start:start + length]
            if fragment in seen:
                continue
            seen.add(fragment)
            accepted.append({"sequence": fragment, "label": 0,
                             "source": "UniProtKB/Swiss-Prot reviewed parent; length-matched putative-negative fragment",
                             "source_id": parent["accession"], "protein_name": parent["name"],
                             "organism": parent["organism"], "uniprot_keywords": parent["keywords"],
                             "fragment_start_1based": start + 1, "fragment_end_1based": start + length,
                             "curation_note": "Random contiguous fragment from reviewed E. coli parent with no obvious AMP/antibiotic/defense/toxin annotation. Not experimentally tested; putative negative only."})
            break

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sequence", "label", "source", "source_id", "protein_name", "organism", "uniprot_keywords",
              "fragment_start_1based", "fragment_end_1based", "curation_note"]
    with OUT.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(accepted)
    manifest = {"retrieved_utc": datetime.now(timezone.utc).isoformat(), "source": "UniProtKB REST API",
                "query": QUERY, "fields": FIELDS, "source_url": source_url, "parent_records_returned": len(records),
                "eligible_parent_proteins": len(eligible), "positive_sequences_for_length_matching": len(positives),
                "negative_fragments_written": len(accepted), "excluded_parent_records": excluded,
                "fragment_rule": "Deterministic random seed 42; fragment lengths sampled from local positive sequences; fragment sequence exact matches to positives excluded.",
                "filter_rule": "Reviewed E. coli (taxon 562) parent proteins 101-1000 aa, standard amino acids; omit AMP/antibiotic/bacteriocin/defense/immunity/toxin/hemolysin-related names, keywords, or function text.",
                "label_semantics": "Weak/putative negative only. UniProt annotation absence and parental protein function do not prove that an isolated fragment lacks antimicrobial activity."}
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("retrieved_utc", "parent_records_returned", "eligible_parent_proteins", "positive_sequences_for_length_matching", "negative_fragments_written", "label_semantics")}, indent=2))


if __name__ == "__main__": main()
