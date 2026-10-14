#!/usr/bin/env python3
"""
14_check_arabidopsis_anchors.py

Check the Arabidopsis anchor identifiers in the curated candidate list against UniProt.

Reads data/genotypes/candidate_genes_curated.csv and, for every AGI in it (the current
`arabidopsis_agi` and the `arabidopsis_agi_as_given` before the 2026-10-14 correction), asks
UniProt (reviewed Swiss-Prot entries, Araport cross-reference) for the gene names that entry
carries. 03_candidate_genes.py accepts an anchor only when the intended Arabidopsis gene symbol is
one of those names, and reports what each replaced AGI actually named.

Writes data/genotypes/arabidopsis_anchor_check.csv. This needs network access, so it is not a
pipeline step: run it by hand and commit its output.

Usage:
    python code/14_check_arabidopsis_anchors.py

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
"""

from __future__ import annotations

import csv
import datetime
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

CURATED = Path("data/genotypes/candidate_genes_curated.csv")
OUT = Path("data/genotypes/arabidopsis_anchor_check.csv")
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
HEADERS = {"User-Agent": "brachypodium-gwas-spaceflight/anchor-check"}


def uniprot_lookup(agi: str) -> tuple[list[dict], str]:
    """Return reviewed UniProt entries cross-referenced to an Araport AGI, and the UniProt release."""
    query = f"xref:araport-{agi} AND organism_id:3702 AND reviewed:true"
    params = urllib.parse.urlencode({"query": query, "fields": "accession,gene_names,protein_name", "format": "tsv"})
    req = urllib.request.Request(f"{UNIPROT}?{params}", headers=HEADERS)
    with urllib.request.urlopen(req, timeout=60) as r:
        release = r.headers.get("X-UniProt-Release", "")
        lines = r.read().decode().strip().split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line], release


def main() -> None:
    curated = list(csv.DictReader(CURATED.open()))
    agis = {r["arabidopsis_agi"] for r in curated} | {r["arabidopsis_agi_as_given"] for r in curated}
    today = datetime.date.today().isoformat()
    rows = []
    for agi in sorted(agis):
        entries, release = uniprot_lookup(agi)
        if len(entries) != 1:
            sys.exit(f"{agi}: expected one reviewed UniProt entry, got {len(entries)}")
        e = entries[0]
        rows.append({
            "arabidopsis_agi": agi,
            "uniprot_accession": e["Entry"],
            "uniprot_gene_names": e["Gene Names"],
            "uniprot_protein_name": e["Protein names"],
            "uniprot_release": release,
            "checked_on": today,
        })
        print(f"{agi}\t{e['Entry']}\t{e['Gene Names']}")
        time.sleep(0.2)

    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"Saved {len(rows)} anchor checks to {OUT}")


if __name__ == "__main__":
    main()
