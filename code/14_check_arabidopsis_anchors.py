#!/usr/bin/env python3
"""
14_check_arabidopsis_anchors.py

Check the identifiers in the curated candidate list.

Reads data/genotypes/candidate_genes_curated.csv and:
  1. for every AGI given there, asks UniProt (reviewed Swiss-Prot entries, Araport
     cross-reference) for the gene names that entry carries. An anchor counts as verified only
     when the intended Arabidopsis gene symbol is one of those names.
     Writes data/genotypes/arabidopsis_anchor_check.csv.
  2. looks up each withdrawn BRADI identifier (the Brachypodium IDs the list gave before
     2026-10-14) in Ensembl, to record which of them exist, and, where the Arabidopsis anchor
     was verified in step 1, whether the withdrawn ID is among the anchor's Brachypodium
     orthologs in Ensembl Compara (plants).
     Writes data/genotypes/withdrawn_bradi_id_check.csv.

This needs network access, so it is not a pipeline step: run it by hand and commit its outputs.
03_candidate_genes.py reads the anchor check.

Usage:
    python code/14_check_arabidopsis_anchors.py

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
"""

from __future__ import annotations

import csv
import datetime
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

CURATED = Path("data/genotypes/candidate_genes_curated.csv")
OUT = Path("data/genotypes/arabidopsis_anchor_check.csv")
OUT_BRADI = Path("data/genotypes/withdrawn_bradi_id_check.csv")
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
ENSEMBL = "https://rest.ensembl.org"
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


def ensembl_get(path: str, tries: int = 5) -> dict | None:
    """GET an Ensembl REST path as JSON; None when Ensembl answers 400 (unknown identifier).

    Ensembl REST returns intermittent 429/5xx errors; those are retried, never read as "not found".
    """
    req = urllib.request.Request(f"{ENSEMBL}{path}", headers={**HEADERS, "Content-Type": "application/json"})
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 400:
                return None
            if e.code in (429, 500, 502, 503, 504) and attempt < tries - 1:
                time.sleep(2 * (attempt + 1))
                continue
            raise


def compara_orthologs(agi: str) -> list[str]:
    """Brachypodium orthologs of an Arabidopsis gene in Ensembl Compara (plants)."""
    d = ensembl_get(f"/homology/id/arabidopsis_thaliana/{agi}?compara=plants;"
                    "target_species=brachypodium_distachyon;type=orthologues;format=condensed")
    if not d or not d.get("data"):
        return []
    return sorted(h["id"] for h in d["data"][0]["homologies"])


def check_withdrawn_bradi_ids(curated: list[dict], verified: dict[str, bool], today: str) -> None:
    release = ensembl_get("/info/eg_version")["version"]
    rows = []
    for r in curated:
        gid = r["withdrawn_bradi_id"]
        hit = ensembl_get(f"/lookup/id/{gid}")
        agi = r["arabidopsis_agi_as_given"]
        if verified[r["candidate"]]:
            orthologs = compara_orthologs(agi)
            is_ortholog = "yes" if gid in orthologs else "no"
        else:
            orthologs, is_ortholog = [], "not checked (anchor not verified)"
        rows.append({
            "candidate": r["candidate"],
            "withdrawn_bradi_id": gid,
            "found_in_ensembl": "yes" if hit else "no",
            "anchor_agi": agi if verified[r["candidate"]] else "",
            "compara_brachypodium_orthologs": " ".join(orthologs),
            "withdrawn_id_is_compara_ortholog": is_ortholog,
            "ensembl_genomes_release": release,
            "checked_on": today,
        })
        time.sleep(0.2)
    with OUT_BRADI.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    n_found = sum(r["found_in_ensembl"] == "yes" for r in rows)
    n_orth = sum(r["withdrawn_id_is_compara_ortholog"] == "yes" for r in rows)
    n_checked = sum(r["withdrawn_id_is_compara_ortholog"] in ("yes", "no") for r in rows)
    print(f"Saved {len(rows)} withdrawn BRADI ID lookups to {OUT_BRADI}: {n_found} found in Ensembl, "
          f"{len(rows) - n_found} not found; {n_orth} of {n_checked} with a verified anchor are a Compara ortholog of it")


def main() -> None:
    curated = list(csv.DictReader(CURATED.open()))
    today = datetime.date.today().isoformat()
    rows = []
    for agi in sorted({r["arabidopsis_agi_as_given"] for r in curated}):
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

    names = {r["arabidopsis_agi"]: r["uniprot_gene_names"].upper().split() for r in rows}
    verified = {r["candidate"]: r["arabidopsis_gene"].upper() in names[r["arabidopsis_agi_as_given"]]
                for r in curated}
    check_withdrawn_bradi_ids(curated, verified, today)


if __name__ == "__main__":
    main()
