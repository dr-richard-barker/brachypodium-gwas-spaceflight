#!/usr/bin/env python3
"""
15_map_orthologs.py

Map each verified Arabidopsis anchor to its orthologs in Brachypodium distachyon, rice
(Oryza sativa) and bread wheat (Triticum aestivum) with Ensembl Compara (plants), and check the
withdrawn Brachypodium identifiers against those orthologs.

Inputs:
  - data/genotypes/candidate_genes_curated.csv
  - data/genotypes/arabidopsis_anchor_check.csv (from 14_check_arabidopsis_anchors.py)

Outputs:
  - data/genotypes/compara_orthologs.csv: one row per anchor x ortholog, with the Compara homology
    type (one2one, one2many, many2many) and percent identities. Anchors whose UniProt check fails
    are skipped.
  - data/genotypes/withdrawn_bradi_id_check.csv: for each BRADI identifier the list gave before
    2026-10-14, whether it exists in Ensembl and whether it is among the anchor's Compara orthologs.

This needs network access, so it is not a pipeline step: run it by hand and commit its outputs.
Ensembl REST returns intermittent 429/5xx errors; those are retried, never read as "not found".

Usage:
    python code/15_map_orthologs.py

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
"""

from __future__ import annotations

import csv
import datetime
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

CURATED = Path("data/genotypes/candidate_genes_curated.csv")
ANCHOR_CHECK = Path("data/genotypes/arabidopsis_anchor_check.csv")
OUT_ORTHOLOGS = Path("data/genotypes/compara_orthologs.csv")
OUT_WITHDRAWN = Path("data/genotypes/withdrawn_bradi_id_check.csv")
ENSEMBL = "https://rest.ensembl.org"
HEADERS = {"User-Agent": "brachypodium-gwas-spaceflight/ortholog-map", "Content-Type": "application/json"}
TARGETS = ["brachypodium_distachyon", "oryza_sativa", "triticum_aestivum"]


def ensembl_get(path: str, tries: int = 5) -> dict | None:
    """GET an Ensembl REST path as JSON; None when Ensembl answers 400 (unknown identifier)."""
    req = urllib.request.Request(f"{ENSEMBL}{path}", headers=HEADERS)
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


def compara_orthologs(agi: str) -> list[dict]:
    """Orthologs of an Arabidopsis gene in the target species, from Ensembl Compara (plants)."""
    targets = ";".join(f"target_species={t}" for t in TARGETS)
    d = ensembl_get(f"/homology/id/arabidopsis_thaliana/{agi}?compara=plants;type=orthologues;{targets};sequence=none")
    if not d or not d.get("data"):
        return []
    return [{
        "target_species": h["target"]["species"],
        "target_gene_id": h["target"]["id"],
        "homology_type": h["type"],
        "target_perc_id": round(h["target"]["perc_id"], 1),
        "anchor_perc_id": round(h["source"]["perc_id"], 1),
    } for h in d["data"][0]["homologies"]]


def main() -> None:
    curated = list(csv.DictReader(CURATED.open()))
    names = {r["arabidopsis_agi"]: r["uniprot_gene_names"].upper().split() for r in csv.DictReader(ANCHOR_CHECK.open())}
    verified = {r["arabidopsis_agi"] for r in curated if r["arabidopsis_gene"].upper() in names[r["arabidopsis_agi"]]}
    release = ensembl_get("/info/eg_version")["version"]
    today = datetime.date.today().isoformat()

    orthologs: dict[str, list[dict]] = {}
    rows = []
    for agi in sorted(verified):
        orthologs[agi] = compara_orthologs(agi)
        for o in sorted(orthologs[agi], key=lambda o: (TARGETS.index(o["target_species"]), o["target_gene_id"])):
            rows.append({"arabidopsis_agi": agi, **o, "ensembl_genomes_release": release, "checked_on": today})
        n_bd = sum(o["target_species"] == TARGETS[0] for o in orthologs[agi])
        print(f"{agi}: {len(orthologs[agi])} orthologs ({n_bd} Brachypodium)")
        time.sleep(0.3)
    with OUT_ORTHOLOGS.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"Saved {len(rows)} ortholog rows for {len(verified)} anchors to {OUT_ORTHOLOGS}")

    checks = []
    for r in curated:
        gid = r["withdrawn_bradi_id"]
        agi = r["arabidopsis_agi"]
        bd = sorted(o["target_gene_id"] for o in orthologs.get(agi, []) if o["target_species"] == TARGETS[0])
        checks.append({
            "candidate": r["candidate"],
            "withdrawn_bradi_id": gid,
            "found_in_ensembl": "yes" if ensembl_get(f"/lookup/id/{gid}") else "no",
            "anchor_agi": agi if agi in verified else "",
            "compara_brachypodium_orthologs": " ".join(bd),
            "withdrawn_id_is_compara_ortholog": ("yes" if gid in bd else "no") if agi in verified
                                                else "not checked (anchor not verified)",
            "ensembl_genomes_release": release,
            "checked_on": today,
        })
        time.sleep(0.2)
    with OUT_WITHDRAWN.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(checks[0]))
        w.writeheader()
        w.writerows(checks)
    n_found = sum(c["found_in_ensembl"] == "yes" for c in checks)
    n_orth = sum(c["withdrawn_id_is_compara_ortholog"] == "yes" for c in checks)
    n_checked = sum(c["withdrawn_id_is_compara_ortholog"] in ("yes", "no") for c in checks)
    print(f"Saved {len(checks)} withdrawn-ID checks to {OUT_WITHDRAWN}: {n_found} found in Ensembl; "
          f"{n_orth} of {n_checked} with a verified anchor are a Compara ortholog of it")


if __name__ == "__main__":
    main()
