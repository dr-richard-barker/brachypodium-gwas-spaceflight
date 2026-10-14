#!/usr/bin/env python3
"""
03_candidate_genes.py

Build the gravitropism candidate-gene table from the curated list, the anchor check and the
Compara ortholog map.

Inputs:
  - data/genotypes/candidate_genes_curated.csv: the hand-curated list. Each candidate is named
    after the Arabidopsis gravitropism gene it follows (`arabidopsis_gene`), with its AGI
    (`arabidopsis_agi`) and the AGI given before the 2026-10-14 correction (`arabidopsis_agi_as_given`).
  - data/genotypes/arabidopsis_anchor_check.csv: UniProt gene names for each AGI
    (14_check_arabidopsis_anchors.py).
  - data/genotypes/compara_orthologs.csv: Brachypodium, rice and wheat orthologs of each anchor
    from Ensembl Compara (15_map_orthologs.py).

Output:
  - data/genotypes/gravitropism_candidate_genes.csv

An anchor is kept only when UniProt lists the intended gene symbol for its AGI. Orthologs are
listed as Compara gives them; one-to-many and many-to-many cases are kept and labelled, and a
Brachypodium gene shared by several candidates is flagged.

Usage:
    python code/03_candidate_genes.py

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

CURATED = Path("data/genotypes/candidate_genes_curated.csv")
CHECK = Path("data/genotypes/arabidopsis_anchor_check.csv")
ORTHOLOGS = Path("data/genotypes/compara_orthologs.csv")
OUT = Path("data/genotypes/gravitropism_candidate_genes.csv")


def gene_names(check: pd.DataFrame, agi: str) -> list[str]:
    return check.loc[agi, "uniprot_gene_names"].split()


def short_name(names: list[str], agi: str) -> str:
    """First one or two UniProt gene names, skipping the AGI-style locus name."""
    return "/".join(n for n in names[:2] if not n.upper().startswith(agi[:2]))


def main() -> None:
    curated = pd.read_csv(CURATED)
    check = pd.read_csv(CHECK).set_index("arabidopsis_agi")
    orthologs = pd.read_csv(ORTHOLOGS)

    needed = set(curated["arabidopsis_agi"]) | set(curated["arabidopsis_agi_as_given"])
    missing = sorted(needed - set(check.index))
    if missing:
        sys.exit(f"No UniProt check for {missing}: run code/14_check_arabidopsis_anchors.py")

    by_anchor: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for _, o in orthologs.iterrows():
        by_anchor[o["arabidopsis_agi"]][o["target_species"]].append(o)

    rows = []
    for _, r in curated.iterrows():
        agi, given = r["arabidopsis_agi"], r["arabidopsis_agi_as_given"]
        names = gene_names(check, agi)
        verified = r["arabidopsis_gene"].upper() in (n.upper() for n in names)
        if verified:
            status = f"verified: UniProt {check.loc[agi, 'uniprot_accession']} ({' '.join(names[:3])})"
            if given != agi:
                status += (f"; corrected on 2026-10-14 from {given}, which is "
                           f"{short_name(gene_names(check, given), given)}")
        else:
            status = (f"TODO: {agi} is {short_name(names, agi)} (UniProt {check.loc[agi, 'uniprot_accession']}), "
                      f"not {r['arabidopsis_gene']}")

        orth = by_anchor.get(agi, {}) if verified else {}
        bd = orth.get("brachypodium_distachyon", [])
        rice = orth.get("oryza_sativa", [])
        wheat = orth.get("triticum_aestivum", [])
        rows.append({
            "candidate": r["candidate"],
            "pathway": r["pathway"],
            "arabidopsis_gene": r["arabidopsis_gene"],
            "arabidopsis_agi": agi if verified else "",
            "anchor_status": status,
            "brachypodium_orthologs": " ".join(o["target_gene_id"] for o in bd),
            "brachypodium_homology": " ".join(sorted({o["homology_type"].replace("ortholog_", "") for o in bd})),
            "rice_orthologs": " ".join(o["target_gene_id"] for o in rice),
            "wheat_orthologs": " ".join(o["target_gene_id"] for o in wheat),
            "ortholog_note": "" if bd else ("TODO: no Brachypodium ortholog in Ensembl Compara (plants)"
                                            if verified else "TODO: anchor not verified"),
        })

    # Flag Brachypodium genes that are orthologs of more than one candidate's anchor.
    shared: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        for gid in row["brachypodium_orthologs"].split():
            shared[gid].append(row["candidate"])
    for row in rows:
        others = sorted({c for gid in row["brachypodium_orthologs"].split() for c in shared[gid]} - {row["candidate"]})
        if others:
            row["ortholog_note"] = "Brachypodium ortholog(s) shared with " + ", ".join(others)

    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    n_verified = (out["arabidopsis_agi"] != "").sum()
    n_bd = (out["brachypodium_orthologs"] != "").sum()
    print(f"Saved {len(out)} candidates to {OUT}: {n_verified} anchors verified, "
          f"{n_bd} with a Brachypodium ortholog in Compara")


if __name__ == "__main__":
    main()
