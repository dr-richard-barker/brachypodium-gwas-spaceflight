#!/usr/bin/env python3
"""
03_candidate_genes.py

Build the gravitropism candidate-gene table from the curated list.

Inputs:
  - data/genotypes/candidate_genes_curated.csv: the hand-curated list. Each candidate is named
    after the Arabidopsis gravitropism gene it is meant to follow, with the AGI given for it.
  - data/genotypes/arabidopsis_anchor_check.csv: the UniProt names for each of those AGIs,
    written by 14_check_arabidopsis_anchors.py.

Output:
  - data/genotypes/gravitropism_candidate_genes.csv

An Arabidopsis anchor is kept only when UniProt lists the intended gene symbol for that AGI.
Otherwise the AGI is left blank and the row says which gene the given AGI actually is.
Brachypodium gene IDs are not assigned yet: the earlier BRADI IDs were withdrawn because most
did not exist in Ensembl Plants. They are TODO until orthologs are mapped (Ensembl Compara).

Usage:
    python code/03_candidate_genes.py

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

CURATED = Path("data/genotypes/candidate_genes_curated.csv")
CHECK = Path("data/genotypes/arabidopsis_anchor_check.csv")
OUT = Path("data/genotypes/gravitropism_candidate_genes.csv")


def main() -> None:
    curated = pd.read_csv(CURATED)
    check = pd.read_csv(CHECK).set_index("arabidopsis_agi")

    missing = sorted(set(curated["arabidopsis_agi_as_given"]) - set(check.index))
    if missing:
        sys.exit(f"No UniProt check for {missing}: run code/14_check_arabidopsis_anchors.py")

    rows = []
    for _, r in curated.iterrows():
        agi = r["arabidopsis_agi_as_given"]
        c = check.loc[agi]
        names = c["uniprot_gene_names"].split()
        if r["arabidopsis_gene"].upper() in (n.upper() for n in names):
            anchor = agi
            status = f"verified: UniProt {c['uniprot_accession']} ({' '.join(names[:3])})"
        else:
            anchor = ""
            actual = "/".join(n for n in names[:2] if not n.upper().startswith(agi[:2]))
            status = (f"TODO: the AGI given earlier, {agi}, is {actual} "
                      f"(UniProt {c['uniprot_accession']}), not {r['arabidopsis_gene']}")
        rows.append({
            "candidate": r["candidate"],
            "pathway": r["pathway"],
            "arabidopsis_gene": r["arabidopsis_gene"],
            "arabidopsis_agi": anchor,
            "anchor_status": status,
            "brachypodium_gene_id": "TODO",
        })

    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    n_ok = (out["arabidopsis_agi"] != "").sum()
    print(f"Saved {len(out)} candidates to {OUT}: {n_ok} anchors verified, {len(out) - n_ok} TODO")


if __name__ == "__main__":
    main()
