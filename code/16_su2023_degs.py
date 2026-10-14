#!/usr/bin/env python3
"""
16_su2023_degs.py

Extract the spaceflight DEG lists published by Su et al. (2023, Life 13:626, CC BY 4.0) and look
up the gravitropism candidate genes in them.

Inputs:
  - data/su2023/su2023_supplementary_tables.xlsx (Su et al. 2023 Supplementary Tables; see SOURCE.md)
  - data/genotypes/candidate_genes_curated.csv
  - data/genotypes/compara_orthologs.csv (from 15_map_orthologs.py)

Outputs:
  - data/su2023/su2023_degs.csv: every row of Supplemental Tables S2 (Bd21), S3 (Bd21-3) and S4
    (Gaz8), with the Ensembl ID for each Phytozome v3.1 gene ID (Bradi3g01250 -> BRADI_3g01250v3)
  - tables/su2023_deg_counts.csv: DEGs per accession, organ and direction
  - tables/candidate_spaceflight_degs.csv: candidate genes found in those lists, either as a Compara
    ortholog of the candidate's Arabidopsis anchor or as a gene whose best-hit Arabidopsis ortholog
    in Su et al.'s own annotation is the anchor

No new differential-expression analysis is done: every value is Su et al.'s.

Options (network; not used by the pipeline):
  --download   fetch the supplement again from Europe PMC and check its SHA-256 hashes
  --check-ids  check that every Su et al. gene ID exists in Ensembl under its converted name, and
               compare Su et al.'s best-hit Arabidopsis orthologs with Ensembl Compara for a seeded
               random sample; writes data/su2023/id_mapping_check.csv

Usage:
    python code/16_su2023_degs.py [--download] [--check-ids]

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import random
import re
import sys
import time
import urllib.error
import urllib.request
import zipfile
from collections import Counter
from pathlib import Path

import openpyxl

XLSX = Path("data/su2023/su2023_supplementary_tables.xlsx")
CURATED = Path("data/genotypes/candidate_genes_curated.csv")
ORTHOLOGS = Path("data/genotypes/compara_orthologs.csv")
OUT_DEGS = Path("data/su2023/su2023_degs.csv")
OUT_IDCHECK = Path("data/su2023/id_mapping_check.csv")
OUT_COUNTS = Path("tables/su2023_deg_counts.csv")
OUT_CANDIDATES = Path("tables/candidate_spaceflight_degs.csv")

SHEETS = {"Bd21": "Supplemental Table S2", "Bd21-3": "Supplemenetal Table S3", "Gaz8": "Supplemental Table S4"}
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10058394/supplementaryFiles"
ZIP_SHA256 = "9961849b0ced6c7d571ae2cd307ccbfca06836c9a48f5d237232134d64f4ecd6"
XLSX_SHA256 = "8a26964bba458f22381e54ccb9d627587cfb9c11190cae9ec5de9155732ca872"
ENSEMBL = "https://rest.ensembl.org"
HEADERS = {"User-Agent": "brachypodium-gwas-spaceflight/su2023", "Content-Type": "application/json",
           "Accept": "application/json"}
SU_ID = re.compile(r"Bradi(\d)g(\d{5})")


def to_ensembl(su_id: str) -> str:
    """Phytozome v3.1 locus name -> Ensembl Plants gene ID (Bradi3g01250 -> BRADI_3g01250v3)."""
    chrom, num = SU_ID.fullmatch(su_id).groups()
    return f"BRADI_{chrom}g{num}v3"


def download() -> None:
    req = urllib.request.Request(EPMC, headers={"User-Agent": HEADERS["User-Agent"]})
    with urllib.request.urlopen(req, timeout=120) as r:
        outer = zipfile.ZipFile(io.BytesIO(r.read()))
    inner_bytes = outer.read("life-13-00626-s001.zip")
    if hashlib.sha256(inner_bytes).hexdigest() != ZIP_SHA256:
        sys.exit("life-13-00626-s001.zip: SHA-256 differs from SOURCE.md")
    xlsx_bytes = zipfile.ZipFile(io.BytesIO(inner_bytes)).read("life-2159103-supplementary/Supplementary Tables.xlsx")
    if hashlib.sha256(xlsx_bytes).hexdigest() != XLSX_SHA256:
        sys.exit("Supplementary Tables.xlsx: SHA-256 differs from SOURCE.md")
    XLSX.write_bytes(xlsx_bytes)
    print(f"Downloaded and verified {XLSX}")


def parse_tables() -> list[dict]:
    if hashlib.sha256(XLSX.read_bytes()).hexdigest() != XLSX_SHA256:
        sys.exit(f"{XLSX}: SHA-256 differs from SOURCE.md")
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    rows = []
    for accession, sheet in SHEETS.items():
        organ = None
        for r in wb[sheet].iter_rows(values_only=True):
            if r[0] in ("Root", "Shoot"):
                organ = r[0].lower()
                continue
            gene = str(r[1]).strip() if r[1] is not None else ""
            if organ is None or not SU_ID.fullmatch(gene):
                continue
            best_at = str(r[6]).strip() if r[6] not in (None, 0) else ""
            best_os = str(r[8]).strip() if r[8] not in (None, 0) else ""
            rows.append({
                "accession": accession,
                "organ": organ,
                "su_gene_id": gene,
                "ensembl_gene_id": to_ensembl(gene),
                "deseq_log2fc": r[2],
                "deseq_padj": r[3],
                "edger_log2fc": r[4],
                "edger_q": r[5],
                "su_best_hit_arabidopsis": best_at,
                "su_best_hit_rice": best_os,
            })
    return rows


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def deg_counts(degs: list[dict]) -> list[dict]:
    out = []
    for accession in SHEETS:
        for organ in ("root", "shoot"):
            sub = [d for d in degs if d["accession"] == accession and d["organ"] == organ]
            out.append({
                "accession": accession,
                "organ": organ,
                "up": sum(d["deseq_log2fc"] > 0 for d in sub),
                "down": sum(d["deseq_log2fc"] < 0 for d in sub),
                "total": len(sub),
            })
    return out


def candidate_hits(degs: list[dict]) -> list[dict]:
    curated = list(csv.DictReader(CURATED.open()))
    bd_orthologs: dict[str, dict[str, str]] = {}
    for o in csv.DictReader(ORTHOLOGS.open()):
        if o["target_species"] == "brachypodium_distachyon":
            bd_orthologs.setdefault(o["arabidopsis_agi"], {})[o["target_gene_id"]] = o["homology_type"]
    hits = []
    for c in curated:
        agi = c["arabidopsis_agi"]
        compara = bd_orthologs.get(agi, {})
        for d in degs:
            if d["ensembl_gene_id"] in compara:
                match, homology = "Compara ortholog of the anchor", compara[d["ensembl_gene_id"]]
            elif d["su_best_hit_arabidopsis"].upper() == agi:
                match, homology = "Su et al. best-hit Arabidopsis ortholog is the anchor", ""
            else:
                continue
            hits.append({
                "candidate": c["candidate"],
                "arabidopsis_gene": c["arabidopsis_gene"],
                "arabidopsis_agi": agi,
                "brachypodium_gene_id": d["ensembl_gene_id"],
                "match": match,
                "compara_homology_type": homology,
                "accession": d["accession"],
                "organ": d["organ"],
                "deseq_log2fc": d["deseq_log2fc"],
                "deseq_padj": d["deseq_padj"],
                "edger_log2fc": d["edger_log2fc"],
                "edger_q": d["edger_q"],
                "source": "Su et al. 2023, Life 13:626, Supplemental Tables S2-S4 (CC BY 4.0)",
            })
    return hits


def ensembl(path: str, body: dict | None = None, tries: int = 5):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(f"{ENSEMBL}{path}", data=data, headers=HEADERS)
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 400:
                return None
            if e.code in (429, 500, 502, 503, 504) and attempt < tries - 1:
                time.sleep(3 * (attempt + 1))
                continue
            raise


def check_ids(degs: list[dict]) -> None:
    ids = sorted({d["ensembl_gene_id"] for d in degs})
    found = {}
    for i in range(0, len(ids), 900):
        found.update(ensembl("/lookup/id", {"ids": ids[i:i + 900]}))
        time.sleep(1)
    n_found = sum(1 for v in found.values() if v)

    with_hit = sorted({(d["ensembl_gene_id"], d["su_best_hit_arabidopsis"].upper()) for d in degs
                       if d["su_best_hit_arabidopsis"].upper().startswith("AT")})
    sample = random.Random(42).sample(with_hit, 100)
    compara_at = []
    for gid, _ in sample:
        h = ensembl(f"/homology/id/brachypodium_distachyon/{gid}?compara=plants;type=orthologues;"
                    "target_species=arabidopsis_thaliana;format=condensed")
        compara_at.append({x["id"].upper() for x in h["data"][0]["homologies"]} if h and h.get("data") else set())
        time.sleep(0.2)
    agree = sum(at in orth for (_, at), orth in zip(sample, compara_at))
    # Null: the same comparison with Su et al.'s best hits shuffled between the sampled genes.
    shuffled = [at for _, at in sample]
    random.Random(43).shuffle(shuffled)
    agree_null = sum(at in orth for at, orth in zip(shuffled, compara_at))

    release = ensembl("/info/eg_version")["version"]
    rows = [
        {"check": "Su et al. gene IDs found in Ensembl under the converted name", "n": len(ids),
         "passed": n_found, "ensembl_genomes_release": release},
        {"check": "Su et al. best-hit Arabidopsis ortholog is among the gene's Ensembl Compara orthologs "
                  "(random sample of genes with a best hit, seed 42)", "n": len(sample),
         "passed": agree, "ensembl_genomes_release": release},
        {"check": "Same comparison with the best hits shuffled between the sampled genes (null, seed 43)",
         "n": len(sample), "passed": agree_null, "ensembl_genomes_release": release},
    ]
    write_csv(OUT_IDCHECK, rows, list(rows[0]))
    for r in rows:
        print(f"{r['check']}: {r['passed']} of {r['n']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    parser.add_argument("--download", action="store_true", help="re-download the supplement (network)")
    parser.add_argument("--check-ids", action="store_true", help="check the ID mapping against Ensembl (network)")
    args = parser.parse_args()

    if args.download:
        download()
    degs = parse_tables()
    write_csv(OUT_DEGS, degs, list(degs[0]))
    counts = deg_counts(degs)
    write_csv(OUT_COUNTS, counts, list(counts[0]))
    for c in counts:
        print(f"{c['accession']:7s} {c['organ']:5s} up {c['up']:5d} down {c['down']:5d} total {c['total']:5d}")
    not_sig = sum(not (d["deseq_padj"] < 0.05 and d["edger_q"] < 0.05) for d in degs)
    print(f"{len(degs)} DEG rows; {not_sig} with DESeq padj or edgeR q >= 0.05")

    hits = candidate_hits(degs)
    fields = ["candidate", "arabidopsis_gene", "arabidopsis_agi", "brachypodium_gene_id", "match",
              "compara_homology_type", "accession", "organ", "deseq_log2fc", "deseq_padj", "edger_log2fc",
              "edger_q", "source"]
    write_csv(OUT_CANDIDATES, hits, fields)
    by_match = Counter(h["match"] for h in hits)
    print(f"Saved {len(hits)} candidate DEG rows to {OUT_CANDIDATES}: {dict(by_match)}")

    if args.check_ids:
        check_ids(degs)


if __name__ == "__main__":
    main()
