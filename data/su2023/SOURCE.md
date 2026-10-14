# Su et al. (2023) supplementary tables

`su2023_supplementary_tables.xlsx` is the file `Supplementary Tables.xlsx` from the supplementary archive of:

> Su S-H, Levine HG, Masson PH (2023). Brachypodium distachyon Seedlings Display Accession-Specific Morphological and
> Transcriptomic Responses to the Microgravity Environment of the International Space Station. *Life* 13(3):626.
> DOI [10.3390/life13030626](https://doi.org/10.3390/life13030626). PMID 36983782, PMCID PMC10058394.

**Licence:** the article is distributed under the Creative Commons Attribution 4.0 licence (CC BY 4.0,
<https://creativecommons.org/licenses/by/4.0/>), © 2023 by the authors, as stated in the article's PMC record.
The workbook is reproduced unchanged; the values derived from it in this repository are credited to Su et al. (2023).

**Where it came from (14 October 2026):** the Europe PMC supplementary-files service for PMC10058394
(`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10058394/supplementaryFiles`), which returns a zip containing
the article's `life-13-00626-s001.zip`; that archive contains `life-2159103-supplementary/Supplementary Tables.xlsx`.

| File | SHA-256 |
|---|---|
| `life-13-00626-s001.zip` | `9961849b0ced6c7d571ae2cd307ccbfca06836c9a48f5d237232134d64f4ecd6` |
| `Supplementary Tables.xlsx` (this file) | `8a26964bba458f22381e54ccb9d627587cfb9c11190cae9ec5de9155732ca872` |

`code/16_su2023_degs.py --download` fetches the archive again and checks both hashes.

**What is used:** Supplemental Tables S2 (Bd21), S3 (Bd21-3) and S4 (Gaz8), the per-organ lists of genes Su et al.
report as differentially expressed between ISS and ground control (DESeq log2FC and padj, edgeR log2FC and q, and
their best-hit Arabidopsis and rice orthologs). Their gene IDs are Phytozome v3.1 locus names (`Bradi3g01250`); the
Ensembl Plants equivalent is `BRADI_3g01250v3` (same JGI annotation of the Bd21 v3.0 assembly). See
`id_mapping_check.csv` for the check of that correspondence.

**Notes on the tables as published:**
- The tables list 1,352 Bd21 rows (325 root, 1,027 shoot), 569 Bd21-3 rows (224, 345) and 4,661 Gaz8 rows
  (475, 4,186); the article's abstract gives 1,353, 570 and 4,527 DEGs. The tables are used as published.
- Every listed row has DESeq padj < 0.05 and edgeR q < 0.05. The table captions call padj Bonferroni-corrected;
  the article's Methods say Benjamini–Hochberg. Values are reported as given.
- Su et al. excluded one Gaz8 root ground-control sample (low mapping rate) and its matched flight sample.
