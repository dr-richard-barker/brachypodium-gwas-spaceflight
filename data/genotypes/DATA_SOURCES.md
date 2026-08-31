# Brachypodium distachyon Genotype Data Sources

This document details the data sources for accessing *Brachypodium distachyon* SNP and variation data for GWAS analysis.

## 1. Ensembl Plants
- **Release Info**: Current release (check https://plants.ensembl.org)
- **VCF Source URL**: `ftp://ftp.ensemblgenomes.org/pub/plants/current/vcf/brachypodium_distachyon/`
- **API**: REST API at `https://rest.ensembl.org` is used to query specific variants for target candidate genes.

## 2. BrachyPan Pan-Genome
- **Reference**: Gordon, S.P., et al. (2017). "Extensive gene content variation in the Brachypodium distachyon pan-genome correlates with population structure." *Nature Communications* 8:2042.
- **Access**: Available through JGI Phytozome and specific BrachyPan data portals.
- **Usage**: Critical for understanding structural variants between ecotypes like Bd21 and Gaz8.

## 3. JGI Phytozome
- **URL**: https://phytozome-next.jgi.doe.gov/
- **Dataset**: *Brachypodium distachyon* v3.1
- **Access**: Requires JGI SSO login. Users can download VCF files for diverse panels, including APEX-06 accessions.

## 4. OSD-375 (APEX-06) Metadata
- **URL**: https://osdr.nasa.gov/osdr/data/osd/files/375
- **Reference**: Su et al. (2023). "Spaceflight and Simulated Microgravity Induce Changes in Brachypodium distachyon Gene Expression." *Life* 13(3):633. DOI: 10.3390/life13030633

## Instructions for full VCF download
For full genome-wide association studies, download the comprehensive VCF files directly via FTP or globus from Ensembl or JGI:
```bash
wget ftp://ftp.ensemblgenomes.org/pub/plants/current/vcf/brachypodium_distachyon/brachypodium_distachyon.vcf.gz
```
