"""
Parse script for NASA OSDR study OSD-375 metadata.

Context:
This project analyzes data from the APEX-06 spaceflight experiment, which investigated
the transcriptomic response of three Brachypodium distachyon ecotypes (Bd21, Bd21-3, Gaz8)
to spaceflight in both root and shoot tissues.
Published paper: Su et al. (2023) Life 13(3):633, DOI: 10.3390/life13030633

Usage:
    python 02_parse_osdr_metadata.py [--input-dir DIR]
"""

import os
import glob
import zipfile
import argparse
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def unzip_metadata(input_dir: str):
    """Unzip any metadata zip files found in the input directory."""
    zip_files = glob.glob(os.path.join(input_dir, "**", "*.zip"), recursive=True)
        
    for zip_file in zip_files:
        extract_dir = os.path.dirname(zip_file)
        logger.info(f"Extracting {zip_file} to {extract_dir}")
        try:
            with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
        except Exception as e:
            logger.error(f"Failed to extract {zip_file}: {e}")

def parse_metadata(input_dir: str):
    """Parse sample and assay metadata files."""
    # Find the sample file
    sample_file_pattern = os.path.join(input_dir, "**", "s_*.txt")
    sample_files = glob.glob(sample_file_pattern, recursive=True)
    
    if not sample_files:
        logger.error("No sample metadata files (s_*.txt) found.")
        return
        
    sample_file = sample_files[0]
    logger.info(f"Parsing sample metadata: {sample_file}")
    
    try:
        s_df = pd.read_csv(sample_file, sep="\t")
    except Exception as e:
        logger.error(f"Failed to read {sample_file}: {e}")
        return
        
    # Attempt to locate key columns conceptually
    col_mapping = {}
    for col in s_df.columns:
        col_lower = col.lower()
        if "sample name" in col_lower:
            col_mapping["Sample_Name"] = col
        elif "organism" in col_lower and "part" not in col_lower:
            col_mapping["Organism"] = col
        elif "ecotype" in col_lower or "accession" in col_lower or "strain" in col_lower:
            col_mapping["Ecotype"] = col
        elif "tissue" in col_lower or "organism part" in col_lower:
            col_mapping["Tissue"] = col
        elif "spaceflight" in col_lower or "condition" in col_lower or "factor" in col_lower:
            if "Condition" not in col_mapping:
                col_mapping["Condition"] = col
        elif "replicate" in col_lower:
            col_mapping["Replicate"] = col
            
    extracted_s_df = pd.DataFrame()
    for nice_name, orig_name in col_mapping.items():
        extracted_s_df[nice_name] = s_df[orig_name]
        
    # Find the assay file
    assay_file_pattern = os.path.join(input_dir, "**", "a_*.txt")
    assay_files = glob.glob(assay_file_pattern, recursive=True)
    
    a_df = None
    if assay_files:
        assay_file = assay_files[0]
        logger.info(f"Parsing assay metadata: {assay_file}")
        try:
            a_df = pd.read_csv(assay_file, sep="\t")
        except Exception as e:
            logger.error(f"Failed to read {assay_file}: {e}")
            
    # Output files
    out_sample = os.path.join(input_dir, "sample_metadata.csv")
    out_summary = os.path.join(input_dir, "experimental_design_summary.csv")
    out_seq = os.path.join(input_dir, "sequencing_parameters.csv")
    
    # Save the cleaned mapping version as sample_metadata.csv
    if not extracted_s_df.empty:
        extracted_s_df.to_csv(out_sample, index=False)
        logger.info(f"Saved sample metadata to {out_sample}")
        
        # Generate summary
        summary_cols = [c for c in ["Ecotype", "Tissue", "Condition"] if c in extracted_s_df.columns]
        if summary_cols:
            summary = extracted_s_df.groupby(summary_cols).size().reset_index(name="Count")
            summary.to_csv(out_summary, index=False)
            logger.info(f"Saved experimental design summary to {out_summary}")
            print("\n--- Experimental Design Summary ---")
            print(summary.to_string(index=False))
            print("-----------------------------------\n")
    else:
        s_df.to_csv(out_sample, index=False)
        logger.info(f"Saved full sample metadata to {out_sample} (could not extract clean columns)")
            
    if a_df is not None:
        seq_cols = [c for c in a_df.columns if any(x in c.lower() for x in ["library", "instrument", "read", "platform"])]
        if seq_cols:
            seq_df = a_df[["Sample Name"] + seq_cols if "Sample Name" in a_df.columns else seq_cols]
            seq_df = seq_df.drop_duplicates()
            seq_df.to_csv(out_seq, index=False)
            logger.info(f"Saved sequencing parameters to {out_seq}")
        else:
            a_df.to_csv(out_seq, index=False)
            logger.info(f"Saved full assay metadata as sequencing parameters to {out_seq}")


def main():
    parser = argparse.ArgumentParser(description="Parse OSDR OSD-375 Metadata.")
    parser.add_argument("--input-dir", type=str, default="data/osdr", help="Directory containing downloaded OSDR data")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    
    unzip_metadata(input_dir)
    parse_metadata(input_dir)


if __name__ == "__main__":
    main()
