"""
Download script for NASA OSDR study OSD-375 (Brachypodium distachyon APEX-06 spaceflight experiment).

Context:
This project analyzes data from the APEX-06 spaceflight experiment, which investigated
the transcriptomic response of three Brachypodium distachyon ecotypes (Bd21, Bd21-3, Gaz8)
to spaceflight in both root and shoot tissues.
Published paper: Su et al. (2023) Life 13(3):633, DOI: 10.3390/life13030633

Usage:
    python 01_download_osdr_osd375.py [--output-dir DIR] [--dry-run]
"""

import os
import csv
import json
import logging
import argparse
from urllib.parse import urljoin
from typing import List, Dict, Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm

OSDR_API_URL = "https://osdr.nasa.gov/osdr/data/osd/files/375"
OSDR_BASE_URL = "https://osdr.nasa.gov"

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def get_session_with_retries() -> requests.Session:
    """Create a requests Session with retry logic."""
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.mount("http://", HTTPAdapter(max_retries=retries))
    return session


def fetch_file_manifest() -> List[Dict[str, Any]]:
    """Fetch the file manifest from the OSDR API."""
    logger.info(f"Fetching OSDR file manifest from {OSDR_API_URL}")
    session = get_session_with_retries()
    response = session.get(OSDR_API_URL, timeout=30)
    response.raise_for_status()
    data = response.json()
    
    # OSDR API response parsing
    files = data.get("study", {}).get("375", {}).get("files", [])
    
    if not files:
        logger.warning("No files found in the API response under expected keys. Parsing tree...")
        def extract_files(node):
            extracted = []
            if isinstance(node, dict):
                if "file_name" in node and "remote_url" in node:
                    extracted.append(node)
                if "children" in node:
                    for child in node["children"]:
                        extracted.extend(extract_files(child))
            elif isinstance(node, list):
                for item in node:
                    extracted.extend(extract_files(item))
            return extracted
        files = extract_files(data)

    logger.info(f"Found {len(files)} files in the manifest.")
    return files


def download_file(url: str, dest_path: str, session: requests.Session) -> bool:
    """Download a file with a progress bar."""
    if os.path.exists(dest_path):
        logger.info(f"File already exists, skipping: {dest_path}")
        return True

    logger.info(f"Downloading {url} to {dest_path}")
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    
    try:
        response = session.get(url, stream=True, timeout=30)
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))
        
        with open(dest_path, "wb") as f, tqdm(
            desc=os.path.basename(dest_path),
            total=total_size,
            unit="iB",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                size = f.write(chunk)
                bar.update(size)
        return True
    except Exception as e:
        logger.error(f"Failed to download {url}: {e}")
        if os.path.exists(dest_path):
            os.remove(dest_path)
        return False


def main():
    parser = argparse.ArgumentParser(description="Download OSDR OSD-375 Data.")
    parser.add_argument("--output-dir", type=str, default="data/osdr", help="Output directory")
    parser.add_argument("--dry-run", action="store_true", help="Generate manifest without downloading")
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    manifest_csv = os.path.join(output_dir, "file_manifest.csv")
    
    files = fetch_file_manifest()
    
    session = get_session_with_retries()
    
    manifest_rows = []
    
    for file_info in files:
        file_name = file_info.get("file_name")
        remote_url = file_info.get("remote_url")
        file_size = file_info.get("file_size", 0)
        category = file_info.get("category", "Unknown")
        
        if not file_name or not remote_url:
            continue
            
        full_url = f"{OSDR_BASE_URL}{remote_url}"
        
        # Categorize
        cat_lower = category.lower()
        if "fastq" in file_name.lower():
            cat = "Raw Data"
        elif "metadata" in file_name.lower() or file_name.endswith(".txt") or file_name.endswith(".zip"):
            cat = "Study Metadata Files"
        elif "processed" in cat_lower or "counts" in file_name.lower():
            cat = "Processed Data"
        else:
            cat = category

        local_dir = os.path.join(output_dir, cat.replace(" ", "_"))
        local_path = os.path.join(local_dir, file_name)
        downloaded = False
        
        if cat == "Raw Data":
            logger.info(f"Skipping download for Raw Data: {file_name}. S3/Download URL: {full_url}")
            local_path = ""
        else:
            if not args.dry_run:
                downloaded = download_file(full_url, local_path, session)
            else:
                logger.info(f"[DRY-RUN] Would download: {file_name} -> {local_path}")
                
        manifest_rows.append({
            "file_name": file_name,
            "category": cat,
            "file_size_bytes": file_size,
            "download_url": full_url,
            "local_path": local_path,
            "downloaded": downloaded
        })
        
    with open(manifest_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["file_name", "category", "file_size_bytes", "download_url", "local_path", "downloaded"])
        writer.writeheader()
        writer.writerows(manifest_rows)
        
    logger.info(f"Manifest written to {manifest_csv}")


if __name__ == "__main__":
    main()
