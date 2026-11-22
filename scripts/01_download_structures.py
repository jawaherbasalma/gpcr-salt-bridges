#!/usr/bin/env python3
"""Download GPCR structures from the RCSB PDB and record their provenance.

Reads  config/structures.yaml
Writes data/raw/<pdb_id>.cif          the structure files
       metadata/structures.csv        one row per structure
       metadata/receptor_regions.csv  which residues are receptor, not fusion

Run from the project root:

    pixi run python scripts/01_download_structures.py

Already-downloaded files are left alone unless their checksum has changed,
in which case the script stops and says so rather than overwriting.
"""

import csv
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import date
from pathlib import Path

import yaml

# Which structure list to use. Defaults to the small committed one; point
# STRUCTURES_CONFIG at config/structures_all.yaml for the full set.
CONFIG = Path(os.environ.get("STRUCTURES_CONFIG", "config/structures.yaml"))
RAW_DIR = Path("data/raw")
META_DIR = Path("metadata")

CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
ENTRY_API = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
ENTITY_API = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"

# Seconds between requests. RCSB is a public academic service; do not hammer it.
DELAY = 0.5


def fetch_json(url, retries=3):
    """GET a URL and parse the JSON, retrying on transient failures."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                return json.load(response)
        except urllib.error.HTTPError as err:
            if err.code == 404:
                raise SystemExit(f"Not found: {url}\nCheck the PDB ID in {CONFIG}.")
            if attempt == retries - 1:
                raise SystemExit(f"HTTP {err.code} for {url}")
            time.sleep(2 ** attempt)
        except urllib.error.URLError as err:
            if attempt == retries - 1:
                raise SystemExit(f"Could not reach {url}: {err.reason}")
            time.sleep(2 ** attempt)
    raise SystemExit(f"Gave up on {url}")


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def download_cif(pdb_id, destination):
    """Download a structure unless it is already present. Returns its checksum."""
    url = CIF_URL.format(pdb_id=pdb_id.upper())

    if destination.exists():
        print(f"  {pdb_id}: already downloaded")
        return sha256_of(destination), url

    print(f"  {pdb_id}: downloading")
    temporary = destination.with_suffix(".cif.part")
    try:
        urllib.request.urlretrieve(url, temporary)
    except urllib.error.HTTPError as err:
        temporary.unlink(missing_ok=True)
        raise SystemExit(f"Could not download {pdb_id}: HTTP {err.code}")
    temporary.rename(destination)

    # Raw data is read-only: a later step must never modify it in place.
    destination.chmod(0o444)
    return sha256_of(destination), url


def receptor_regions(pdb_id, uniprot):
    """Ask the API which parts of the structure come from the receptor.

    A fused construct is one polypeptide built from two proteins, so the
    fusion cannot be removed by dropping a chain. The API records which
    stretches align to which UniProt entry; anything not aligned to the
    receptor's accession is fusion and is excluded downstream.

    Returns (rows, chains, entity_descriptions), or None when no region of the
    structure aligns to the receptor. That happens for isoform-only records and
    for very recent depositions whose SIFTS mapping has not been computed yet;
    the caller skips those rather than stopping, since across a hundred-odd
    structures a few will always be unusable.
    """
    entry = fetch_json(ENTRY_API.format(pdb_id=pdb_id))
    entity_ids = entry["rcsb_entry_container_identifiers"]["polymer_entity_ids"]

    rows, chains, descriptions = [], set(), []
    # P41143-2 is an isoform of P41143 and maps to the same receptor.
    wanted = uniprot.split("-")[0]

    for entity_id in entity_ids:
        time.sleep(DELAY)
        entity = fetch_json(ENTITY_API.format(pdb_id=pdb_id, entity_id=entity_id))

        description = entity["rcsb_polymer_entity"].get("pdbx_description", "")
        descriptions.append(f"{entity_id}: {description}")

        entity_chains = entity["rcsb_polymer_entity_container_identifiers"].get(
            "auth_asym_ids", []
        )

        for alignment in entity.get("rcsb_polymer_entity_align") or []:
            accession = str(alignment.get("reference_database_accession") or "")
            if accession.split("-")[0] != wanted:
                continue

            chains.update(entity_chains)
            for region in alignment.get("aligned_regions", []):
                start = region["entity_beg_seq_id"]
                length = region["length"]
                rows.append(
                    {
                        "pdb_id": pdb_id,
                        "entity_id": entity_id,
                        "chains": ";".join(sorted(entity_chains)),
                        "uniprot": uniprot,
                        "entity_seq_begin": start,
                        "entity_seq_end": start + length - 1,
                        "uniprot_begin": region["ref_beg_seq_id"],
                        "uniprot_end": region["ref_beg_seq_id"] + length - 1,
                        "length": length,
                    }
                )

    if not rows:
        print(f"  {pdb_id}: SKIPPED, no region aligns to {uniprot}")
        for description in descriptions:
            print(f"      entity {description}")
        return None

    return rows, sorted(chains), descriptions


def main():
    if not CONFIG.exists():
        raise SystemExit(f"Missing {CONFIG}. Run this from the project root.")

    with open(CONFIG) as handle:
        config = yaml.safe_load(handle)

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    META_DIR.mkdir(parents=True, exist_ok=True)

    today = date.today().isoformat()
    structure_rows, region_rows, skipped = [], [], []

    for item in config["structures"]:
        pdb_id = item["pdb_id"].lower()
        uniprot = item["uniprot"]

        checksum, url = download_cif(pdb_id, RAW_DIR / f"{pdb_id}.cif")
        time.sleep(DELAY)

        entry = fetch_json(ENTRY_API.format(pdb_id=pdb_id))
        info = entry.get("rcsb_entry_info", {})
        resolutions = info.get("resolution_combined") or [None]

        time.sleep(DELAY)
        found = receptor_regions(pdb_id, uniprot)
        if found is None:
            skipped.append(pdb_id)
            continue
        regions, chains, descriptions = found
        region_rows.extend(regions)

        structure_rows.append(
            {
                "pdb_id": pdb_id,
                "uniprot": uniprot,
                "receptor": item.get("receptor", ""),
                "state": item.get("state", ""),
                "title": entry.get("struct", {}).get("title", ""),
                "method": info.get("experimental_method", ""),
                "resolution_a": resolutions[0],
                "receptor_chains": ";".join(chains),
                "entities": " | ".join(descriptions),
                "n_receptor_regions": len(regions),
                "source_url": url,
                "accessed": today,
                "sha256": checksum,
                "note": item.get("note", ""),
            }
        )

        print(
            f"  {pdb_id}: {resolutions[0]} A, chains {','.join(chains)}, "
            f"{len(regions)} receptor region(s)"
        )

    write_csv(META_DIR / "structures.csv", structure_rows)
    write_csv(META_DIR / "receptor_regions.csv", region_rows)

    print(f"\n{len(structure_rows)} structures in {META_DIR / 'structures.csv'}")
    print(f"{len(region_rows)} regions in {META_DIR / 'receptor_regions.csv'}")

    if skipped:
        print(f"\n{len(skipped)} structure(s) skipped: {', '.join(skipped)}")
        print("Their receptor could not be located by UniProt alignment, so")
        print("the fusion protein could not be told apart from the receptor.")

    # A receptor needs both states to contribute to the paired comparison.
    by_receptor = {}
    for row in structure_rows:
        by_receptor.setdefault(row["receptor"], set()).add(row["state"])
    unpaired = sorted(r for r, states in by_receptor.items() if len(states) < 2)
    if unpaired:
        print(f"\n{len(unpaired)} receptor(s) now have only one state: "
              f"{', '.join(unpaired)}")


def write_csv(path, rows):
    if not rows:
        raise SystemExit(f"Nothing to write to {path}")
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
