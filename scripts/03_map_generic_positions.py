#!/usr/bin/env python3
"""Work out which residue in each structure is 3.50, 6.30 and so on.

Generic numbering (Ballesteros-Weinstein) names the same functional position
across every class A GPCR: 3.50 is the conserved arginine of the DRY motif
whether it is residue 131 or 135 or 139. GPCRdb publishes that mapping against
UniProt sequences.

Structures, though, use author residue numbers, which are not always UniProt
numbers. So this script chains two mappings:

    generic number -> UniProt position     (from GPCRdb)
    UniProt position -> author number      (from the regions table + the file)

Getting this wrong shifts every position silently, so the output records the
amino acid found at each mapped position. If 3.50 is not an arginine,
something is wrong.

Reads  config/structures.yaml           for the positions of interest
       metadata/structures.csv          written by script 01
       metadata/receptor_regions.csv    written by script 01
       data/raw/<pdb_id>.cif

Writes data/raw/gpcrdb/<entry_name>.json   cached GPCRdb residue tables
       metadata/key_positions.csv          one row per structure per position

Run from the project root:

    pixi run python scripts/03_map_generic_positions.py
"""

import csv
import json
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

import gemmi
import yaml

CONFIG = Path("config/structures.yaml")
RAW_DIR = Path("data/raw")
GPCRDB_DIR = RAW_DIR / "gpcrdb"
META_DIR = Path("metadata")

PROTEIN_API = "https://gpcrdb.org/services/protein/accession/{accession}/"
RESIDUES_API = "https://gpcrdb.org/services/residues/{entry_name}/"

DELAY = 0.5


def fetch_json(url):
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            return json.load(response)
    except urllib.error.HTTPError as err:
        raise SystemExit(f"HTTP {err.code} for {url}")
    except urllib.error.URLError as err:
        raise SystemExit(f"Could not reach {url}: {err.reason}")


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def gpcrdb_residues(accession):
    """Fetch a receptor's residue table, caching it so reruns work offline."""
    cache = GPCRDB_DIR / f"{accession}.json"

    if cache.exists():
        with open(cache) as handle:
            return json.load(handle)

    print(f"  {accession}: fetching from GPCRdb")
    protein = fetch_json(PROTEIN_API.format(accession=accession))
    entry_name = protein.get("entry_name")
    if not entry_name:
        raise SystemExit(f"GPCRdb has no entry for {accession}")

    time.sleep(DELAY)
    residues = fetch_json(RESIDUES_API.format(entry_name=entry_name))

    payload = {"accession": accession, "entry_name": entry_name, "residues": residues}
    GPCRDB_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache, "w") as handle:
        json.dump(payload, handle, indent=1)
    return payload


def generic_to_uniprot(payload, wanted):
    """{'3.50': (131, 'R'), '6.30': (268, 'E'), ...} for the positions wanted.

    GPCRdb writes generic numbers as '3.50x50' — the Ballesteros-Weinstein
    number, then GPCRdb's own structure-based number after the x. Only the
    part before the x is used here.
    """
    found = {}
    for residue in payload["residues"]:
        generic = residue.get("display_generic_number")
        if not generic:
            continue
        bw = str(generic).split("x")[0]
        if bw in wanted:
            found[bw] = (residue["sequence_number"], residue["amino_acid"])
    return found


def uniprot_to_author(pdb_id, regions):
    """{(chain, uniprot_position): author_residue_number}.

    Two steps. The regions table says which stretch of the entity sequence
    corresponds to which stretch of UniProt, so a UniProt position becomes an
    entity position. The structure file then gives the author number for that
    entity position, which is what arpeggio reports.
    """
    structure = gemmi.read_structure(str(RAW_DIR / f"{pdb_id}.cif"))
    structure.setup_entities()

    # entity position -> UniProt position, for receptor regions only
    entity_to_uniprot = {}
    chains_wanted = set()
    for region in regions:
        chains_wanted.update(region["chains"].split(";"))
        entity_begin = int(region["entity_seq_begin"])
        uniprot_begin = int(region["uniprot_begin"])
        length = int(region["length"])
        for offset in range(length):
            entity_to_uniprot[entity_begin + offset] = uniprot_begin + offset

    mapping = {}
    for chain in structure[0]:
        if chain.name not in chains_wanted:
            continue
        for residue in chain:
            label_seq = residue.label_seq
            if label_seq is None:
                continue
            uniprot_position = entity_to_uniprot.get(label_seq)
            if uniprot_position is not None:
                mapping[(chain.name, uniprot_position)] = (
                    residue.seqid.num,
                    gemmi.find_tabulated_residue(residue.name).one_letter_code.upper()
                    if gemmi.find_tabulated_residue(residue.name)
                    else "?",
                )
    return mapping


def main():
    with open(CONFIG) as handle:
        config = yaml.safe_load(handle)

    wanted = [str(p) for p in config.get("generic_positions", ["3.49", "3.50", "5.58", "6.30"])]
    print(f"Positions: {', '.join(wanted)}\n")

    structures = read_csv(META_DIR / "structures.csv")

    regions_by_pdb = defaultdict(list)
    for region in read_csv(META_DIR / "receptor_regions.csv"):
        regions_by_pdb[region["pdb_id"]].append(region)

    # One GPCRdb lookup per receptor, not per structure.
    generic_by_accession = {}
    for structure in structures:
        accession = structure["uniprot"]
        if accession not in generic_by_accession:
            payload = gpcrdb_residues(accession)
            generic_by_accession[accession] = generic_to_uniprot(payload, wanted)
            time.sleep(DELAY)

    rows = []
    for structure in structures:
        pdb_id = structure["pdb_id"]
        accession = structure["uniprot"]
        generic = generic_by_accession[accession]

        author_numbers = uniprot_to_author(pdb_id, regions_by_pdb[pdb_id])
        chains = sorted({chain for chain, _ in author_numbers})

        for position in wanted:
            if position not in generic:
                print(f"  {pdb_id}: GPCRdb has no {position} for {accession}")
                continue
            uniprot_position, expected_aa = generic[position]

            for chain in chains:
                entry = author_numbers.get((chain, uniprot_position))
                if entry is None:
                    # The residue exists in the sequence but not in the model:
                    # unresolved density, or trimmed from the construct.
                    rows.append(
                        {
                            "pdb_id": pdb_id,
                            "receptor": structure.get("receptor", ""),
                            "state": structure.get("state", ""),
                            "chain": chain,
                            "generic_number": position,
                            "uniprot_position": uniprot_position,
                            "expected_aa": expected_aa,
                            "auth_seq_id": "",
                            "observed_aa": "",
                            "resolved": False,
                            "aa_matches": "",
                        }
                    )
                    continue

                author_number, observed_aa = entry
                rows.append(
                    {
                        "pdb_id": pdb_id,
                        "receptor": structure.get("receptor", ""),
                        "state": structure.get("state", ""),
                        "chain": chain,
                        "generic_number": position,
                        "uniprot_position": uniprot_position,
                        "expected_aa": expected_aa,
                        "auth_seq_id": author_number,
                        "observed_aa": observed_aa,
                        "resolved": True,
                        "aa_matches": observed_aa == expected_aa,
                    }
                )

    META_DIR.mkdir(parents=True, exist_ok=True)
    output = META_DIR / "key_positions.csv"
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    resolved = [r for r in rows if r["resolved"]]
    mismatched = [r for r in resolved if not r["aa_matches"]]

    print(f"\n{len(rows)} rows written to {output}")
    print(f"{len(resolved)} positions resolved in the model, "
          f"{len(rows) - len(resolved)} missing")

    if mismatched:
        print(f"\n{len(mismatched)} MISMATCHES — the mapping is wrong for these:")
        for row in mismatched[:20]:
            print(f"  {row['pdb_id']} {row['chain']} {row['generic_number']}: "
                  f"expected {row['expected_aa']}, found {row['observed_aa']} "
                  f"at {row['auth_seq_id']}")
    else:
        print("\nEvery resolved position has the amino acid GPCRdb predicts.")


if __name__ == "__main__":
    sys.exit(main())
