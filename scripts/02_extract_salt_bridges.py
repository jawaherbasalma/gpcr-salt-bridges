#!/usr/bin/env python3
"""Find the salt bridges in each downloaded structure.

Reads  metadata/structures.csv          written by script 01
       metadata/receptor_regions.csv    written by script 01
       data/raw/<pdb_id>.cif            the structures

Writes data/processed/arpeggio/<pdb_id>.json   arpeggio's raw output
       data/processed/salt_bridges.csv         one row per residue pair

Run from the project root:

    pixi run python scripts/02_extract_salt_bridges.py

Arpeggio is slow — roughly a minute per structure. Its output is cached, so
rerunning only redoes structures whose JSON is missing.

What counts as a salt bridge is defined in docs/contact-spec.md. In short:
a residue pair where at least one atom-atom contact is labelled `ionic` by
arpeggio, both partners are amino acids, and both lie in the receptor rather
than in a fusion protein.
"""

import csv
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import gemmi

RAW_DIR = Path("data/raw")
META_DIR = Path("metadata")
PROCESSED_DIR = Path("data/processed")
ARPEGGIO_DIR = PROCESSED_DIR / "arpeggio"

OUTPUT = PROCESSED_DIR / "salt_bridges.csv"


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def run_arpeggio(pdb_id, chains):
    """Run arpeggio on the original structure, restricted to the receptor chains.

    Active-state structures are usually complexes -- receptor plus a G protein,
    often plus a nanobody -- and computing every contact in the whole assembly
    is slow for results the analysis then discards. Arpeggio's own `-s` flag
    narrows the computation to chosen chains.

    The file handed to arpeggio must be the one downloaded from the PDB. An
    earlier version of this script wrote a trimmed copy with gemmi instead,
    which was much faster but silently wrong: the rewritten file lost its
    polymer annotation, arpeggio then treated every residue as a free monomer
    rather than an amino acid, and stopped assigning charged groups entirely.
    Every structure came back with zero salt bridges.
    """
    output_path = ARPEGGIO_DIR / f"{pdb_id}.json"

    if not output_path.exists():
        structure_path = RAW_DIR / f"{pdb_id}.cif"
        if not structure_path.exists():
            raise SystemExit(f"Missing {structure_path}. Run script 01 first.")

        selection = [f"/{chain}//" for chain in sorted(chains)]
        print(f"  {pdb_id}: running arpeggio on chains {','.join(sorted(chains))}")
        result = subprocess.run(
            # The path comes first: `-s` takes any number of selections, so
            # argparse swallows a trailing filename as one more selection.
            ["pdbe-arpeggio", str(structure_path),
             "-o", str(ARPEGGIO_DIR), "-m",
             "-s", *selection],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 or not output_path.exists():
            raise SystemExit(
                f"arpeggio failed on {pdb_id}:\n{result.stderr[-2000:]}"
            )
    else:
        print(f"  {pdb_id}: using cached arpeggio output")

    with open(output_path) as handle:
        return json.load(handle)


def receptor_residues(pdb_id, regions):
    """Which residues belong to the receptor rather than a fusion protein.

    The regions table gives positions in the *entity* sequence, but arpeggio
    reports author residue numbers, and the two differ. gemmi exposes both for
    every residue, so this builds the set of author numbers that fall inside a
    receptor region.

    Returns {(chain, auth_seq_id), ...}.
    """
    structure = gemmi.read_structure(str(RAW_DIR / f"{pdb_id}.cif"))
    structure.setup_entities()

    spans = [(int(r["entity_seq_begin"]), int(r["entity_seq_end"])) for r in regions]
    chains_wanted = set()
    for region in regions:
        chains_wanted.update(region["chains"].split(";"))

    keep = set()
    for chain in structure[0]:
        if chain.name not in chains_wanted:
            continue
        for residue in chain:
            label_seq = residue.label_seq
            if label_seq is None:
                continue
            if any(start <= label_seq <= end for start, end in spans):
                keep.add((chain.name, residue.seqid.num))
    return keep


def salt_bridges(contacts, keep):
    """Collapse arpeggio's atom-atom contacts into residue-level salt bridges."""
    bridges = {}

    for contact in contacts:
        if "ionic" not in contact["contact"]:
            continue

        begin, end = contact["bgn"], contact["end"]

        # Amino acids only. Without this, glutamate-to-zinc is a salt bridge.
        if begin.get("label_comp_type") != "P" or end.get("label_comp_type") != "P":
            continue

        # Same chain only.
        if begin["auth_asym_id"] != end["auth_asym_id"]:
            continue
        chain = begin["auth_asym_id"]

        # Receptor only, not the fusion protein.
        if (chain, begin["auth_seq_id"]) not in keep:
            continue
        if (chain, end["auth_seq_id"]) not in keep:
            continue

        # Canonical order, so the lower residue number is always first.
        first, second = sorted(
            [
                (begin["auth_seq_id"], begin["label_comp_id"]),
                (end["auth_seq_id"], end["label_comp_id"]),
            ]
        )
        key = (chain, first, second)

        distance = contact["distance"]
        if key not in bridges or distance < bridges[key]:
            bridges[key] = distance

    return bridges


def main():
    structures = read_csv(META_DIR / "structures.csv")
    all_regions = read_csv(META_DIR / "receptor_regions.csv")

    regions_by_pdb = defaultdict(list)
    for region in all_regions:
        regions_by_pdb[region["pdb_id"]].append(region)

    ARPEGGIO_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for structure in structures:
        pdb_id = structure["pdb_id"]
        chains = set(filter(None, structure.get("receptor_chains", "").split(";")))

        contacts = run_arpeggio(pdb_id, chains)
        keep = receptor_residues(pdb_id, regions_by_pdb[pdb_id])
        bridges = salt_bridges(contacts, keep)

        for (chain, (num1, name1), (num2, name2)), distance in sorted(bridges.items()):
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "receptor": structure.get("receptor", ""),
                    "state": structure.get("state", ""),
                    "resolution_a": structure.get("resolution_a", ""),
                    "chain": chain,
                    "res1_num": num1,
                    "res1_name": name1,
                    "res2_num": num2,
                    "res2_name": name2,
                    "min_distance_a": distance,
                    "involves_histidine": name1 == "HIS" or name2 == "HIS",
                }
            )

        found_in = sorted({chain for chain, _, _ in bridges})
        print(
            f"  {pdb_id}: {len(bridges)} salt bridges "
            f"across {len(found_in)} chain(s), {len(keep)} receptor residues"
        )
        if not bridges:
            # Zero is possible but unusual. Say why, so a silent failure of
            # the kind described in run_arpeggio cannot pass unnoticed.
            ionic = sum(1 for c in contacts if "ionic" in c["contact"])
            polymer = sum(1 for c in contacts
                          if c["bgn"].get("label_comp_type") == "P")
            print(f"      no bridges: {len(contacts)} contacts, {ionic} ionic, "
                  f"{polymer} on polymer residues")

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n{len(rows)} salt bridges written to {OUTPUT}")


if __name__ == "__main__":
    sys.exit(main())
