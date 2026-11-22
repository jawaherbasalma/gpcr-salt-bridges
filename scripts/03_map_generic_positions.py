#!/usr/bin/env python3
"""Work out which residue in each structure is 3.50, 6.30 and so on.

Generic numbering (Ballesteros-Weinstein) names the same functional position
across every class A GPCR: 3.50 is the conserved arginine of the DRY motif
whether it is residue 131 or 135 or 139. GPCRdb publishes that mapping against
UniProt sequences. Structures, though, use author residue numbers, and those
are not UniProt numbers.

How the mapping is done, and why it is not arithmetic
-----------------------------------------------------
An earlier version of this script converted a UniProt position into an author
number by offset: the regions table says which stretch of the construct
sequence corresponds to which stretch of UniProt, so add the difference. That
is wrong for most crystallised GPCRs. Receptors are rarely solved as the
native protein -- intracellular loop 3 is commonly deleted to make the
receptor crystallise -- and a single aligned region cannot express an internal
deletion. In 4bvn, 28 residues of ICL3 are missing, so every position after
the loop came out 28 too high: position 6.30 was reported as residue 313,
which is elsewhere in the protein, and the distance to 3.50 came out at 34 A
when the two helices cannot be more than about 15 A apart. Neither RCSB nor
PDBe reports the gap -- both collapse the alignment into one straight line --
so there is no authoritative mapping to fetch.

This version places each position by sequence context instead. It takes the
eleven-residue window around the target in the UniProt sequence and finds
where that window sits in the structure, then reads the author number off the
middle residue. That survives internal deletions, N- and C-terminal
truncation, and fusion renumbering, because it never assumes the two sequences
run in step.

Two deliberate tolerances:

    MAX_MISMATCH  Engineered constructs carry point mutations -- 4bvn has
                  Y227A, a thermostabilising substitution right at position
                  5.58. A window is allowed a couple of mismatches so a
                  mutation does not prevent the position being found.

    MIN_INFORMATIVE  Unresolved residues leave gaps in the model. Only
                  positions that lie in a contiguous modelled stretch are
                  compared, and a position is placed only if enough of the
                  window is comparable.

A position is placed only when exactly one window matches best. If two match
equally well, or none matches well enough, the position is recorded as
unplaced rather than guessed at.

Because placement is by context, the amino acid check now means something
different and sharper. Under the old scheme a wrong residue could mean either
a mutation or a broken mapping, and there was no telling which. Here the
surrounding ten residues have already agreed, so a mismatch at the centre can
only be a real substitution in the construct.

Reads  config/structures.yaml           for the positions of interest
       metadata/structures.csv          written by script 01
       metadata/receptor_regions.csv    written by script 01, for chain names
       data/raw/<pdb_id>.cif

Writes data/raw/gpcrdb/<accession>.json   cached GPCRdb residue tables
       metadata/key_positions.csv         one row per structure per position

Run from the project root:

    pixi run python scripts/03_map_generic_positions.py
"""

import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

import gemmi
import yaml

# Which structure list to use. Defaults to the small committed one; point
# STRUCTURES_CONFIG at config/structures_all.yaml for the full set.
CONFIG = Path(os.environ.get("STRUCTURES_CONFIG", "config/structures.yaml"))
RAW_DIR = Path("data/raw")
GPCRDB_DIR = RAW_DIR / "gpcrdb"
META_DIR = Path("metadata")

PROTEIN_API = "https://gpcrdb.org/services/protein/accession/{accession}/"
RESIDUES_API = "https://gpcrdb.org/services/residues/{entry_name}/"

DELAY = 0.5

# Sequence-context anchoring. See the module docstring.
WINDOW = 5           # residues either side of the target, so an 11-mer
MIN_INFORMATIVE = 7  # how many of the 11 must be comparable
MAX_MISMATCH = 2     # engineered constructs carry point mutations


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

    GPCRdb writes generic numbers as '3.50x50' -- the Ballesteros-Weinstein
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


def uniprot_sequence(payload):
    """{position: one_letter} for the whole receptor, from GPCRdb."""
    return {r["sequence_number"]: r["amino_acid"] for r in payload["residues"]}


def one_letter(residue):
    """The one-letter code for an amino acid residue, or None if it is not one."""
    info = gemmi.find_tabulated_residue(residue.name)
    if info is None or not info.is_amino_acid():
        return None
    code = info.one_letter_code
    return code.upper() if code else None


def chain_residues(structure, chain_name):
    """[(author_number, one_letter), ...] for the modelled amino acids, in order."""
    out = []
    for chain in structure[0]:
        if chain.name != chain_name:
            continue
        for residue in chain:
            code = one_letter(residue)
            if code:
                out.append((residue.seqid.num, code))
    return out


def anchor(residues, uniprot_seq, position):
    """Find the author number for a UniProt position by sequence context.

    Returns (author_number, observed_aa, context_mismatches), or None when the
    position cannot be placed confidently.
    """
    offsets = range(-WINDOW, WINDOW + 1)
    expected = {offset: uniprot_seq.get(position + offset) for offset in offsets}

    candidates = []
    for index, (author_here, _) in enumerate(residues):
        informative = matches = 0
        for offset in offsets:
            other = index + offset
            if not 0 <= other < len(residues):
                continue
            author_there, code = residues[other]
            # Compare only across a contiguous modelled stretch. Where
            # residues are missing from the model the author numbering jumps,
            # and comparing across the jump would line the window up against
            # the wrong residues.
            if author_there - author_here != offset:
                continue
            wanted = expected[offset]
            if wanted is None:
                continue
            informative += 1
            matches += code == wanted
        if informative >= MIN_INFORMATIVE and informative - matches <= MAX_MISMATCH:
            candidates.append((matches, informative, index))

    if not candidates:
        return None

    best_score = max(candidate[0] for candidate in candidates)
    best = [candidate for candidate in candidates if candidate[0] == best_score]
    if len(best) != 1:
        # Two places in the sequence fit equally well. Refuse rather than pick.
        return None

    matches, informative, index = best[0]
    author_number, code = residues[index]
    return author_number, code, informative - matches


def main():
    with open(CONFIG) as handle:
        config = yaml.safe_load(handle)

    wanted = [str(p) for p in config.get("generic_positions",
                                         ["3.49", "3.50", "5.58", "6.30"])]
    print(f"Positions: {', '.join(wanted)}\n")

    structures = read_csv(META_DIR / "structures.csv")

    chains_by_pdb = defaultdict(set)
    for region in read_csv(META_DIR / "receptor_regions.csv"):
        chains_by_pdb[region["pdb_id"]].update(region["chains"].split(";"))

    # One GPCRdb lookup per receptor, not per structure.
    generic_by_accession = {}
    sequence_by_accession = {}
    for structure in structures:
        accession = structure["uniprot"]
        if accession not in generic_by_accession:
            payload = gpcrdb_residues(accession)
            generic_by_accession[accession] = generic_to_uniprot(payload, wanted)
            sequence_by_accession[accession] = uniprot_sequence(payload)
            time.sleep(DELAY)

    rows = []
    for structure in structures:
        pdb_id = structure["pdb_id"]
        accession = structure["uniprot"]
        generic = generic_by_accession[accession]
        uniprot_seq = sequence_by_accession[accession]

        model = gemmi.read_structure(str(RAW_DIR / f"{pdb_id}.cif"))
        model.setup_entities()
        chains = sorted(chains_by_pdb[pdb_id])
        residues_by_chain = {name: chain_residues(model, name) for name in chains}

        for position in wanted:
            if position not in generic:
                print(f"  {pdb_id}: GPCRdb has no {position} for {accession}")
                continue
            uniprot_position, expected_aa = generic[position]

            for chain in chains:
                row = {
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
                    "context_mismatches": "",
                }

                placed = anchor(residues_by_chain[chain], uniprot_seq, uniprot_position)
                if placed is not None:
                    author_number, observed_aa, context_mismatches = placed
                    row.update(
                        {
                            "auth_seq_id": author_number,
                            "observed_aa": observed_aa,
                            "resolved": True,
                            "aa_matches": observed_aa == expected_aa,
                            "context_mismatches": context_mismatches,
                        }
                    )
                rows.append(row)

    META_DIR.mkdir(parents=True, exist_ok=True)
    output = META_DIR / "key_positions.csv"
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    placed = [row for row in rows if row["resolved"]]
    unplaced = [row for row in rows if not row["resolved"]]
    substituted = [row for row in placed if not row["aa_matches"]]

    print(f"\n{len(rows)} rows written to {output}")
    print(f"{len(placed)} positions placed, {len(unplaced)} not placed")

    if unplaced:
        print(f"\n{len(unplaced)} position(s) could not be placed. Either the "
              f"residue is absent from the model, or no window matched well "
              f"enough to be sure:")
        for row in unplaced[:20]:
            print(f"  {row['pdb_id']} {row['chain']} {row['generic_number']} "
                  f"(UniProt {row['uniprot_position']}, expected "
                  f"{row['expected_aa']})")

    if substituted:
        print(f"\n{len(substituted)} position(s) carry a different amino acid "
              f"than the native sequence. The surrounding residues matched, so "
              f"these are substitutions in the construct rather than mapping "
              f"errors:")
        for row in substituted[:20]:
            print(f"  {row['pdb_id']} {row['chain']} {row['generic_number']}: "
                  f"{row['expected_aa']}{row['uniprot_position']}"
                  f"{row['observed_aa']} at author residue {row['auth_seq_id']}")
    else:
        print("\nEvery placed position carries the native amino acid.")


if __name__ == "__main__":
    sys.exit(main())
