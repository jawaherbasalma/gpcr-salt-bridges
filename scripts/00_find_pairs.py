#!/usr/bin/env python3
"""Find every class A receptor solved in both an inactive and an active state.

Writes config/structures.yaml, so the structure list is computed from the
GPCRdb catalogue rather than copied out of a paper. The selection rule is:

    every class A receptor with at least one Inactive and at least one Active
    structure; the best-resolution structure in each state.

Intermediate-state structures are ignored -- the comparison is between the two
endpoints.

Usage, from the project root:

    pixi run python scripts/00_find_pairs.py            # all pairs
    pixi run python scripts/00_find_pairs.py --limit 12 # best 12 by resolution

Each structure costs about a minute of arpeggio time later, so --limit is
worth using while the pipeline is still being built.
"""

import argparse
import json
import sys
import time
import urllib.request
from collections import defaultdict
from datetime import date
from pathlib import Path

STRUCTURE_API = "https://gpcrdb.org/services/structure/"
PROTEIN_API = "https://gpcrdb.org/services/protein/{entry_name}/"

CONFIG = Path("config/structures.yaml")
CACHE = Path("data/raw/gpcrdb/structure_catalogue.json")

CLASS_A = "Class A (Rhodopsin)"
DELAY = 0.3


def fetch_json(url, timeout=60):
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return json.load(response)


def catalogue():
    """The full GPCRdb structure list, cached so reruns work offline."""
    if CACHE.exists():
        print(f"Using cached catalogue: {CACHE}")
        with open(CACHE) as handle:
            return json.load(handle)

    print("Fetching the GPCRdb structure catalogue")
    data = fetch_json(STRUCTURE_API)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE, "w") as handle:
        json.dump(data, handle)
    return data


def best_per_state(structures):
    """{'inactive': entry, 'active': entry} keeping the best resolution each."""
    best = {}
    for entry in structures:
        state = str(entry.get("state", "")).lower()
        if state not in ("inactive", "active"):
            continue
        resolution = entry.get("resolution")
        if resolution is None:
            continue
        current = best.get(state)
        if current is None or resolution < current["resolution"]:
            best[state] = entry
    return best


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None,
                        help="keep only the N pairs with the best resolution")
    parser.add_argument("--output", default=str(CONFIG),
                        help="where to write the config file")
    args = parser.parse_args()
    output = Path(args.output)

    data = catalogue()
    print(f"{len(data)} structures in the catalogue")

    by_receptor = defaultdict(list)
    for entry in data:
        if entry.get("class") != CLASS_A:
            continue
        by_receptor[entry["protein"]].append(entry)

    print(f"{len(by_receptor)} class A receptors with at least one structure")

    pairs = {}
    for receptor, structures in by_receptor.items():
        best = best_per_state(structures)
        if "inactive" in best and "active" in best:
            pairs[receptor] = best

    print(f"{len(pairs)} receptors have both an inactive and an active structure")

    # Rank by the worse of the two resolutions, so a pair is only as good as
    # its weaker half.
    ranked = sorted(
        pairs.items(),
        key=lambda item: max(item[1]["inactive"]["resolution"],
                             item[1]["active"]["resolution"]),
    )
    if args.limit:
        ranked = ranked[: args.limit]
        print(f"keeping the best {len(ranked)}")

    # UniProt accessions, one lookup per receptor.
    print("\nLooking up UniProt accessions")
    entries = []
    for receptor, best in ranked:
        protein = fetch_json(PROTEIN_API.format(entry_name=receptor))
        accession = protein.get("accession")
        if not accession:
            print(f"  {receptor}: no accession, skipped")
            continue
        time.sleep(DELAY)

        short = receptor.split("_")[0]
        for state in ("inactive", "active"):
            entry = best[state]
            entries.append(
                {
                    "pdb_id": entry["pdb_code"],
                    "uniprot": accession,
                    "receptor": short,
                    "state": state,
                    "note": (
                        f"{protein.get('name', receptor)}. "
                        f"{entry.get('type', '')}, {entry['resolution']} A, "
                        f"preferred chain {entry.get('preferred_chain', '?')}."
                    ),
                }
            )
        print(f"  {short:10s} {accession}  "
              f"{best['inactive']['pdb_code']} / {best['active']['pdb_code']}")

    lines = [
        "# Structures to analyse -- GENERATED, do not edit by hand.",
        "#",
        "# Written by scripts/00_find_pairs.py from the GPCRdb structure",
        "# catalogue (https://gpcrdb.org/services/structure/).",
        "#",
        "# Selection rule: every class A receptor with at least one Inactive",
        "# and at least one Active structure; the best-resolution structure in",
        "# each state. Intermediate states are ignored.",
        "#",
        f"# Generated {date.today().isoformat()}: {len(entries) // 2} pairs,",
        f"# {len(entries)} structures.",
        "",
        "generic_positions:",
        '  - "3.49"',
        '  - "3.50"',
        '  - "5.58"',
        '  - "6.30"',
        "",
        "structures:",
    ]
    for entry in entries:
        lines += [
            f"  - pdb_id: {entry['pdb_id']}",
            f"    uniprot: {entry['uniprot']}",
            f"    receptor: {entry['receptor']}",
            f"    state: {entry['state']}",
            f"    note: {entry['note']}",
            "",
        ]

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines))
    print(f"\n{len(entries)} structures ({len(entries) // 2} pairs) -> {output}")
    print("Previous contents are in git if you want them back.")


if __name__ == "__main__":
    sys.exit(main())
