#!/usr/bin/env python3
"""Test whether the ionic lock breaks when a GPCR activates.

Hypothesis: a salt bridge between generic positions 3.50 and 6.30 is present
in inactive-state structures and absent in active-state ones.

Two readouts per structure:

  presence   is the 3.50-6.30 pair in the salt bridge table?
  distance   the closest heavy-atom approach between those two residues,
             measured from coordinates whether or not it qualifies as a
             salt bridge. A pair at 4.2 A is not a bridge but is clearly
             different from one at 12 A, and the binary throws that away.

Reads  data/processed/salt_bridges.csv   written by script 02
       metadata/key_positions.csv        written by script 03
       data/raw/<pdb_id>.cif

Writes results/lock_status.csv           one row per structure per chain
       results/lock_test.json            the statistics

Run from the project root:

    pixi run python scripts/04_test_lock.py
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import gemmi
from scipy.stats import fisher_exact, mannwhitneyu

RAW_DIR = Path("data/raw")
META_DIR = Path("metadata")
PROCESSED_DIR = Path("data/processed")
RESULTS_DIR = Path("results")

POSITION_A = "3.50"
POSITION_B = "6.30"


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def min_heavy_atom_distance(pdb_id, chain_name, num_a, num_b):
    """Closest approach between two residues, ignoring hydrogens."""
    structure = gemmi.read_structure(str(RAW_DIR / f"{pdb_id}.cif"))
    structure.remove_hydrogens()

    residue_a = residue_b = None
    for chain in structure[0]:
        if chain.name != chain_name:
            continue
        for residue in chain:
            if residue.seqid.num == num_a:
                residue_a = residue
            elif residue.seqid.num == num_b:
                residue_b = residue

    if residue_a is None or residue_b is None:
        return None

    closest = None
    for atom_a in residue_a:
        for atom_b in residue_b:
            distance = atom_a.pos.dist(atom_b.pos)
            if closest is None or distance < closest:
                closest = distance
    return round(closest, 2)


def main():
    positions = read_csv(META_DIR / "key_positions.csv")
    bridges = read_csv(PROCESSED_DIR / "salt_bridges.csv")

    # Which residue is 3.50 and which is 6.30, per structure and chain.
    lookup = {}
    for row in positions:
        if row["generic_number"] in (POSITION_A, POSITION_B):
            key = (row["pdb_id"], row["chain"], row["generic_number"])
            lookup[key] = row

    # Every salt-bridged residue pair, for quick membership tests.
    bridged = set()
    for row in bridges:
        bridged.add(
            (row["pdb_id"], row["chain"], int(row["res1_num"]), int(row["res2_num"]))
        )

    structures = {}
    for row in positions:
        structures[row["pdb_id"]] = (row["receptor"], row["state"])

    chains_by_pdb = defaultdict(set)
    for row in positions:
        chains_by_pdb[row["pdb_id"]].add(row["chain"])

    results = []
    for pdb_id in sorted(structures):
        receptor, state = structures[pdb_id]

        for chain in sorted(chains_by_pdb[pdb_id]):
            row_a = lookup.get((pdb_id, chain, POSITION_A))
            row_b = lookup.get((pdb_id, chain, POSITION_B))

            if not row_a or not row_b:
                continue

            resolved = row_a["resolved"] == "True" and row_b["resolved"] == "True"

            if not resolved:
                results.append(
                    {
                        "pdb_id": pdb_id,
                        "receptor": receptor,
                        "state": state,
                        "chain": chain,
                        "res_3x50": row_a["auth_seq_id"],
                        "aa_3x50": row_a["observed_aa"],
                        "res_6x30": row_b["auth_seq_id"],
                        "aa_6x30": row_b["observed_aa"],
                        "both_resolved": False,
                        "lock_present": "",
                        "min_distance_a": "",
                    }
                )
                continue

            num_a = int(row_a["auth_seq_id"])
            num_b = int(row_b["auth_seq_id"])
            low, high = sorted([num_a, num_b])

            present = (pdb_id, chain, low, high) in bridged
            distance = min_heavy_atom_distance(pdb_id, chain, num_a, num_b)

            results.append(
                {
                    "pdb_id": pdb_id,
                    "receptor": receptor,
                    "state": state,
                    "chain": chain,
                    "res_3x50": num_a,
                    "aa_3x50": row_a["observed_aa"],
                    "res_6x30": num_b,
                    "aa_6x30": row_b["observed_aa"],
                    "both_resolved": True,
                    "lock_present": present,
                    "min_distance_a": distance,
                }
            )

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output = RESULTS_DIR / "lock_status.csv"
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    # One row per structure for the test, so two copies of the same receptor
    # in one crystal do not count as two independent observations.
    per_structure = {}
    for row in results:
        if row["both_resolved"] and row["pdb_id"] not in per_structure:
            per_structure[row["pdb_id"]] = row
    usable = list(per_structure.values())

    inactive = [r for r in usable if r["state"] == "inactive"]
    active = [r for r in usable if r["state"] == "active"]

    table = [
        [sum(1 for r in inactive if r["lock_present"]),
         sum(1 for r in inactive if not r["lock_present"])],
        [sum(1 for r in active if r["lock_present"]),
         sum(1 for r in active if not r["lock_present"])],
    ]

    odds_ratio, p_fisher = fisher_exact(table)

    distances_inactive = [r["min_distance_a"] for r in inactive if r["min_distance_a"]]
    distances_active = [r["min_distance_a"] for r in active if r["min_distance_a"]]

    p_mannwhitney = None
    if len(distances_inactive) >= 3 and len(distances_active) >= 3:
        _, p_mannwhitney = mannwhitneyu(
            distances_inactive, distances_active, alternative="two-sided"
        )

    stats = {
        "positions": [POSITION_A, POSITION_B],
        "structures_usable": len(usable),
        "inactive_with_lock": table[0][0],
        "inactive_total": len(inactive),
        "active_with_lock": table[1][0],
        "active_total": len(active),
        "fisher_odds_ratio": None if odds_ratio != odds_ratio else odds_ratio,
        "fisher_p": p_fisher,
        "median_distance_inactive_a": (
            sorted(distances_inactive)[len(distances_inactive) // 2]
            if distances_inactive else None
        ),
        "median_distance_active_a": (
            sorted(distances_active)[len(distances_active) // 2]
            if distances_active else None
        ),
        "mannwhitney_p": p_mannwhitney,
    }

    with open(RESULTS_DIR / "lock_test.json", "w") as handle:
        json.dump(stats, handle, indent=2)

    print(f"{POSITION_A} - {POSITION_B}\n")
    print(f"{'receptor':12s} {'state':9s} {'pdb':6s} {'3.50':>6s} {'6.30':>6s} "
          f"{'lock':>6s} {'dist':>6s}")
    for row in sorted(usable, key=lambda r: (r["receptor"], r["state"])):
        print(f"{row['receptor']:12s} {row['state']:9s} {row['pdb_id']:6s} "
              f"{row['aa_3x50']}{row['res_3x50']:<5} {row['aa_6x30']}{row['res_6x30']:<5} "
              f"{str(row['lock_present']):>6s} {row['min_distance_a']:>6}")

    print(f"\ninactive: {table[0][0]}/{len(inactive)} have the lock")
    print(f"active:   {table[1][0]}/{len(active)} have the lock")
    print(f"Fisher exact p = {p_fisher:.4f}")
    if p_mannwhitney is not None:
        print(f"Mann-Whitney on distances p = {p_mannwhitney:.4f}")

    skipped = [r for r in results if not r["both_resolved"]]
    if skipped:
        print(f"\n{len(skipped)} chain(s) excluded, a position not modelled:")
        for row in skipped:
            print(f"  {row['pdb_id']} chain {row['chain']}")


if __name__ == "__main__":
    sys.exit(main())
