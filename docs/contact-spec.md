# Contact and salt-bridge specification

This document defines exactly what counts as a residue contact and a salt bridge
in this project.

It has one job: it is the brief handed to each code-generation session. Two
implementations are built independently from this text and nothing else. Where
they disagree, the disagreement is traced back to an ambiguity here, and this
document is corrected.

Anything left vague here becomes a difference between the two implementations.
That is tolerable — it is part of the result — but the differences should be
interesting ones, not accidents of wording.

Parameter values below are defaults. The values actually used live in
`config/parameters.yaml`.

---

## 1. Scope

Experimentally determined structures of class A GPCRs, downloaded from the RCSB
PDB in mmCIF format. The structures used are listed in
`config/structures.yaml`.

Predicted structures (AlphaFold and similar) are out of scope.

---

## 2. Structure preparation

**Chain.** Only the receptor chain is analysed. Its chain identifier is given
per structure in `config/structures.yaml`, because it varies and cannot be
inferred reliably. Other chains — G protein subunits, nanobodies, antibody
fragments — are ignored.

**Assembly.** The asymmetric unit is used, not the biological assembly.

**Non-protein content.** Waters, ions, lipids, detergents and ligands are
excluded. Only amino acid residues are considered.

**Fusion proteins.** Many GPCR structures replace intracellular loop 3 with a
soluble fusion partner such as BRIL or T4 lysozyme. These are excluded.
Residues are treated as fusion when their sequence number falls outside the
receptor's native numbering range, which is given per structure in
`config/structures.yaml`. Structures where this cannot be determined cleanly
are dropped from the set and the reason recorded.

**Missing residues.** Unresolved residues are simply absent. No interpolation,
no renumbering. Gaps in numbering are expected and are not an error.

**Alternate conformations.** Where a residue is modelled in more than one
position, the conformation with the highest occupancy is used. On a tie, the
alphabetically first altloc identifier is used.

**Modified residues.** Selenomethionine (MSE) is treated as methionine. Other
modified residues are excluded, and each exclusion is recorded in the output
log.

**Hydrogens.** Ignored. Most crystal structures do not resolve them, so all
distances are between heavy atoms.

---

## 3. Contact definition

Two residues are in contact when any heavy atom of one lies within
**4.5 Å** of any heavy atom of the other.

**Sequence separation.** Residues closer than **4 positions apart** in sequence
are excluded, since neighbouring residues in a helix are trivially in contact.
Residues on either side of a gap caused by unresolved density are compared by
their sequence numbers, not by their position in the resolved chain.

**Chains.** Intra-chain only, following section 2.

---

## 4. Salt bridge definition

A salt bridge is a contact between oppositely charged side-chain groups.

**Negatively charged atoms**

| Residue | Atoms |
|---|---|
| Asp | OD1, OD2 |
| Glu | OE1, OE2 |

**Positively charged atoms**

| Residue | Atoms |
|---|---|
| Lys | NZ |
| Arg | NE, NH1, NH2 |

**Criterion.** A salt bridge exists when any nitrogen from the positive list is
within **4.0 Å** of any oxygen from the negative list.

**Histidine.** Excluded by default. Its protonation state depends on local pH
and is not determined by the structure. A configuration flag
(`include_histidine`) allows ND1 and NE2 to be added to the positive list; when
enabled, histidine bridges are marked in a separate column rather than mixed in.

**Termini.** The charged N- and C-termini of the chain are not counted. Only
side chains.

Every salt bridge is also a contact under section 3, so salt bridges are a
labelled subset of the contact table, not a separate calculation.

---

## 5. Output format

One CSV per structure per implementation, written to `data/processed/`.

| Column | Meaning |
|---|---|
| `pdb_id` | Four-character PDB identifier, lowercase |
| `chain` | Chain identifier as in the mmCIF file |
| `res1_num` | Residue number of the first residue |
| `res1_name` | Three-letter residue name |
| `res1_atom` | Atom name of the closest approach |
| `res2_num` | Residue number of the second residue |
| `res2_name` | Three-letter residue name |
| `res2_atom` | Atom name of the closest approach |
| `min_distance` | Closest heavy-atom distance, Ångström, 2 decimal places |
| `is_salt_bridge` | `true` or `false` per section 4 |

**Canonical ordering.** Within a row, the residue with the lower sequence number
is always `res1`. Without this rule the two implementations would emit the same
contact in opposite orders and the comparison would report false disagreements.

**Row ordering.** Sorted by `pdb_id`, then `res1_num`, then `res2_num`.

**Insertion codes.** Where a residue carries an insertion code, it is appended
to the residue number with no separator, for example `52A`.

---

## 6. How the two implementations are compared

Two rows describe the same contact when `pdb_id`, `chain`, `res1_num` and
`res2_num` are all equal.

`min_distance` may differ in the last decimal place without counting as a
disagreement; a tolerance of 0.05 Å is allowed. `res1_atom` and `res2_atom` may
differ when two atom pairs are nearly equidistant, and this is reported
separately rather than counted as a missed contact.

Reported per structure:

- contacts found by both implementations
- contacts found only by A
- contacts found only by B
- rows where `is_salt_bridge` disagrees

---

## 7. Decisions still open

Resolve these before generating code, or accept that the two implementations
will diverge here.

1. **Sequence separation.** Four positions is a guess. Too small and helical
   i,i+4 pairs dominate the output; too large and genuine short-range bridges
   are lost. The helical salt bridges may actually be worth keeping.
2. **Contact cutoff.** 4.5 Å between heavy atoms is a common choice. Some work
   uses 8 Å between Cβ atoms instead, which is coarser and would give a very
   different table.
3. **Histidine.** Excluding it is the conservative choice but loses real
   bridges.
4. **Fusion detection.** Relying on numbering ranges requires knowing each
   receptor's native range. An alternative is to read the entity description
   from the mmCIF file, which is more automatic but less predictable.

---

## 8. Sanity checks

These must hold, and failing one means the implementation is wrong rather than
the biology being surprising.

- **Rhodopsin, inactive state.** R135–E247 is present as a salt bridge.
- **β2-adrenergic receptor, inactive (2RH1).** R131 and E268 are the lock
  partners; check whether the bridge is present, since this structure is known
  not to show it clearly.
- **β2-adrenergic receptor, active (3SN6).** The R131–E268 bridge is absent.
- Every salt bridge also appears as a contact.
- No row has `res1_num` greater than or equal to `res2_num`.
