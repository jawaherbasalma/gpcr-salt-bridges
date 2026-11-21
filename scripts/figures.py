"""Figures for the ionic lock analysis.

Each function takes a dataframe and returns a matplotlib figure, so the Quarto
report can call them at render time rather than embedding pre-made images.

Colours: two categories only, inactive and active. Blue and orange from a
palette validated for colour-vision deficiency (worst-pair CVD dE 24.7, normal
33.6, both well clear of the floors). State is never carried by colour alone --
the x axis, the legend and the direct labels all repeat it.

To see them standalone:

    pixi run python scripts/figures.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# --- palette -------------------------------------------------------------
INACTIVE = "#2a78d6"   # categorical slot 1, blue
ACTIVE = "#eb6834"     # categorical slot 2, orange

INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"
BASELINE = "#c3c2b7"

STATE_COLOUR = {"inactive": INACTIVE, "active": ACTIVE}
STATE_ORDER = ["inactive", "active"]

RESULTS = Path("results")
PROCESSED = Path("data/processed")
META = Path("metadata")


def _style(ax):
    """Recessive chrome: hairline grid, no box, muted ticks."""
    ax.set_facecolor(SURFACE)
    ax.figure.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_MUTED, labelsize=9, length=3, width=0.8)
    ax.yaxis.label.set_color(INK_SECONDARY)
    ax.xaxis.label.set_color(INK_SECONDARY)
    ax.title.set_color(INK)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    return ax


def load():
    """The per-structure results, one row per structure, both positions resolved."""
    lock = pd.read_csv(RESULTS / "lock_status.csv")
    lock = lock[lock["both_resolved"]].copy()
    # One chain per structure, matching the statistical test.
    lock = lock.sort_values(["pdb_id", "chain"]).drop_duplicates("pdb_id")
    # Resolution is recorded in the structure metadata, not the lock table.
    meta = pd.read_csv(META / "structures.csv")[["pdb_id", "resolution_a"]]
    return lock.merge(meta, on="pdb_id", how="left")


def fig_distance_by_state(lock):
    """Every structure as a point. The separation is the finding, so show it.

    A bar of group means would hide that there are six and five observations
    and that they do not overlap.
    """
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    _style(ax)

    # Sorted by distance, then alternating left and right of the group centre.
    # Points close together vertically land on opposite sides, so their labels
    # cannot collide.
    for index, state in enumerate(STATE_ORDER):
        subset = (
            lock[lock["state"] == state]
            .sort_values("min_distance_a")
            .reset_index(drop=True)
        )

        for position, row in subset.iterrows():
            side = -1 if position % 2 == 0 else 1
            x = index + side * 0.12
            bridged = bool(row["lock_present"])

            ax.plot(
                x, row["min_distance_a"],
                marker="o", markersize=9,
                color=STATE_COLOUR[state] if bridged else SURFACE,
                markeredgecolor=STATE_COLOUR[state],
                markeredgewidth=1.8,
                zorder=3,
            )
            ax.annotate(
                row["receptor"],
                (x, row["min_distance_a"]),
                textcoords="offset points",
                xytext=(-12 if side < 0 else 12, -3),
                ha="right" if side < 0 else "left",
                fontsize=8, color=INK_SECONDARY,
            )

    ax.set_xticks(range(len(STATE_ORDER)))
    ax.set_xticklabels([s.capitalize() for s in STATE_ORDER], fontsize=10, color=INK)
    ax.set_xlim(-0.75, 1.75)
    ax.set_ylabel("Closest approach, 3.50 to 6.30 (Å)")
    ax.set_title("TM3 and TM6 separate on activation", fontsize=12, loc="left", pad=14)

    ax.plot([], [], "o", markersize=9, color=INK_MUTED,
            markeredgecolor=INK_MUTED, label="salt bridge detected")
    ax.plot([], [], "o", markersize=9, color=SURFACE,
            markeredgecolor=INK_MUTED, markeredgewidth=1.8, label="no salt bridge")
    legend = ax.legend(frameon=False, fontsize=9, loc="lower right")
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)

    fig.tight_layout()
    return fig


def fig_paired(lock):
    """Each receptor compared to itself. Removes every between-receptor difference."""
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    _style(ax)
    ax.grid(axis="y", color=GRID, linewidth=0.8)

    pivot = lock.pivot_table(
        index="receptor", columns="state", values="min_distance_a", aggfunc="first"
    )
    pivot = pivot.dropna(subset=["inactive", "active"])

    # Stack the end labels downward so receptors with similar active-state
    # distances do not print on top of each other.
    values = pivot[["inactive", "active"]].to_numpy()
    min_gap = (values.max() - values.min()) * 0.055
    label_y, previous = {}, None
    for receptor, row in pivot.sort_values("active", ascending=False).iterrows():
        y = row["active"]
        if previous is not None and previous - y < min_gap:
            y = previous - min_gap
        label_y[receptor] = y
        previous = y

    for receptor, row in pivot.iterrows():
        ax.plot([0, 1], [row["inactive"], row["active"]],
                color=INK_MUTED, linewidth=1.4, zorder=2)
        ax.plot(0, row["inactive"], "o", markersize=9, color=INACTIVE, zorder=3)
        ax.plot(1, row["active"], "o", markersize=9, color=ACTIVE, zorder=3)
        ax.annotate(receptor, (1, label_y[receptor]),
                    textcoords="offset points", xytext=(11, -3),
                    fontsize=8, color=INK_SECONDARY)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Inactive", "Active"], fontsize=10, color=INK)
    ax.set_xlim(-0.25, 1.85)
    ax.set_ylabel("Closest approach, 3.50 to 6.30 (Å)")
    ax.set_title(
        f"Every receptor moves the same way (n = {len(pivot)})",
        fontsize=12, loc="left", pad=14,
    )

    fig.tight_layout()
    return fig


def fig_bridges_per_structure(lock):
    """Total salt bridges in the analysed chain of each structure.

    Counted in the same single chain used for the lock test, so structures with
    two copies of the receptor are not double counted.
    """
    bridges = pd.read_csv(PROCESSED / "salt_bridges.csv")
    counts = (
        bridges.merge(lock[["pdb_id", "chain"]], on=["pdb_id", "chain"])
        .groupby("pdb_id")
        .size()
        .rename("n_bridges")
    )

    data = lock.merge(counts, on="pdb_id", how="left").fillna({"n_bridges": 0})
    data = data.sort_values(["state", "n_bridges"], ascending=[True, False])

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    _style(ax)

    # A gap between the two state blocks, so the grouping reads without
    # relying on colour.
    gap_after = (data["state"].values[:-1] != data["state"].values[1:]).argmax()
    positions = [i + (0.9 if i > gap_after else 0) for i in range(len(data))]

    colours = [STATE_COLOUR[s] for s in data["state"]]
    bars = ax.bar(
        positions, data["n_bridges"],
        color=colours, width=0.68, zorder=3,
        edgecolor=SURFACE, linewidth=2,
    )
    for bar, value in zip(bars, data["n_bridges"]):
        ax.annotate(f"{int(value)}", (bar.get_x() + bar.get_width() / 2, value),
                    textcoords="offset points", xytext=(0, 4),
                    ha="center", fontsize=8, color=INK_SECONDARY)

    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"{r.receptor}\n{r.pdb_id}" for r in data.itertuples()],
        fontsize=7.5, color=INK_SECONDARY,
    )
    ax.set_ylabel("Salt bridges in the analysed chain")
    ax.set_title("Salt bridges per structure", fontsize=12, loc="left", pad=14)

    for state in STATE_ORDER:
        ax.bar(0, 0, color=STATE_COLOUR[state], label=state.capitalize())
    legend = ax.legend(frameon=False, fontsize=9, loc="upper right")
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)

    fig.tight_layout()
    return fig


def fig_resolution_confound(lock):
    """Does detection depend on structure quality rather than on state?

    Side chains are placed less precisely at low resolution, so a structure
    solved at 3.5 A may yield fewer detected contacts than one at 1.8 A for
    reasons that have nothing to do with biology.
    """
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    _style(ax)

    for state in STATE_ORDER:
        subset = lock[lock["state"] == state]
        ax.plot(
            subset["resolution_a"], subset["min_distance_a"],
            "o", markersize=9, color=STATE_COLOUR[state],
            markeredgecolor=SURFACE, markeredgewidth=1.5,
            linestyle="none", label=state.capitalize(), zorder=3,
        )
        for _, row in subset.iterrows():
            ax.annotate(row["receptor"], (row["resolution_a"], row["min_distance_a"]),
                        textcoords="offset points", xytext=(10, -3),
                        fontsize=8, color=INK_SECONDARY)

    ax.set_xlabel("Resolution (Å) — lower is better")
    ax.set_ylabel("Closest approach, 3.50 to 6.30 (Å)")
    ax.set_title("Is the separation explained by structure quality?",
                 fontsize=12, loc="left", pad=14)
    # Room on the right so a label on the last point is not clipped.
    ax.margins(x=0.18)
    legend = ax.legend(frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)

    fig.tight_layout()
    return fig


THREE_LETTER = {
    "D": "Asp", "E": "Glu", "K": "Lys", "R": "Arg", "L": "Leu",
    "A": "Ala", "T": "Thr", "S": "Ser", "H": "His", "N": "Asn",
    "Q": "Gln", "G": "Gly", "V": "Val", "I": "Ile", "M": "Met",
    "F": "Phe", "W": "Trp", "Y": "Tyr", "C": "Cys", "P": "Pro",
}
ACIDIC = {"D", "E"}


def fig_per_receptor(lock):
    """One row per receptor: both states, and what sits at 6.30.

    The row label carries the explanation. A receptor with a leucine at 6.30
    has no charge to pair with 3.50, so no salt bridge is possible however
    close the two come.
    """
    rows = []
    for receptor, group in lock.groupby("receptor"):
        inactive = group[group["state"] == "inactive"]
        active = group[group["state"] == "active"]
        residue = group["aa_6x30"].iloc[0]
        rows.append(
            {
                "receptor": receptor,
                "aa_6x30": residue,
                "acidic": residue in ACIDIC,
                "inactive": inactive["min_distance_a"].iloc[0] if len(inactive) else None,
                "active": active["min_distance_a"].iloc[0] if len(active) else None,
                "lock_inactive": bool(inactive["lock_present"].iloc[0]) if len(inactive) else False,
            }
        )

    data = pd.DataFrame(rows).sort_values("inactive", na_position="last")
    data = data.reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    _style(ax)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.grid(axis="y", visible=False)

    # pandas stores a missing float as NaN, and `NaN is not None` is True,
    # so test with notna or the missing points plot silently as nothing.
    for y, row in data.iterrows():
        has_inactive = pd.notna(row["inactive"])
        has_active = pd.notna(row["active"])

        if has_inactive and has_active:
            ax.plot([row["inactive"], row["active"]], [y, y],
                    color=INK_MUTED, linewidth=1.4, zorder=2)
        if has_inactive:
            ax.plot(row["inactive"], y, "o", markersize=9, zorder=3,
                    color=INACTIVE if row["lock_inactive"] else SURFACE,
                    markeredgecolor=INACTIVE, markeredgewidth=1.8)
        if has_active:
            ax.plot(row["active"], y, "o", markersize=9, zorder=3,
                    color=SURFACE, markeredgecolor=ACTIVE, markeredgewidth=1.8)
        elif has_inactive:
            ax.annotate("6.30 not resolved in the active structure",
                        (row["inactive"], y), textcoords="offset points",
                        xytext=(14, -3), fontsize=8, color=INK_MUTED)

    labels = [
        f"{r.receptor}\n{THREE_LETTER.get(r.aa_6x30, r.aa_6x30)} at 6.30"
        + ("" if r.acidic else " — no charge")
        for r in data.itertuples()
    ]
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(labels, fontsize=8, color=INK_SECONDARY)
    ax.set_ylim(-0.6, len(data) - 0.4)
    ax.set_xlabel("Closest approach, 3.50 to 6.30 (Å)")
    ax.set_title("Each receptor, inactive to active",
                 fontsize=12, loc="left", pad=34)
    ax.margins(x=0.08)

    ax.plot([], [], "o", markersize=9, color=INACTIVE,
            markeredgecolor=INACTIVE, label="inactive, salt bridge")
    ax.plot([], [], "o", markersize=9, color=SURFACE,
            markeredgecolor=INACTIVE, markeredgewidth=1.8, label="inactive, none")
    ax.plot([], [], "o", markersize=9, color=SURFACE,
            markeredgecolor=ACTIVE, markeredgewidth=1.8, label="active")
    # Above the plot, so it cannot sit on top of a row.
    legend = ax.legend(frameon=False, fontsize=8, ncol=3,
                       loc="lower left", bbox_to_anchor=(0, 1.0))
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)

    fig.tight_layout()
    return fig


def main():
    lock = load()
    RESULTS.mkdir(parents=True, exist_ok=True)

    figures = {
        "fig1_distance_by_state": fig_distance_by_state(lock),
        "fig2_paired": fig_paired(lock),
        "fig3_bridges_per_structure": fig_bridges_per_structure(lock),
        "fig4_resolution": fig_resolution_confound(lock),
        "fig5_per_receptor": fig_per_receptor(lock),
    }

    for name, figure in figures.items():
        path = RESULTS / f"{name}.png"
        figure.savefig(path, dpi=200, facecolor=SURFACE)
        print(f"  {path}")


if __name__ == "__main__":
    main()
