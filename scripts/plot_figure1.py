from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import configure_style, panel_label, save_figure


GUIDE_ORDER = [
    "Site1_A4_C19",
    "Site1_A4_U19",
    "Site1_C4_C19",
    "Site1_C4_U19",
    "Site2",
]


def validate_columns(data: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def plot(source_dir: Path, output_dir: Path) -> None:
    """Recreate only the quantitative Fig. 1 panels that are appropriate to redraw.

    Fig. 1a is intentionally excluded. Fig. 1e is also excluded because the
    published panel contains RNA secondary-structure diagrams; an MFE bar chart
    would not be a faithful substitute.
    """
    configure_style()
    source = source_dir / "Source_data_fig1.xlsx"
    composition = pd.read_excel(source, sheet_name="Fig.1a-b_base_comp", header=2)
    pairs = pd.read_excel(source, sheet_name="Fig.1f_pair_prob")

    validate_columns(composition, {"Position", "A", "T", "C", "G"}, "Fig. 1b")
    validate_columns(
        pairs,
        {"guide_id", "guide_label", "i", "j", "pair_probability"},
        "Fig. 1f",
    )

    fig = plt.figure(figsize=(7.2, 4.6))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[1.12, 1.0],
        left=0.09,
        right=0.98,
        bottom=0.10,
        top=0.93,
        hspace=0.55,
    )

    axis_b = fig.add_subplot(outer[0, 0])
    nucleotide_colors = {
        "A": "#2388B8",
        "C": "#2CA02C",
        "G": "#EC168C",
        "T": "#E6A600",
    }
    snp = composition.loc[
        composition["Position"].between(2180, 2870)
        & ((composition[["A", "T", "C", "G"]] > 0.01).sum(axis=1) > 1)
    ].copy()
    for start, end in ((2222, 2243), (2446, 2467)):
        axis_b.axvspan(start, end, color="#23D9E1", alpha=0.72, linewidth=0, zorder=0)
    for base in ("A", "C", "G", "T"):
        values = snp[base].where(snp[base] > 0.01)
        axis_b.scatter(
            snp["Position"],
            values,
            s=10,
            color=nucleotide_colors[base],
            label=base,
            linewidths=0,
            zorder=2,
        )
    axis_b.set(
        xlabel=r"SNP positions in the $sxtA4$ gene (bp)",
        ylabel="Relative abundance (%)",
        xlim=(2180, 2870),
        ylim=(-4, 104),
        xticks=[2180, 2480, 2780],
        yticks=[0, 20, 40, 60, 80, 100],
    )
    axis_b.legend(
        ncol=4,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        columnspacing=1.4,
        handletextpad=0.35,
    )
    panel_label(axis_b, "b")

    pair_grid = outer[1, 0].subgridspec(1, 5, wspace=0.38)
    for index, guide in enumerate(GUIDE_ORDER):
        axis = fig.add_subplot(pair_grid[0, index])
        subset = pairs.loc[
            (pairs["guide_id"] == guide) & (pairs["pair_probability"] >= 0.01)
        ].copy()
        if subset.empty:
            raise ValueError(f"Fig. 1f has no base-pair probabilities for {guide}.")
        probability = subset["pair_probability"].to_numpy(dtype=float)
        axis.scatter(
            subset["i"],
            subset["j"],
            s=4 + 35 * probability,
            c=probability,
            cmap="Blues",
            vmin=0,
            vmax=1,
            alpha=0.85,
            edgecolors="none",
        )
        axis.plot([1, 42], [1, 42], color="#D9E1E5", lw=0.45)
        for boundary in (21.5, 29.5):
            axis.axvline(boundary, color="#D9E1E5", lw=0.45, ls=":")
            axis.axhline(boundary, color="#D9E1E5", lw=0.45, ls=":")
        axis.set(
            xlim=(0.5, 42.5),
            ylim=(42.5, 0.5),
            aspect="equal",
            title=str(subset["guide_label"].iloc[0]),
            xticks=[1, 21, 42],
            yticks=[1, 21, 42] if index == 0 else [],
            xlabel="i",
        )
        if index == 0:
            axis.set_ylabel("j")
            panel_label(axis, "f")

    save_figure(fig, output_dir, "Figure1_data_panels")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recreate current Fig. 1b and Fig. 1f; Fig. 1a and Fig. 1e are intentionally excluded."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
