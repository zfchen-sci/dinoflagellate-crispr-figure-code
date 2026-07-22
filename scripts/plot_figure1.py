from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import COLORS, add_bracket, configure_style, panel_label, save_figure, significance_label


def plot(source_dir: Path, output_dir: Path) -> None:
    configure_style()
    source = source_dir / "Source_data_fig1.xlsx"
    tip = pd.read_excel(source, sheet_name="Fig.1d")
    uvrd = pd.read_excel(source, sheet_name="Fig.1f")
    comp = pd.read_excel(source, sheet_name="Fig.1g_base_comp", header=2)

    fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])

    ax = fig.add_subplot(grid[0, 0])
    groups = [500, 50]
    x = np.arange(2)
    width = 0.34
    for offset, outcome, label, color in (
        (-width / 2, "germination_rate", "Germination", COLORS["blue"]),
        (width / 2, "viability_rate", "Viable", COLORS["teal"]),
    ):
        means = [tip.loc[tip.tip_size_nm == group, outcome].mean() for group in groups]
        sds = [tip.loc[tip.tip_size_nm == group, outcome].std(ddof=1) for group in groups]
        ax.bar(x + offset, means, width, yerr=sds, capsize=2, color=color, edgecolor=COLORS["dark"], label=label)
        for i, group in enumerate(groups):
            values = tip.loc[tip.tip_size_nm == group, outcome].to_numpy()
            ax.scatter(np.full(values.size, x[i] + offset), values, s=10, c=COLORS["dark"], zorder=3)
        p = stats.ttest_ind(
            tip.loc[tip.tip_size_nm == 500, outcome], tip.loc[tip.tip_size_nm == 50, outcome], equal_var=False
        ).pvalue
        y = 79 if outcome == "germination_rate" else 88
        add_bracket(ax, x[0] + offset, x[1] + offset, y, significance_label(p, tiers=1), height=2)
    ax.set(xticks=x, xticklabels=["500", "50"], xlabel="Micropipette tip size (nm)", ylabel="Percentage (%)", ylim=(0, 100))
    ax.legend(loc="upper left")
    panel_label(ax, "d")

    ax = fig.add_subplot(grid[0, 1])
    groups = [0.003, 0.012]
    x = np.arange(2)
    for offset, outcome, label, color in (
        (-width / 2, "germination_rate", "Germination", COLORS["blue"]),
        (width / 2, "viability_rate", "Viable", COLORS["teal"]),
    ):
        means = [uvrd.loc[np.isclose(uvrd.uvrd_uM, group), outcome].mean() for group in groups]
        sds = [uvrd.loc[np.isclose(uvrd.uvrd_uM, group), outcome].std(ddof=1) for group in groups]
        ax.bar(x + offset, means, width, yerr=sds, capsize=2, color=color, edgecolor=COLORS["dark"], label=label)
        for i, group in enumerate(groups):
            values = uvrd.loc[np.isclose(uvrd.uvrd_uM, group), outcome].to_numpy()
            ax.scatter(np.full(values.size, x[i] + offset), values, s=10, c=COLORS["dark"], zorder=3)
        p = stats.ttest_ind(
            uvrd.loc[np.isclose(uvrd.uvrd_uM, 0.003), outcome],
            uvrd.loc[np.isclose(uvrd.uvrd_uM, 0.012), outcome],
            equal_var=False,
        ).pvalue
        y = 73 if outcome == "germination_rate" else 82
        add_bracket(ax, x[0] + offset, x[1] + offset, y, significance_label(p, tiers=1), height=2)
    ax.set(xticks=x, xticklabels=["0.003", "0.012"], xlabel="UvrD concentration (µM)", ylabel="Percentage (%)", ylim=(0, 95))
    panel_label(ax, "f")

    ax = fig.add_subplot(grid[1, :])
    nucleotide_colors = {"A": "#2A85B8", "C": "#3A9B42", "G": "#D42A7B", "T": "#E5A60A"}
    window_start, window_end = 2180, 2870
    snp = comp.loc[
        comp["Position"].between(window_start, window_end)
        & ((comp[["A", "T", "C", "G"]] > 0.01).sum(axis=1) > 1)
    ].copy()
    crRNA_intervals = {
        "crRNA-site1": (2222, 2243),
        "crRNA-site2": (2446, 2467),
    }
    for index, (label, (start, end)) in enumerate(crRNA_intervals.items()):
        ax.axvspan(
            start,
            end,
            color="#23D9E1",
            alpha=0.9,
            linewidth=0,
            label="crRNA" if index == 0 else None,
            zorder=0,
        )
        ax.text(
            (start + end) / 2,
            1.035,
            label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=6,
            fontweight="bold",
            clip_on=False,
        )
    for base in ("A", "C", "G", "T"):
        y = snp[base].where(snp[base] > 0.01)
        ax.scatter(snp.Position, y, s=8, color=nucleotide_colors[base], label=base, linewidths=0, zorder=2)
    major_ticks = [2180, 2480, 2780]
    minor_ticks = [2280, 2380, 2580, 2680, 2870]
    ax.set(
        xlabel="SNP positions in the sxtA4 gene (bp)",
        ylabel="Relative abundance (%)",
        xlim=(window_start, window_end),
        ylim=(-5, 105),
        xticks=major_ticks,
    )
    ax.set_xticks(minor_ticks, minor=True)
    ax.tick_params(axis="x", which="minor", length=2.5, labelbottom=False)
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(base) for base in ("A", "C", "G", "T")] + [labels.index("crRNA")]
    ax.legend(
        [handles[index] for index in order],
        ["A (blue)", "C (green)", "G (magenta)", "T (yellow)", "crRNA"],
        ncol=5,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.19),
        columnspacing=1.0,
        handletextpad=0.35,
    )
    panel_label(ax, "h")

    save_figure(fig, output_dir, "Figure1_data_panels")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
