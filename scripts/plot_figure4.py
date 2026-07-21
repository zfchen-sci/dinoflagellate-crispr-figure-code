from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from common import COLORS, add_bracket, configure_style, panel_label, save_figure, significance_label


WT = "Wild-type strain"
EDITED = "sxtA4-edited strain"


def plot(source_dir: Path, output_dir: Path) -> None:
    configure_style()
    source = source_dir / "Source_data_fig4.xlsx"
    annotation = pd.read_excel(source, sheet_name="Fig.4a")
    sxta = pd.read_excel(source, sheet_name="Fig.4b")
    pst = pd.read_excel(source, sheet_name="Fig.4c")
    deg = pd.read_excel(source, sheet_name="Fig.4d")
    persistent = pd.read_excel(source, sheet_name="Fig.4e")

    fig = plt.figure(figsize=(7.2, 8.4), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.35])

    ax = fig.add_subplot(grid[0, 0])
    labels = annotation.iloc[:, 0].astype(str).tolist()
    values = annotation.iloc[:, 1].astype(float).to_numpy()
    colors = [COLORS["teal"], "#50545A"]
    wedges, _ = ax.pie(values, colors=colors, startangle=90, counterclock=False, wedgeprops={"width": 0.42, "edgecolor": "white"})
    ax.text(0, 0, f"{int(values.sum()):,}\ntotal genes", ha="center", va="center", fontsize=8)
    legend = [f"{label}: {int(value):,} ({value / values.sum() * 100:.1f}%)" for label, value in zip(labels, values)]
    ax.legend(wedges, legend, loc="lower center", bbox_to_anchor=(0.5, -0.14))
    ax.set_title("Overall annotation status")
    panel_label(ax, "a")

    ax = fig.add_subplot(grid[0, 1])
    days = sorted(sxta.day.unique())
    x = np.arange(len(days))
    width = 0.34
    for offset, group, color in ((-width / 2, WT, COLORS["wt"]), (width / 2, EDITED, COLORS["edited"])):
        means = [sxta.loc[(sxta.group == group) & (sxta.day == day), "fpkm"].mean() for day in days]
        sds = [sxta.loc[(sxta.group == group) & (sxta.day == day), "fpkm"].std(ddof=1) for day in days]
        ax.bar(x + offset, means, width, yerr=sds, capsize=2, color=color, edgecolor=COLORS["dark"], label=group)
        for i, day in enumerate(days):
            vals = sxta.loc[(sxta.group == group) & (sxta.day == day), "fpkm"].to_numpy()
            ax.scatter(x[i] + offset + np.linspace(-0.04, 0.04, len(vals)), vals, s=8, color=COLORS["dark"], zorder=3)
    raw_p = []
    for day in days:
        sub = sxta[sxta.day == day]
        raw_p.append(stats.ttest_ind(sub.loc[sub.group == WT, "fpkm"], sub.loc[sub.group == EDITED, "fpkm"], equal_var=False).pvalue)
    adjusted = multipletests(raw_p, method="holm")[1]
    for i, (day, p) in enumerate(zip(days, adjusted)):
        sub = sxta[sxta.day == day]
        fold = sub.loc[sub.group == WT, "fpkm"].mean() / sub.loc[sub.group == EDITED, "fpkm"].mean()
        y = sub.fpkm.max() * 1.08
        add_bracket(ax, x[i] - width / 2, x[i] + width / 2, y, significance_label(float(p)), height=2)
        ax.text(x[i], y + 8, f"{fold:.1f}×", ha="center", va="bottom")
    ax.set(xticks=x, xticklabels=[f"{day}d" for day in days], xlabel="Time", ylabel="sxtA expression (FPKM)", ylim=(0, 150))
    ax.legend(loc="upper left")
    panel_label(ax, "b")

    ax = fig.add_subplot(grid[1, 0])
    genes = list(pst.pst_gene.drop_duplicates())
    x = np.arange(len(genes))
    raw_p = []
    for gene in genes:
        sub = pst[pst.pst_gene == gene]
        raw_p.append(stats.ttest_ind(sub.loc[sub.group == WT, "fpkm"], sub.loc[sub.group == EDITED, "fpkm"], equal_var=False).pvalue)
    adjusted = multipletests(raw_p, method="holm")[1]
    for offset, group, color in ((-width / 2, WT, COLORS["wt"]), (width / 2, EDITED, COLORS["edited"])):
        means = [pst.loc[(pst.pst_gene == gene) & (pst.group == group), "fpkm"].mean() for gene in genes]
        sds = [pst.loc[(pst.pst_gene == gene) & (pst.group == group), "fpkm"].std(ddof=1) for gene in genes]
        ax.bar(x + offset, means, width, yerr=sds, capsize=1.5, color=color, edgecolor=COLORS["dark"], label=group)
    ax.set_yscale("symlog", linthresh=1)
    for i, (gene, p) in enumerate(zip(genes, adjusted)):
        sub = pst[pst.pst_gene == gene]
        y = max(sub.fpkm.max() * 1.15, 0.25)
        add_bracket(ax, x[i] - width / 2, x[i] + width / 2, y, significance_label(float(p)), height=max(y * 0.08, 0.05))
    ax.set(xticks=x, xticklabels=[f"${gene}$" for gene in genes], ylabel="Gene expression (FPKM)")
    ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.02), loc="lower center")
    panel_label(ax, "c")

    ax = fig.add_subplot(grid[1, 1])
    x = np.arange(len(deg))
    ax.bar(x - width / 2, deg.iloc[:, 1], width, color=COLORS["navy"], label="Up in edited")
    ax.bar(x + width / 2, deg.iloc[:, 2], width, color=COLORS["wt"], label="Down in edited")
    ax.set(xticks=x, xticklabels=[str(v) for v in deg.iloc[:, 0]], xlabel="Time (d)", ylabel="Number of DEGs")
    ax.legend(loc="upper left")
    panel_label(ax, "d")

    ax = fig.add_subplot(grid[2, :])
    columns = [
        "plot_WT_5d_vs_M_5d_log2fc",
        "plot_WT_9d_vs_M_9d_log2fc",
        "plot_WT_13d_vs_M_13d_log2fc",
        "plot_WT_15d_vs_M_15d_log2fc",
    ]
    order = persistent.sort_values(["response_class", "mean_log2fc"], ascending=[True, False])
    matrix = order[columns].to_numpy(dtype=float)
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-6, vmax=6, interpolation="nearest", rasterized=True)
    ax.set(xticks=np.arange(4), xticklabels=["5", "9", "13", "15"], xlabel="Time (d)", ylabel="Persistent DEGs")
    ax.set_yticks([])
    colorbar = fig.colorbar(image, ax=ax, fraction=0.018, pad=0.02)
    colorbar.set_label("log2 fold change M/WT")
    dark = int(persistent.dark_gene.sum())
    ax.set_title(f"626 shared DEGs; {dark} unannotated")
    panel_label(ax, "e")

    save_figure(fig, output_dir, "Figure4_data_panels")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
