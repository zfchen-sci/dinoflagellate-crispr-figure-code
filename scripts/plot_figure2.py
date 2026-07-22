from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from common import COLORS, add_bracket, configure_style, games_howell, panel_label, save_figure, significance_label


def injection_data(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Injection_data", header=21)
    df = df[df.crRNA.isin(["Site1", "Site2"])].copy()
    for label, denominator in {
        "per_injected": "Injected cells",
        "per_germinated": "Germinated cells",
        "per_viable": "Viable cells",
    }.items():
        df[label] = 100 * df["Positive cells"] / df[denominator].replace(0, np.nan)
    df["viability"] = 100 * df["Viable cells"] / df["Injected cells"]
    return df


def bar_with_points(ax, labels, groups, value, colors, ylabel, ylim=None):
    x = np.arange(len(labels))
    means = [np.nanmean(group[value]) for group in groups]
    sds = [np.nanstd(group[value], ddof=1) for group in groups]
    ax.bar(x, means, yerr=sds, capsize=2, color=colors, edgecolor=COLORS["dark"])
    for i, group in enumerate(groups):
        vals = group[value].dropna().to_numpy()
        jitter = np.linspace(-0.07, 0.07, len(vals)) if len(vals) > 1 else np.array([0.0])
        ax.scatter(x[i] + jitter, vals, s=9, c=COLORS["dark"], zorder=3)
    ax.set_xticks(x, labels)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    return x, means, sds


def plot(source_dir: Path, output_dir: Path) -> None:
    configure_style()
    source = source_dir / "Source_data_fig2.xlsx"
    raw = injection_data(source)
    edit = pd.read_excel(source, sheet_name="Fig.2d-f")
    edit["lineage"] = edit["sample"].str.rsplit("-", n=1).str[0]
    metric = "corrected_editing_efficiency_pct"

    fig, axes = plt.subplots(5, 2, figsize=(7.2, 12.0), constrained_layout=True)
    axes = axes.ravel()

    ax = axes[0]
    strategies = ["Cytosol", "UvrD+Cytosol", "Nuclei", "UvrD+Nuclei"]
    strategy_colors = ["#2679A9", "#E77E24", "#3C9A48", "#D93636"]
    doses = [5, 30, 1200]
    for strategy, color in zip(strategies, strategy_colors):
        sub = raw[(raw.crRNA == "Site2") & (raw["Injection strategy"] == strategy)]
        means = [sub.loc[sub["RNP concentration (nM)"] == dose, "viability"].mean() for dose in doses]
        sds = [sub.loc[sub["RNP concentration (nM)"] == dose, "viability"].std(ddof=1) for dose in doses]
        ax.plot(doses, means, marker="o", ms=3, color=color, label=strategy)
        ax.fill_between(doses, np.array(means) - sds, np.array(means) + sds, color=color, alpha=0.14, linewidth=0)
    ax.set_xscale("log")
    ax.set(xticks=doses, xticklabels=["5", "30", "1,200"], xlabel="RNP concentration (nM)", ylabel="Viable cells (% injected)", ylim=(0, 100))
    ax.legend(ncol=2, loc="upper left")
    panel_label(ax, "a")

    ax = axes[1]
    sub = raw[(raw.crRNA == "Site2") & (raw["Injection strategy"] == "UvrD+Nuclei")]
    x = np.arange(3)
    width = 0.24
    outcomes = [("per_injected", "per injected", "#575B60"), ("per_germinated", "per germinated", "#5C8BC9"), ("per_viable", "per viable", "#E9A15E")]
    for j, (outcome, label, color) in enumerate(outcomes):
        means = [sub.loc[sub["RNP concentration (nM)"] == dose, outcome].mean() for dose in doses]
        sds = [sub.loc[sub["RNP concentration (nM)"] == dose, outcome].std(ddof=1) for dose in doses]
        positions = x + (j - 1) * width
        ax.bar(positions, means, width, yerr=sds, capsize=2, color=color, edgecolor=COLORS["dark"], label=label)
        for i, dose in enumerate(doses):
            vals = sub.loc[sub["RNP concentration (nM)"] == dose, outcome].dropna().to_numpy()
            ax.scatter(np.full(vals.size, positions[i]), vals, s=8, c=COLORS["dark"], zorder=3)
    tukey = pairwise_tukeyhsd(sub.per_viable, sub["RNP concentration (nM)"])
    for level, ((g1, g2), p) in enumerate(zip(combinations(tukey.groupsunique, 2), tukey.pvalues)):
        i1, i2 = doses.index(int(g1)), doses.index(int(g2))
        add_bracket(ax, x[i1] + width, x[i2] + width, 78 + level * 8, significance_label(float(p)))
    ax.set(xticks=x, xticklabels=["5", "30", "1,200"], xlabel="RNP concentration (nM)", ylabel="Editing recovery (%)", ylim=(0, 108))
    ax.legend(loc="upper left")
    panel_label(ax, "b")

    ax = axes[2]
    sub = raw[(raw["Injection strategy"] == "UvrD+Nuclei") & (raw["RNP concentration (nM)"] == 1200)]
    x = np.arange(2)
    for j, (outcome, label, color) in enumerate(outcomes):
        positions = x + (j - 1) * width
        for i, site in enumerate(("Site1", "Site2")):
            vals = sub.loc[sub.crRNA == site, outcome].dropna()
            ax.bar(positions[i], vals.mean(), width, yerr=vals.std(ddof=1), capsize=2, color=color, edgecolor=COLORS["dark"])
            ax.scatter(np.full(vals.size, positions[i]), vals, s=8, c=COLORS["dark"], zorder=3)
        p = stats.ttest_ind(sub.loc[sub.crRNA == "Site1", outcome], sub.loc[sub.crRNA == "Site2", outcome], equal_var=False).pvalue
        add_bracket(ax, positions[0], positions[1], 106 + j * 10, significance_label(p))
    ax.set(xticks=x, xticklabels=["crRNA-site1", "crRNA-site2"], ylabel="Editing recovery (%)", ylim=(0, 143))
    panel_label(ax, "c")

    ax = axes[3]
    lineage = edit.groupby(["Target_site", "lineage"], as_index=False)[metric].mean()
    groups = [lineage[lineage.Target_site == site] for site in ("crRNA_site1", "crRNA_site2")]
    x, means, sds = bar_with_points(ax, ["crRNA-site1", "crRNA-site2"], groups, metric, [COLORS["blue"], COLORS["orange"]], "Editing efficiency (%)", (0, 100))
    result = stats.ttest_ind(groups[0][metric], groups[1][metric], equal_var=False)
    add_bracket(ax, 0, 1, 75, significance_label(result.pvalue))
    panel_label(ax, "d")

    ax = axes[4]
    site1 = edit[edit.Target_site == "crRNA_site1"]
    labels = list(site1.lineage.drop_duplicates())
    groups = [site1[site1.lineage == label] for label in labels]
    bar_with_points(ax, labels, groups, metric, plt.cm.Blues(np.linspace(0.25, 0.7, len(labels))), "Editing efficiency (%)", (0, 100))
    comparisons = games_howell({label: site1.loc[site1.lineage == label, metric] for label in labels})
    significant = [row for row in comparisons if row["p_adjusted"] < 0.05]
    for level, row in enumerate(significant):
        add_bracket(
            ax,
            labels.index(row["group_1"]),
            labels.index(row["group_2"]),
            77 + level * 7.0,
            significance_label(row["p_adjusted"]),
            height=1.7,
        )
    ax.tick_params(axis="x", rotation=30)
    panel_label(ax, "e1")

    ax = axes[5]
    site2 = edit[edit.Target_site == "crRNA_site2"]
    labels = list(site2.lineage.drop_duplicates())
    groups = [site2[site2.lineage == label] for label in labels]
    bar_with_points(ax, labels, groups, metric, plt.cm.Oranges(np.linspace(0.2, 0.75, len(labels))), "Editing efficiency (%)", (0, 120))
    comparisons = games_howell({label: site2.loc[site2.lineage == label, metric] for label in labels})
    significant = [row for row in comparisons if row["p_adjusted"] < 0.05]
    for level, row in enumerate(significant):
        add_bracket(
            ax,
            labels.index(row["group_1"]),
            labels.index(row["group_2"]),
            70 + level * 7.0,
            significance_label(row["p_adjusted"]),
            height=1.7,
        )
    ax.tick_params(axis="x", rotation=30)
    panel_label(ax, "e2")

    ax = axes[6]
    seq = edit.groupby(["Target_site", "lineage"])[["normalized_wt_pct", "normalized_mutant_like_pct", "normalized_mutant_pct"]].mean().reset_index()
    labels = seq.lineage.tolist()
    x = np.arange(len(labels))
    bottom = np.zeros(len(seq))
    for col, label, color in (
        ("normalized_mutant_pct", "Mutant", "#638CB6"),
        ("normalized_mutant_like_pct", "Mutant-like", "#D49A64"),
        ("normalized_wt_pct", "WT", "#D8D4CB"),
    ):
        ax.bar(x, seq[col], bottom=bottom, color=color, edgecolor="white", linewidth=0.4, label=label)
        bottom += seq[col].to_numpy()
    ax.axvline(4.5, color="#D0D0D0", lw=0.7)
    ax.set(xticks=x, xticklabels=labels, ylabel="Normalized sequence type (%)", ylim=(0, 100))
    ax.tick_params(axis="x", rotation=35)
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    panel_label(ax, "f")

    ax = axes[7]
    repair = pd.read_excel(source, sheet_name="Fig.2g")
    pivot = repair.pivot(index="Site", columns="repair_type", values="percent_of_mutant_reads").fillna(0).reindex(["site1", "site2"])
    bottom = np.zeros(len(pivot))
    for label, color in (("NHEJ", "#55A9D8"), ("MMEJ", "#E4A300"), ("Other", "#BABABA")):
        values = pivot[label] if label in pivot else np.zeros(len(pivot))
        ax.bar(np.arange(len(pivot)), values, bottom=bottom, color=color, edgecolor=COLORS["dark"], linewidth=0.5, label=label)
        bottom += np.asarray(values)
    ax.set(xticks=[0, 1], xticklabels=["site1", "site2"], ylabel="Repair-type reads (%)", ylim=(0, 100))
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    panel_label(ax, "g")

    ax = axes[8]
    candidates = pd.read_excel(source, sheet_name="Fig.2h")
    ax.bar(candidates.Category, candidates.candidate_items_n, color=["#58A9D7", "#E4A300", "#13A47A"], edgecolor=COLORS["dark"])
    ax.set_ylabel("Annotated candidate proteins (n)")
    ax.tick_params(axis="x", rotation=20)
    panel_label(ax, "h")

    ax = axes[9]
    motif = pd.read_excel(source, sheet_name="Fig.2i")
    pivot = motif.pivot(index=["Site", "sample"], columns="motif_length_bp", values="percent_of_sample_mutant_reads").fillna(0).reset_index()
    x = np.arange(len(pivot))
    bottom = np.zeros(len(pivot))
    for length, color in ((2, "#1079A9"), (3, "#E2A400"), (4, "#149C73")):
        ax.bar(x, pivot[length], bottom=bottom, color=color, width=0.82, label=f"{length} bp")
        bottom += pivot[length].to_numpy()
    split = int((pivot.Site == "site1").sum()) - 0.5
    ax.axvline(split, color="#BFBFBF", lw=0.7)
    ax.set(xticks=x, xticklabels=pivot["sample"], ylabel="MMEJ reads (%)")
    ax.tick_params(axis="x", rotation=90, labelsize=5)
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    panel_label(ax, "i")

    save_figure(fig, output_dir, "Figure2_data_panels")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
