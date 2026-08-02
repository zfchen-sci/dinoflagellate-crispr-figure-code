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
    data = pd.read_excel(path, sheet_name="Fig.3f-g_raw", header=21)
    data = data[data["crRNA"].isin(["Site1", "Site2"])].copy()
    for label, denominator in {
        "per_injected": "Injected cells",
        "per_germinated": "Germinated cells",
        "per_viable": "Viable cells",
    }.items():
        data[label] = 100 * data["Positive cells"] / data[denominator].replace(0, np.nan)
    return data


def bar_with_points(ax, labels, groups, value, colors, ylabel, ylim=None):
    x = np.arange(len(labels))
    means = [np.nanmean(group[value]) for group in groups]
    sds = [np.nanstd(group[value], ddof=1) for group in groups]
    ax.bar(x, means, yerr=sds, capsize=2, color=colors, edgecolor=COLORS["dark"])
    for i, group in enumerate(groups):
        values = group[value].dropna().to_numpy(dtype=float)
        jitter = np.linspace(-0.07, 0.07, len(values)) if len(values) > 1 else np.array([0.0])
        ax.scatter(x[i] + jitter, values, s=9, c=COLORS["dark"], zorder=3)
    ax.set_xticks(x, labels)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    return x


def cell_recovery_panel(ax: plt.Axes, raw: pd.DataFrame, strategy: str, panel: str) -> None:
    subset = raw.loc[raw["Injection strategy"] == strategy].copy()
    doses = [5, 30, 1200]
    x = np.arange(3)
    jitter = np.array([-0.035, 0.0, 0.035])
    for offset, column, label, color in (
        (-0.10, "Germination (% injected)", "Germinated", "#4E8FDF"),
        (0.10, "Viable-cell recovery (% injected)", "Viable", "#F3A055"),
    ):
        for dose_index, dose in enumerate(doses):
            values = subset.loc[
                subset["RNP concentration (nM)"] == dose, column
            ].dropna().to_numpy(dtype=float)
            if len(values) != 3:
                raise ValueError(
                    f"Fig. 3{panel}: {strategy}, {dose} nM has {len(values)} values; expected 3."
                )
            position = x[dose_index] + offset
            ax.scatter(
                position + jitter,
                values,
                s=13,
                facecolor=color,
                edgecolor=COLORS["dark"],
                linewidth=0.65,
                alpha=0.92,
                zorder=4,
                label=label if dose_index == 0 else None,
            )
            ax.errorbar(
                position,
                values.mean(),
                yerr=values.std(ddof=1),
                fmt="_",
                markersize=8,
                markeredgewidth=1.4,
                color=color,
                capsize=2,
                elinewidth=0.85,
                zorder=5,
            )
    display = {
        "Cytosol": "Cytosolic RNP",
        "Nuclei": "Nuclear RNP",
        "UvrD+Cytosol": "Cytosolic RNP + UvrD",
        "UvrD+Nuclei": "Nuclear RNP + UvrD",
    }[strategy]
    ax.set(
        xticks=x,
        xticklabels=["5", "30", "1,200"],
        xlabel="Cas12a RNP concentration (nM)",
        ylabel="Recovery (% of injected cells)",
        ylim=(0, 105),
        yticks=[0, 25, 50, 75, 100],
        title=display,
    )
    panel_label(ax, panel)


def plot(source_dir: Path, output_dir: Path) -> None:
    configure_style()
    source = source_dir / "Source_data_fig3.xlsx"
    cell_raw = pd.read_excel(source, sheet_name="Fig.3b-e_raw", header=3)
    required_cell_columns = {"Injection strategy", "RNP concentration (nM)", "Germination (% injected)", "Viable-cell recovery (% injected)"}
    if not required_cell_columns.issubset(cell_raw.columns):
        raise ValueError("Fig. 3b–e source-data columns are incomplete.")
    raw = injection_data(source)
    edit = pd.read_excel(source, sheet_name="Fig.3h-j")
    edit["lineage"] = edit["sample"].str.rsplit("-", n=1).str[0]
    metric = "corrected_editing_efficiency_pct"

    fig, axes = plt.subplots(7, 2, figsize=(7.2, 16.0), constrained_layout=True)
    axes = axes.ravel()

    for ax, strategy, panel in zip(axes[:4], ["Cytosol", "Nuclei", "UvrD+Cytosol", "UvrD+Nuclei"], ["b", "c", "d", "e"]):
        cell_recovery_panel(ax, cell_raw, strategy, panel)
    axes[0].legend(loc="upper left")

    ax = axes[4]
    subset = raw[(raw["crRNA"] == "Site2") & (raw["Injection strategy"] == "UvrD+Nuclei")]
    doses = [5, 30, 1200]
    x = np.arange(3)
    width = 0.24
    outcomes = [
        ("per_injected", "per injected", "#575B60"),
        ("per_germinated", "per germinated", "#5C8BC9"),
        ("per_viable", "per viable", "#E9A15E"),
    ]
    for j, (outcome, label, color) in enumerate(outcomes):
        positions = x + (j - 1) * width
        groups = [subset.loc[subset["RNP concentration (nM)"] == dose, outcome].dropna() for dose in doses]
        ax.bar(positions, [group.mean() for group in groups], width, yerr=[group.std(ddof=1) for group in groups], capsize=2, color=color, edgecolor=COLORS["dark"], label=label)
        for position, group in zip(positions, groups):
            ax.scatter(np.full(len(group), position), group, s=8, c=COLORS["dark"], zorder=3)
    tukey = pairwise_tukeyhsd(subset["per_viable"], subset["RNP concentration (nM)"])
    for level, ((g1, g2), p_value) in enumerate(zip(combinations(tukey.groupsunique, 2), tukey.pvalues)):
        add_bracket(ax, x[doses.index(int(g1))] + width, x[doses.index(int(g2))] + width, 78 + level * 8, significance_label(float(p_value)))
    ax.set(xticks=x, xticklabels=["5", "30", "1,200"], xlabel="RNP concentration (nM)", ylabel="Edited-lineage recovery (%)", ylim=(0, 108))
    ax.legend(loc="upper left")
    panel_label(ax, "f")

    ax = axes[5]
    subset = raw[(raw["Injection strategy"] == "UvrD+Nuclei") & (raw["RNP concentration (nM)"] == 1200)]
    x = np.arange(2)
    for j, (outcome, label, color) in enumerate(outcomes):
        positions = x + (j - 1) * width
        for i, site in enumerate(("Site1", "Site2")):
            values = subset.loc[subset["crRNA"] == site, outcome].dropna()
            ax.bar(positions[i], values.mean(), width, yerr=values.std(ddof=1), capsize=2, color=color, edgecolor=COLORS["dark"], label=label if i == 0 else None)
            ax.scatter(np.full(len(values), positions[i]), values, s=8, c=COLORS["dark"], zorder=3)
        p_value = stats.ttest_ind(subset.loc[subset["crRNA"] == "Site1", outcome], subset.loc[subset["crRNA"] == "Site2", outcome], equal_var=False).pvalue
        add_bracket(ax, positions[0], positions[1], 106 + j * 10, significance_label(p_value))
    ax.set(xticks=x, xticklabels=["crRNA-site1", "crRNA-site2"], ylabel="Edited-lineage recovery (%)", ylim=(0, 143))
    ax.legend(loc="upper left")
    panel_label(ax, "g")

    ax = axes[6]
    lineage = edit.groupby(["Target_site", "lineage"], as_index=False)[metric].mean()
    groups = [lineage[lineage["Target_site"] == site] for site in ("crRNA_site1", "crRNA_site2")]
    bar_with_points(ax, ["crRNA-site1", "crRNA-site2"], groups, metric, [COLORS["blue"], COLORS["orange"]], "Copy-level editing (%)", (0, 100))
    result = stats.ttest_ind(groups[0][metric], groups[1][metric], equal_var=False)
    add_bracket(ax, 0, 1, 75, significance_label(result.pvalue))
    panel_label(ax, "h")

    for axis_index, site, cmap in ((7, "crRNA_site1", plt.cm.Blues), (8, "crRNA_site2", plt.cm.Oranges)):
        ax = axes[axis_index]
        site_data = edit[edit["Target_site"] == site]
        labels = list(site_data["lineage"].drop_duplicates())
        groups = [site_data[site_data["lineage"] == label] for label in labels]
        bar_with_points(ax, labels, groups, metric, cmap(np.linspace(0.25, 0.75, len(labels))), "Copy-level editing (%)", (0, 110))
        comparisons_result = games_howell({label: site_data.loc[site_data["lineage"] == label, metric] for label in labels})
        for level, comparison in enumerate(row for row in comparisons_result if row["p_adjusted"] < 0.05):
            add_bracket(ax, labels.index(comparison["group_1"]), labels.index(comparison["group_2"]), 77 + level * 7, significance_label(comparison["p_adjusted"]), height=1.7)
        ax.tick_params(axis="x", rotation=30)
        if axis_index == 7:
            panel_label(ax, "i")

    ax = axes[9]
    sequence = edit.groupby(["Target_site", "lineage"])[["normalized_wt_pct", "normalized_mutant_like_pct", "normalized_mutant_pct"]].mean().reset_index()
    labels = sequence["lineage"].tolist()
    x = np.arange(len(labels))
    bottom = np.zeros(len(sequence))
    for column, label, color in (
        ("normalized_mutant_pct", "Mutant", "#638CB6"),
        ("normalized_mutant_like_pct", "Mutant-like", "#D49A64"),
        ("normalized_wt_pct", "WT", "#D8D4CB"),
    ):
        ax.bar(x, sequence[column], bottom=bottom, color=color, edgecolor="white", linewidth=0.4, label=label)
        bottom += sequence[column].to_numpy()
    ax.axvline(4.5, color="#D0D0D0", lw=0.7)
    ax.set(xticks=x, xticklabels=labels, ylabel="Normalized sequence type (%)", ylim=(0, 100))
    ax.tick_params(axis="x", rotation=35)
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    panel_label(ax, "j")

    ax = axes[10]
    repair = pd.read_excel(source, sheet_name="Fig.3k")
    pivot = repair.pivot(index="Site", columns="repair_type", values="percent_of_mutant_reads").fillna(0).reindex(["site1", "site2"])
    bottom = np.zeros(len(pivot))
    for label, color in (("NHEJ", "#55A9D8"), ("MMEJ", "#E4A300"), ("Other", "#BABABA")):
        values = pivot[label] if label in pivot else np.zeros(len(pivot))
        ax.bar(np.arange(len(pivot)), values, bottom=bottom, color=color, edgecolor=COLORS["dark"], linewidth=0.5, label=label)
        bottom += np.asarray(values)
    ax.set(xticks=[0, 1], xticklabels=["crRNA-site1", "crRNA-site2"], ylabel="Junction class (%)", ylim=(0, 100))
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    panel_label(ax, "k")

    ax = axes[11]
    candidates = pd.read_excel(source, sheet_name="Fig.3l")
    ax.bar(candidates["Category"], candidates["candidate_items_n"], color=["#58A9D7", "#E4A300", "#13A47A"], edgecolor=COLORS["dark"])
    ax.set_ylabel("Annotated candidate proteins (n)")
    ax.tick_params(axis="x", rotation=20)
    panel_label(ax, "l")

    ax = axes[12]
    motif = pd.read_excel(source, sheet_name="Fig.3m")
    motif_pivot = motif.pivot(index=["Site", "sample"], columns="motif_length_bp", values="percent_of_sample_mutant_reads").fillna(0).reset_index()
    x = np.arange(len(motif_pivot))
    bottom = np.zeros(len(motif_pivot))
    for length, color in ((2, "#1079A9"), (3, "#E2A400"), (4, "#149C73")):
        ax.bar(x, motif_pivot[length], bottom=bottom, color=color, width=0.82, label=f"{length} bp")
        bottom += motif_pivot[length].to_numpy()
    ax.axvline(int((motif_pivot["Site"] == "site1").sum()) - 0.5, color="#BFBFBF", lw=0.7)
    ax.set(xticks=x, xticklabels=motif_pivot["sample"], ylabel="Mutant reads with microhomology (%)")
    ax.tick_params(axis="x", rotation=90, labelsize=5)
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    panel_label(ax, "m")
    axes[13].axis("off")

    save_figure(fig, output_dir, "Figure3_data_panels")


def main() -> None:
    parser = argparse.ArgumentParser(description="Recreate data-driven panels for current Figure 3.")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
