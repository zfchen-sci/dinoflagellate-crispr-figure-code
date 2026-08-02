from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path as MplPath
from scipy import stats

from common import COLORS, configure_style, panel_label, save_figure


def plot_ed1b(source: Path, output_dir: Path) -> None:
    data = pd.read_excel(source, sheet_name="ED_Fig.1b", header=5)
    required = {"log10(sxtA4-copies)", "Ct Value"}
    if not required.issubset(data.columns):
        raise ValueError("ED_Fig.1b source-data columns are incomplete.")
    means = (
        data.dropna(subset=list(required))
        .groupby("log10(sxtA4-copies)", as_index=False)["Ct Value"]
        .mean()
        .sort_values("log10(sxtA4-copies)")
    )
    if len(means) < 3:
        raise ValueError("ED_Fig.1b requires Ct values for at least three copy-number levels.")

    x = means["log10(sxtA4-copies)"].to_numpy(dtype=float)
    y = means["Ct Value"].to_numpy(dtype=float)
    regression = stats.linregress(x, y)
    fitted_x = np.linspace(x.min(), x.max(), 300)
    fitted_y = regression.intercept + regression.slope * fitted_x
    residuals = y - (regression.intercept + regression.slope * x)
    degrees = len(x) - 2
    residual_sd = np.sqrt(np.sum(residuals**2) / degrees)
    distance = (fitted_x - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2)
    critical = stats.t.ppf(0.975, degrees)
    confidence = critical * residual_sd * np.sqrt(1 / len(x) + distance)
    prediction = critical * residual_sd * np.sqrt(1 + 1 / len(x) + distance)

    fig, ax = plt.subplots(figsize=(3.5, 2.8), constrained_layout=True)
    prediction_band = ax.fill_between(fitted_x, fitted_y - prediction, fitted_y + prediction, color="#F9DADA", alpha=0.65, linewidth=0)
    confidence_band = ax.fill_between(fitted_x, fitted_y - confidence, fitted_y + confidence, color="#F1A5A5", alpha=0.70, linewidth=0)
    ax.plot(fitted_x, fitted_y, color="#C83E3E", linewidth=1.1)
    ax.scatter(x, y, marker="s", s=14, color="black", linewidths=0, zorder=3)
    ax.legend(
        [confidence_band, prediction_band],
        ["95% confidence band", "95% prediction band"],
        loc="upper right",
        fontsize=5.8,
        handlelength=1.8,
        labelspacing=0.25,
    )
    ax.text(
        0.10,
        0.18,
        rf"$y = {regression.slope:.2f}x + {regression.intercept:.2f}$" + "\n" + rf"$R^2 = {regression.rvalue**2:.3f}$",
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
    )
    ax.set(
        xlim=(-0.15, 9.15),
        ylim=(5, 42),
        xticks=np.arange(10),
        yticks=[10, 20, 30, 40],
        xlabel=r"log$_{10}$(copies)",
        ylabel="Ct value",
    )
    panel_label(ax, "b")
    save_figure(fig, output_dir, "Extended_Data_Figure1b")


def plot_ed1c(source: Path, output_dir: Path) -> None:
    data = pd.read_excel(source, sheet_name="ED_Fig.1c", header=2)
    biological = data.groupby(["Number of cells", "Samples"], as_index=False)["Ct Value"].mean()
    biological["copies_per_cell"] = 10 ** ((39.67 - biological["Ct Value"]) / 3.59) * 80 / biological["Number of cells"]
    groups = sorted(biological["Number of cells"].unique())
    fig, ax = plt.subplots(figsize=(3.5, 2.8), constrained_layout=True)
    x = np.arange(len(groups))
    means = [biological.loc[biological["Number of cells"] == group, "copies_per_cell"].mean() for group in groups]
    sds = [biological.loc[biological["Number of cells"] == group, "copies_per_cell"].std(ddof=1) for group in groups]
    ax.bar(x, means, yerr=sds, capsize=2, color=["#EAB07A", "#8FC58F", "#B2A2CB"], edgecolor=COLORS["dark"])
    for i, group in enumerate(groups):
        values = biological.loc[biological["Number of cells"] == group, "copies_per_cell"].to_numpy()
        ax.scatter(x[i] + np.linspace(-0.05, 0.05, len(values)), values, s=10, color=COLORS["dark"], zorder=3)
    ax.set(xticks=x, xticklabels=[f"{int(v):,}" for v in groups], xlabel="Number of cells used for analysis", ylabel=r"$sxtA4$ copy number per cell")
    panel_label(ax, "c")
    save_figure(fig, output_dir, "Extended_Data_Figure1c")


def validate_ed1g_lane_map(source: Path) -> None:
    lane_map = pd.read_excel(source, sheet_name="ED_Fig.1g", header=3)
    required = {"crRNA", "Lane", "Reaction medium / control", "RNP present", "Qualitative observation", "Displayed source file"}
    if not required.issubset(lane_map.columns) or len(lane_map) != 9:
        raise ValueError("ED Fig. 1g lane map is incomplete.")
    expected = {"crRNA-site1": 4, "crRNA-site2": 5}
    if lane_map.groupby("crRNA").size().to_dict() != expected:
        raise ValueError("ED Fig. 1g must contain four crRNA-site1 lanes and five crRNA-site2 lanes.")


def plot_ed1d_f(source: Path, output_dir: Path) -> None:
    """Recreate ED Fig. 1d–f using the palette and geometry of the formal artwork."""
    accessibility = pd.read_excel(source, sheet_name="ED_Fig.1d")
    opening = pd.read_excel(source, sheet_name="ED_Fig.1e")
    position = pd.read_excel(source, sheet_name="ED_Fig.1f")
    guide_order = [
        "Site1_A4_C19",
        "Site1_A4_U19",
        "Site1_C4_C19",
        "Site1_C4_U19",
        "Site2",
    ]
    guide_labels = [
        "Site1 A4/C19",
        "Site1 A4/U19",
        "Site1 C4/C19",
        "Site1 C4/U19",
        "Site2",
    ]
    short_labels = ["S1\nA4/C19", "S1\nA4/U19", "S1\nC4/C19", "S1\nC4/U19", "Site2"]
    if set(accessibility["guide_id"]) != set(guide_order) or set(opening["guide_id"]) != set(guide_order):
        raise ValueError("ED Fig. 1d–e requires all four crRNA-site1 expansions and crRNA-site2.")
    expected_positions = {guide: 42 for guide in guide_order}
    if set(position["guide_id"]) != set(guide_order) or position.groupby("guide_id").size().to_dict() != expected_positions:
        raise ValueError("ED Fig. 1f requires 42 positions for each of five mature-crRNA models.")

    fig = plt.figure(figsize=(7.2, 6.2))
    grid = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.0, 1.0, 1.18],
        left=0.09,
        right=0.96,
        bottom=0.09,
        top=0.96,
        hspace=0.58,
    )
    x = np.arange(len(guide_order))
    width = 0.24
    full_color = "#1D9383"
    seed5_color = "#C23892"
    seed8_color = "#7B5AB6"

    axis_d = fig.add_subplot(grid[0, 0])
    ordered = accessibility.set_index("guide_id").reindex(guide_order)
    for offset, column, label, color in (
        (-width, "full_spacer_aup", "full spacer", full_color),
        (0, "seed_1_5_aup", "seed 1–5", seed5_color),
        (width, "seed_1_8_aup", "seed 1–8", seed8_color),
    ):
        axis_d.bar(
            x + offset,
            ordered[column],
            width,
            color=color,
            edgecolor="white",
            linewidth=0.4,
            label=label,
        )
    axis_d.set(
        xticks=x,
        xticklabels=short_labels,
        ylabel="Average unpaired probability (AUP)",
        ylim=(0, 1.02),
        xlabel="crRNA",
    )
    axis_d.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.62, 1.02))
    panel_label(axis_d, "d")

    axis_e = fig.add_subplot(grid[1, 0])
    ordered = opening.set_index("guide_id").reindex(guide_order)
    for offset, column, label, color in (
        (-width, "seed_1_5_opening_energy_kcal_mol", "seed 1–5", seed5_color),
        (0, "seed_1_8_opening_energy_kcal_mol", "seed 1–8", seed8_color),
        (width, "full_spacer_opening_energy_kcal_mol", "full", full_color),
    ):
        axis_e.bar(
            x + offset,
            ordered[column],
            width,
            color=color,
            edgecolor="white",
            linewidth=0.4,
            label=label,
        )
    axis_e.set(
        xticks=x,
        xticklabels=short_labels,
        ylabel=r"Opening penalty, $\Delta G$ (kcal mol$^{-1}$)",
        xlabel="crRNA",
    )
    axis_e.legend(ncol=1, loc="upper left", bbox_to_anchor=(0.00, 1.01))
    panel_label(axis_e, "e")

    axis_f = fig.add_subplot(grid[2, 0])
    matrix = (
        position.pivot(index="guide_id", columns="position", values="unpaired_probability")
        .reindex(guide_order)
    )
    color_map = LinearSegmentedColormap.from_list(
        "ed1_accessibility", ["#F6FAFB", "#9BD1C9", "#0B6E66"]
    )
    image = axis_f.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        cmap=color_map,
        vmin=0,
        vmax=1,
        interpolation="nearest",
        rasterized=True,
    )
    axis_f.axvline(20.5, color="white", lw=1.0)
    axis_f.axvline(28.5, color=seed8_color, lw=0.75, ls="--")
    axis_f.set(
        yticks=np.arange(5),
        yticklabels=guide_labels,
        xticks=np.arange(42),
        xticklabels=np.arange(1, 43),
        xlabel="Mature-model position",
    )
    axis_f.tick_params(axis="x", labelsize=4.6)
    axis_f.text(7.5, -1.03, "canonical handle", color="#3D87BE", fontsize=5.4, fontweight="bold")
    axis_f.text(23.0, -1.03, "seed 1–8", color=seed8_color, fontsize=5.4, fontweight="bold")
    axis_f.text(32.0, -1.03, "spacer remainder", color=full_color, fontsize=5.4, fontweight="bold")
    colorbar = fig.colorbar(image, ax=axis_f, fraction=0.02, pad=0.02)
    colorbar.set_label("Unpaired probability")
    panel_label(axis_f, "f")

    save_figure(fig, output_dir, "Extended_Data_Figure1d_f")


ED2_BLUE_CMAP = LinearSegmentedColormap.from_list(
    "ed2_blue",
    ["#F7FAFC", "#A3C9DC", "#6599C0", "#2F6395", "#113B65"],
)
ED2_ORANGE_CMAP = LinearSegmentedColormap.from_list(
    "ed2_orange",
    ["#F7F4F0", "#F0C9A8", "#D78350", "#7F3025"],
)


def ed2_matrix(
    data: pd.DataFrame,
    site: str,
    event_type: str,
    category_column: str,
    samples: list[str],
    categories: list[int] | None = None,
) -> pd.DataFrame:
    subset = data[(data["site"] == site) & (data["type"] == event_type)]
    if subset.duplicated(["sample", category_column]).any():
        raise ValueError(f"Duplicate ED Fig. 2 values for {site} {event_type} {category_column}.")
    matrix = subset.pivot(
        index="sample",
        columns=category_column,
        values="percent_of_sample_mutant_reads",
    )
    if categories is None:
        categories = sorted(int(value) for value in matrix.columns)
    missing_categories = sorted(set(categories) - set(matrix.columns))
    if missing_categories:
        raise ValueError(f"Missing ED Fig. 2 categories: {missing_categories}")
    return matrix.reindex(index=samples, columns=categories).fillna(0)


def draw_ed2_heatmap(
    ax: plt.Axes,
    matrix: pd.DataFrame,
    cmap,
    vmax: float,
    show_samples: bool,
    y_label: str | None = None,
):
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        interpolation="nearest",
        rasterized=True,
    )
    ax.set_xticks(np.arange(len(matrix.columns)), [str(int(value)) for value in matrix.columns])
    ax.tick_params(axis="x", labelrotation=90, labelsize=5.2, length=1.8, width=0.55, pad=1.2)
    if show_samples:
        ax.set_yticks(np.arange(len(matrix.index)), matrix.index, fontsize=5.8)
        ax.tick_params(axis="y", length=1.8, width=0.55, pad=2)
        if y_label:
            ax.set_ylabel(y_label, fontsize=6.5, labelpad=4)
    else:
        ax.set_yticks([])
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    return image


def ed2_scale_maximum(matrix: pd.DataFrame) -> int:
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("ED Fig. 2 heatmap data contain non-finite values.")
    if (values < 0).any() or (values > 100).any():
        raise ValueError("ED Fig. 2 heatmap percentages must be between 0 and 100.")
    maximum = float(values.max())
    if maximum <= 0:
        raise ValueError("ED Fig. 2 heatmap data must contain at least one positive value.")
    return int(np.ceil(maximum))


def ed2_colorbar_ticks(maximum: int) -> tuple[list[int], list[str]]:
    if maximum <= 5:
        ticks = list(range(maximum + 1))
    elif maximum <= 20:
        ticks = list(range(0, maximum, 5)) + [maximum]
    elif maximum <= 40:
        ticks = list(range(0, maximum, 10)) + [maximum]
    else:
        ticks = list(range(0, maximum, 25)) + [maximum]
    labels = [str(value) for value in ticks]
    return ticks, labels


def plot_ed2(source: Path, output_dir: Path) -> None:
    lengths = pd.read_excel(source, sheet_name="ED_Fig.2a")
    positions = pd.read_excel(source, sheet_name="ED_Fig.2b")
    required_length_columns = {"site", "sample", "type", "length_bp", "percent_of_sample_mutant_reads"}
    required_position_columns = {"site", "sample", "type", "position", "percent_of_sample_mutant_reads"}
    if not required_length_columns.issubset(lengths.columns) or not required_position_columns.issubset(positions.columns):
        raise ValueError("ED Fig. 2 source-data columns are incomplete.")
    if lengths[list(required_length_columns)].isna().any().any() or positions[list(required_position_columns)].isna().any().any():
        raise ValueError("ED Fig. 2 source data contain missing values in required columns.")

    sites = ("site1", "site2")
    sample_order = {
        site: sorted(
            set(lengths.loc[lengths["site"] == site, "sample"])
            | set(positions.loc[positions["site"] == site, "sample"])
        )
        for site in sites
    }
    if any(len(samples) != 15 for samples in sample_order.values()):
        raise ValueError("ED Fig. 2 must contain 15 samples for each crRNA site.")

    site1_insertion_positions = [128, 131, 132, 133, 134, 136, 137, 139, 140, 142]
    matrices = {}
    for site in sites:
        samples = sample_order[site]
        matrices[("a", site, "del")] = ed2_matrix(lengths, site, "del", "length_bp", samples)
        matrices[("a", site, "ins")] = ed2_matrix(lengths, site, "ins", "length_bp", samples)
        matrices[("b", site, "del")] = ed2_matrix(positions, site, "del", "position", samples)
        position_categories = site1_insertion_positions if site == "site1" else None
        matrices[("b", site, "ins")] = ed2_matrix(
            positions,
            site,
            "ins",
            "position",
            samples,
            position_categories,
        )

    expected_shapes = {
        ("a", "site1", "del"): (15, 46),
        ("a", "site1", "ins"): (15, 13),
        ("a", "site2", "del"): (15, 27),
        ("a", "site2", "ins"): (15, 5),
        ("b", "site1", "del"): (15, 47),
        ("b", "site1", "ins"): (15, 10),
        ("b", "site2", "del"): (15, 35),
        ("b", "site2", "ins"): (15, 7),
    }
    observed_shapes = {key: matrix.shape for key, matrix in matrices.items()}
    if observed_shapes != expected_shapes:
        raise ValueError(f"Unexpected ED Fig. 2 matrix shapes: {observed_shapes}")

    fig = plt.figure(figsize=(7.2, 8.35), facecolor="white")
    axes_layout = {
        ("a", "site1", "del"): [0.120, 0.775, 0.400, 0.209],
        ("a", "site1", "ins"): [0.570, 0.775, 0.180, 0.209],
        ("a", "site2", "del"): [0.120, 0.538, 0.400, 0.209],
        ("a", "site2", "ins"): [0.570, 0.538, 0.180, 0.209],
        ("b", "site1", "del"): [0.115, 0.276, 0.470, 0.191],
        ("b", "site1", "ins"): [0.645, 0.276, 0.205, 0.191],
        ("b", "site2", "del"): [0.115, 0.060, 0.470, 0.191],
        ("b", "site2", "ins"): [0.645, 0.060, 0.205, 0.191],
    }
    colorbar_layout = {
        ("a", "site1", "del"): [0.530, 0.792, 0.009, 0.175],
        ("a", "site1", "ins"): [0.760, 0.792, 0.009, 0.175],
        ("a", "site2", "del"): [0.530, 0.555, 0.009, 0.175],
        ("a", "site2", "ins"): [0.760, 0.555, 0.009, 0.175],
        ("b", "site1", "del"): [0.595, 0.292, 0.009, 0.160],
        ("b", "site1", "ins"): [0.860, 0.292, 0.009, 0.160],
        ("b", "site2", "del"): [0.595, 0.075, 0.009, 0.160],
        ("b", "site2", "ins"): [0.860, 0.075, 0.009, 0.160],
    }
    colour_scales = {
        (panel, event_type): max(
            ed2_scale_maximum(matrix)
            for (matrix_panel, _, matrix_event_type), matrix in matrices.items()
            if matrix_panel == panel and matrix_event_type == event_type
        )
        for panel in ("a", "b")
        for event_type in ("del", "ins")
    }
    images = {}
    for panel in ("a", "b"):
        for site in sites:
            for event_type in ("del", "ins"):
                ax = fig.add_axes(axes_layout[(panel, site, event_type)])
                cmap = ED2_BLUE_CMAP if event_type == "del" else ED2_ORANGE_CMAP
                site_number = site[-1]
                plural = "Samples" if panel == "a" else "Sample"
                key = (panel, site, event_type)
                images[key] = draw_ed2_heatmap(
                    ax,
                    matrices[key],
                    cmap,
                    colour_scales[(panel, event_type)],
                    show_samples=event_type == "del",
                    y_label=f"{plural} for crRNA-site{site_number}" if event_type == "del" else None,
                )
                if site == "site2":
                    if panel == "a":
                        label = "Total deleted bases per mutant sequence (bp)" if event_type == "del" else "Total inserted bases per mutant sequence (bp)"
                    else:
                        label = "Mutant-sequence deletion spectrum by reference sequence position" if event_type == "del" else "Mutant-sequence insertion spectrum\nby reference sequence position"
                    ax.set_xlabel(label, fontsize=6.5, labelpad=4)

    for key, position in colorbar_layout.items():
        ticks, tick_labels = ed2_colorbar_ticks(colour_scales[(key[0], key[2])])
        cax = fig.add_axes(position)
        colorbar = fig.colorbar(images[key], cax=cax, ticks=ticks)
        colorbar.ax.set_yticklabels(tick_labels)
        colorbar.ax.tick_params(labelsize=5.1, length=1.8, width=0.55, pad=1.6)
        colorbar.outline.set_linewidth(0.6)

    for y in (0.879, 0.642, 0.372, 0.155):
        fig.text(
            0.913,
            y,
            "Mutant sequence proportion (%)",
            rotation=90,
            ha="center",
            va="center",
            fontsize=6.5,
        )

    fig.text(0.012, 0.995, "a", ha="left", va="top", fontsize=9, fontweight="bold")
    fig.text(0.012, 0.505, "b", ha="left", va="top", fontsize=9, fontweight="bold")
    save_figure(fig, output_dir, "Extended_Data_Figure2")


ED6_SAMPLE_ORDER = [
    f"{condition}_{day}d_{replicate}"
    for day in (5, 9, 13, 15)
    for condition in ("WT", "M")
    for replicate in (1, 2, 3)
]

ED6_ANNOTATION_COLORS = {
    "PFAM": "#2C948D",
    "eggNOG": "#427796",
    "NR": "#EB9D5F",
    "Swissprot": "#DE6C52",
    "KEGG": "#897F7D",
    "GO": "#64498C",
}


def plot_ed6(source: Path, output_dir: Path) -> None:
    pca = pd.read_excel(source, sheet_name="ED_Fig.6a")
    corr = pd.read_excel(source, sheet_name="ED_Fig.6b", index_col=0)
    missing_samples = sorted(set(ED6_SAMPLE_ORDER) - set(corr.index) | set(ED6_SAMPLE_ORDER) - set(corr.columns))
    if missing_samples:
        raise ValueError(f"Missing ED Fig. 6b samples: {missing_samples}")
    corr = corr.reindex(index=ED6_SAMPLE_ORDER, columns=ED6_SAMPLE_ORDER)
    annot = pd.read_excel(source, sheet_name="ED_Fig.6c")
    support = pd.read_excel(source, sheet_name="ED_Fig.6d")
    intersections = pd.read_excel(source, sheet_name="ED_Fig.6e")

    fig = plt.figure(figsize=(7.2, 7.2), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.25])

    ax = fig.add_subplot(grid[0, 0])
    palette = plt.cm.tab10(np.linspace(0, 1, pca.condition.nunique()))
    for color, (condition, sub) in zip(palette, pca.groupby("condition")):
        ax.scatter(sub.PC1, sub.PC2, s=18, color=color, label=condition)
    ax.axhline(0, color="#D0D0D0", lw=0.6)
    ax.axvline(0, color="#D0D0D0", lw=0.6)
    ax.set(xlabel="PC1 (66%)", ylabel="PC2 (20%)")
    ax.legend(ncol=2, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=5.5)
    panel_label(ax, "a")

    ax = fig.add_subplot(grid[0, 1])
    image = ax.imshow(corr.to_numpy(dtype=float), vmin=0.9, vmax=1, cmap="viridis", interpolation="nearest")
    ax.set_xticks(np.arange(len(corr.columns)), corr.columns, rotation=90, fontsize=4.2)
    ax.set_yticks(np.arange(len(corr.index)), corr.index, fontsize=4.2)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    colorbar.ax.set_title("r value", fontsize=6.5, pad=4)
    ax.set_xlabel("Pearson correlation", labelpad=5)
    panel_label(ax, "b")

    ax = fig.add_subplot(grid[1, 0])
    annot = annot.sort_values("percent")
    bar_colors = [ED6_ANNOTATION_COLORS[source] for source in annot.source]
    ax.barh(annot.source, annot.percent, color=bar_colors)
    for y, (_, row) in enumerate(annot.iterrows()):
        ax.text(row.percent + 0.3, y, f"{int(row['count']):,} ({row.percent:.1f}%)", va="center", fontsize=5.5)
    ax.set_xlabel("Genes with annotation (%)")
    panel_label(ax, "c")

    ax = fig.add_subplot(grid[1, 1])
    ax.bar(support.annotation_source_count.astype(str), support.percent, color="#659B87")
    for i, row in support.iterrows():
        ax.text(i, row.percent + 1.2, f"{int(row['count']):,}\n{row.percent:.1f}%", ha="center", fontsize=5.2)
    ax.set(xlabel="Number of annotation sources per gene", ylabel="Genes (%)")
    panel_label(ax, "d")

    ax = fig.add_subplot(grid[2, :])
    x = np.arange(len(intersections))
    ax.bar(x, intersections["Gene count"], color=COLORS["navy"])
    ax.set_ylabel("Genes")
    ax.set_xticks(x, [str(v) for v in intersections["Gene count"]], rotation=90, fontsize=5)
    ax2 = ax.inset_axes([0, -0.55, 1, 0.42])
    membership = intersections[["5 d", "9 d", "13 d", "15 d"]].to_numpy(dtype=int).T
    for col in range(membership.shape[1]):
        present = np.where(membership[:, col] == 1)[0]
        if len(present):
            ax2.plot([col, col], [present.min(), present.max()], color=COLORS["dark"], lw=0.7)
        ax2.scatter(np.full(4, col), np.arange(4), s=8, facecolors=np.where(membership[:, col] == 1, COLORS["dark"], "white"), edgecolors="#B0B0B0", linewidths=0.4)
    ax2.set(yticks=np.arange(4), yticklabels=["5 d", "9 d", "13 d", "15 d"], xticks=[], xlim=(-0.5, len(x) - 0.5), ylim=(-0.5, 3.5))
    ax2.spines[["top", "right", "bottom", "left"]].set_visible(False)
    panel_label(ax, "e")
    save_figure(fig, output_dir, "Extended_Data_Figure6")


def half_circle_marker(side: str) -> MplPath:
    if side == "left":
        angles = np.linspace(np.pi / 2, 3 * np.pi / 2, 60)
    elif side == "right":
        angles = np.linspace(-np.pi / 2, np.pi / 2, 60)
    else:
        raise ValueError("side must be 'left' or 'right'")
    vertices = np.column_stack([np.cos(angles), np.sin(angles)])
    vertices = np.vstack([[0.0, 0.0], vertices, [0.0, 0.0]])
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * len(angles) + [MplPath.CLOSEPOLY]
    return MplPath(vertices, codes)


def fdr_alpha(value: float) -> float:
    if value <= 0.001:
        return 1.0
    if value <= 0.01:
        return 0.72
    if value <= 0.05:
        return 0.45
    return 0.16


def plot_ed7(source: Path, output_dir: Path) -> None:
    data = pd.read_excel(source, sheet_name="ED_Fig.7")
    times = ["5d", "9d", "13d", "15d"]
    up_color = "#CB493D"
    down_color = "#2E69A5"
    size_scale = 14.0
    left_marker = half_circle_marker("left")
    right_marker = half_circle_marker("right")

    module_labels = {
        "Protein turnover": "Protein\nturnover",
        "Protein folding and translation": "Protein folding\nand translation",
        "Lipid activation and synthase": "Lipid activation\nand synthase",
        "Fatty-acid biosynthesis": "Fatty-acid\nbiosynthesis",
    }
    row_metadata = (
        data[["order", "module", "displayLabel"]]
        .drop_duplicates()
        .sort_values("order")
    )
    display_rows = [
        (int(row.order), str(row.displayLabel), module_labels[str(row.module)])
        for row in row_metadata.itertuples(index=False)
    ]
    if len(display_rows) != 20 or set(data.groupby("order").size()) != {4}:
        raise ValueError("ED_Fig.7 must contain 20 GO terms with four time points per term.")

    fig = plt.figure(figsize=(7.2, 6.5), facecolor="white")
    ax = fig.add_axes([0.53, 0.15, 0.20, 0.80])
    legend_ax = fig.add_axes([0.76, 0.20, 0.22, 0.73])

    for xi in range(4):
        ax.axvline(xi, color="#D9D9D9", lw=0.45, zorder=0)

    for yi, (source_order, _, _) in enumerate(display_rows):
        sub = data[data["order"] == source_order].set_index("timepoint")
        for xi, time in enumerate(times):
            if time not in sub.index:
                continue
            row = sub.loc[time]
            alpha = fdr_alpha(float(row.FDR))
            up = float(row.up)
            down = float(row.down)
            if up > 0:
                ax.scatter(xi, yi, s=up * size_scale, marker=left_marker, color=up_color, alpha=alpha, linewidths=0, zorder=3)
            if down > 0:
                ax.scatter(xi, yi, s=down * size_scale, marker=right_marker, color=down_color, alpha=alpha, linewidths=0, zorder=3)

    labels = [row[1] for row in display_rows]
    ax.set_xlim(-0.45, 3.45)
    ax.set_ylim(len(display_rows) - 0.35, -0.65)
    ax.set_xticks(np.arange(4), ["5", "9", "13", "15"])
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.set_xlabel("Time (d)", labelpad=3)
    ax.tick_params(axis="x", length=3, width=0.7, pad=2, labelsize=7.5)
    ax.tick_params(axis="y", length=0, pad=7, labelsize=7.2)
    for tick in ax.get_yticklabels():
        tick.set_horizontalalignment("right")
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)

    axes_box = ax.get_position()
    group_centres = [2, 7, 12, 17]
    group_labels = [
        "Protein\nturnover",
        "Protein folding\nand translation",
        "Lipid activation\nand synthase",
        "Fatty-acid\nbiosynthesis",
    ]
    for centre, label in zip(group_centres, group_labels):
        fy = axes_box.y1 - (centre + 0.5) / len(display_rows) * axes_box.height
        fig.text(0.047, fy, label, rotation=90, ha="center", va="center", fontsize=7.5, fontweight="bold", linespacing=1.05)
    fig.add_artist(Line2D([0.006, 0.006], [axes_box.y0, axes_box.y1], transform=fig.transFigure, color="#B9D9D5", lw=3.0))
    for boundary in (5, 10, 15):
        fy = axes_box.y1 - boundary / len(display_rows) * axes_box.height
        fig.add_artist(Line2D([0.006, axes_box.x1], [fy, fy], transform=fig.transFigure, color="#D9D9D9", lw=0.45, zorder=0))

    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.axis("off")
    legend_ax.text(0.00, 0.98, "Direction", ha="left", va="top", fontsize=8, fontweight="bold")
    legend_ax.scatter(0.08, 0.91, s=9 * size_scale, marker=left_marker, color=up_color, linewidths=0)
    legend_ax.scatter(0.08, 0.91, s=9 * size_scale, marker=right_marker, color=down_color, linewidths=0)
    legend_ax.text(0.18, 0.91, "Up / Down", ha="left", va="center", fontsize=7.5)

    legend_ax.text(0.00, 0.81, "Transparency (FDR)", ha="left", va="top", fontsize=8, fontweight="bold")
    fdr_items = [("≤ 0.001", 1.0), ("0.001–0.01", 0.72), ("0.01–0.05", 0.45), ("> 0.05", 0.16)]
    for index, (label, alpha) in enumerate(fdr_items):
        y = 0.74 - index * 0.075
        legend_ax.scatter(0.08, y, s=32, color="#5E6B75", alpha=alpha, linewidths=0)
        legend_ax.text(0.18, y, label, ha="left", va="center", fontsize=7.2)

    legend_ax.text(0.00, 0.43, "DEG count", ha="left", va="top", fontsize=8, fontweight="bold")
    for index, count in enumerate((5, 10, 15)):
        y = 0.35 - index * 0.105
        legend_ax.scatter(0.08, y, s=count * size_scale, facecolors="none", edgecolors="#303030", linewidths=0.75)
        legend_ax.text(0.18, y, str(count), ha="left", va="center", fontsize=7.2)

    fig.text(0.49, 0.045, "GO enrichment terms linked to growth in the mutant strain", ha="center", va="center", fontsize=8)
    save_figure(fig, output_dir, "Extended_Data_Figure7")


def plot(source_dir: Path, output_dir: Path) -> None:
    configure_style()
    source = source_dir / "Source_data_Extended_data_fig.xlsx"
    plot_ed1b(source, output_dir)
    plot_ed1c(source, output_dir)
    plot_ed1d_f(source, output_dir)
    validate_ed1g_lane_map(source)
    plot_ed2(source, output_dir)
    plot_ed6(source, output_dir)
    plot_ed7(source, output_dir)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
