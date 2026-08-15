from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
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
    save_figure(fig, output_dir, "Supplementary_Figure_S1b")


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
    save_figure(fig, output_dir, "Supplementary_Figure_S1c")


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

    save_figure(fig, output_dir, "Supplementary_Figure_S1d_f")


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
    save_figure(fig, output_dir, "Supplementary_Figure_S2")


def plot(source_dir: Path, output_dir: Path) -> None:
    configure_style()
    supplementary_1 = source_dir / "Source_Data_Supplementary_Figure_S1.xlsx"
    supplementary_2 = source_dir / "Source_Data_Supplementary_Figure_S2.xlsx"
    plot_ed1b(supplementary_1, output_dir)
    plot_ed1c(supplementary_1, output_dir)
    plot_ed1d_f(supplementary_1, output_dir)
    plot_ed2(supplementary_2, output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recreate data-driven panels for Supplementary Figures S1 and S2."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
