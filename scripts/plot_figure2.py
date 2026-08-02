from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import COLORS, add_bracket, configure_style, panel_label, save_figure, significance_label


def grouped_recovery_panel(
    ax: plt.Axes,
    data: pd.DataFrame,
    group_column: str,
    groups: list[float],
    tick_labels: list[str],
    xlabel: str,
    panel: str,
) -> None:
    required = {group_column, "germination_rate", "viability_rate"}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Figure 2{panel} is missing columns: {', '.join(missing)}")
    x = np.arange(len(groups))
    width = 0.34
    for offset, outcome, label, color in (
        (-width / 2, "germination_rate", "Germination", COLORS["blue"]),
        (width / 2, "viability_rate", "Viable-cell recovery", COLORS["teal"]),
    ):
        group_values = [
            data.loc[np.isclose(data[group_column], group), outcome].dropna().to_numpy(dtype=float)
            for group in groups
        ]
        means = [values.mean() for values in group_values]
        sds = [values.std(ddof=1) for values in group_values]
        ax.bar(x + offset, means, width, yerr=sds, capsize=2, color=color, edgecolor=COLORS["dark"], label=label)
        for i, values in enumerate(group_values):
            ax.scatter(np.full(values.size, x[i] + offset), values, s=10, c=COLORS["dark"], zorder=3)
        p_value = stats.ttest_ind(group_values[0], group_values[1], equal_var=False).pvalue
        add_bracket(ax, x[0] + offset, x[1] + offset, 87 if outcome == "viability_rate" else 78, significance_label(p_value, tiers=1), height=2)
    ax.set(xticks=x, xticklabels=tick_labels, xlabel=xlabel, ylabel="Percentage (%)", ylim=(0, 100))
    panel_label(ax, panel)


def plot(source_dir: Path, output_dir: Path) -> None:
    configure_style()
    source = source_dir / "Source_data_fig2.xlsx"
    tip = pd.read_excel(source, sheet_name="Fig.2d")
    uvrd = pd.read_excel(source, sheet_name="Fig.2f")

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), constrained_layout=True)
    grouped_recovery_panel(axes[0], tip, "tip_size_nm", [500, 50], ["500", "50"], "Micropipette tip size (nm)", "d")
    grouped_recovery_panel(axes[1], uvrd, "uvrd_uM", [0.003, 0.012], ["0.003", "0.012"], "UvrD concentration (µM)", "f")
    axes[0].legend(loc="upper left")
    save_figure(fig, output_dir, "Figure2_data_panels")


def main() -> None:
    parser = argparse.ArgumentParser(description="Recreate data-driven panels for current Figure 2.")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
