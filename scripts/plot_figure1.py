from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import COLORS, configure_style, panel_label, save_figure


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
    configure_style()
    source = source_dir / "Source_data_fig1.xlsx"
    composition = pd.read_excel(source, sheet_name="Fig.1a-b_base_comp", header=2)
    metrics = pd.read_excel(source, sheet_name="Fig.1e_metrics")
    pairs = pd.read_excel(source, sheet_name="Fig.1f_pair_prob")

    validate_columns(composition, {"Position", "A", "T", "C", "G", "GC content"}, "Fig. 1a–b")
    validate_columns(metrics, {"guide_id", "guide_label", "sequence_type", "mfe_structure", "mfe_dg_kcal_mol"}, "Fig. 1e")
    validate_columns(pairs, {"guide_id", "guide_label", "i", "j", "pair_probability", "present_in_mfe"}, "Fig. 1f")

    fig = plt.figure(figsize=(7.2, 6.3), constrained_layout=True)
    outer = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.9, 1.25])

    ax = fig.add_subplot(outer[0, 0])
    ax.plot(composition["Position"], composition["GC content"], color=COLORS["navy"], lw=0.8)
    ax.axvspan(2222, 2243, color="#23D9E1", alpha=0.35, linewidth=0)
    ax.axvspan(2446, 2467, color="#23D9E1", alpha=0.35, linewidth=0)
    ax.set(xlabel=r"Position in $sxtA$ (bp)", ylabel="GC content (%)", xlim=(1, 2870), ylim=(0, 100))
    panel_label(ax, "a")

    ax = fig.add_subplot(outer[0, 1])
    nucleotide_colors = {"A": "#2A85B8", "C": "#3A9B42", "G": "#D42A7B", "T": "#E5A60A"}
    snp = composition.loc[
        composition["Position"].between(2180, 2870)
        & ((composition[["A", "T", "C", "G"]] > 0.01).sum(axis=1) > 1)
    ].copy()
    for start, end in ((2222, 2243), (2446, 2467)):
        ax.axvspan(start, end, color="#23D9E1", alpha=0.55, linewidth=0, zorder=0)
    for base in ("A", "C", "G", "T"):
        values = snp[base].where(snp[base] > 0.01)
        ax.scatter(snp["Position"], values, s=7, color=nucleotide_colors[base], label=base, linewidths=0)
    ax.set(xlabel=r"SNP position in $sxtA4$ (bp)", ylabel="Relative abundance (%)", xlim=(2180, 2870), ylim=(-3, 103))
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.17), columnspacing=0.8, handletextpad=0.3)
    panel_label(ax, "b")

    ax = fig.add_subplot(outer[1, :])
    mature = metrics.loc[metrics["sequence_type"] == "mature canonical-handle model"].copy()
    mature["guide_id"] = pd.Categorical(mature["guide_id"], GUIDE_ORDER, ordered=True)
    mature = mature.sort_values("guide_id")
    if list(mature["guide_id"].astype(str)) != GUIDE_ORDER:
        raise ValueError("Fig. 1e requires five mature-crRNA models in the expected guide order.")
    labels = mature["guide_label"].str.replace("Site", "crRNA-site", regex=False)
    y = np.arange(len(mature))
    bars = ax.barh(y, -mature["mfe_dg_kcal_mol"], color=COLORS["teal"], edgecolor=COLORS["dark"], height=0.64)
    for bar, value, structure in zip(bars, mature["mfe_dg_kcal_mol"], mature["mfe_structure"]):
        ax.text(bar.get_width() + 0.18, bar.get_y() + bar.get_height() / 2, f"{value:.1f} kcal mol$^{{-1}}$", va="center", fontsize=6.2)
        ax.text(0.12, bar.get_y() + bar.get_height() / 2, structure, va="center", ha="left", fontsize=5.4, family="monospace", color="white")
    ax.set(yticks=y, yticklabels=labels, xlabel="Magnitude of predicted MFE (kcal mol$^{-1}$)")
    ax.invert_yaxis()
    panel_label(ax, "e")

    pair_grid = outer[2, :].subgridspec(1, 5, wspace=0.12)
    for index, guide in enumerate(GUIDE_ORDER):
        ax = fig.add_subplot(pair_grid[0, index])
        subset = pairs.loc[pairs["guide_id"] == guide].copy()
        if subset.empty:
            raise ValueError(f"Fig. 1f has no base-pair probabilities for {guide}.")
        probability = subset["pair_probability"].to_numpy(dtype=float)
        colors = np.where(subset["present_in_mfe"].astype(bool), COLORS["red"], COLORS["blue"])
        ax.scatter(subset["i"], subset["j"], s=4 + 34 * probability, c=colors, alpha=0.75, linewidths=0)
        ax.plot([1, 42], [1, 42], color="#BDBDBD", lw=0.5)
        ax.axvline(21.5, color="#888888", lw=0.5, ls="--")
        ax.axhline(21.5, color="#888888", lw=0.5, ls="--")
        label = str(subset["guide_label"].iloc[0]).replace("Site", "crRNA-site")
        ax.set(xlim=(0.5, 42.5), ylim=(0.5, 42.5), aspect="equal", title=label)
        ax.set_xticks([1, 21, 42])
        ax.set_yticks([1, 21, 42] if index == 0 else [])
        if index == 0:
            ax.set_ylabel("Position j")
            panel_label(ax, "f")
        ax.set_xlabel("Position i")

    save_figure(fig, output_dir, "Figure1_data_panels")


def main() -> None:
    parser = argparse.ArgumentParser(description="Recreate data-driven panels for current Figure 1.")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
