from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import COLORS, add_bracket, configure_style, panel_label, save_figure, significance_label


WT = "Wild-type strain"
EDITED = "sxtA4-edited strain"


def parse_density(value: object) -> float:
    text = str(value).replace("×", "x")
    match = re.search(r"([0-9.]+)\s*x\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹]+)", text)
    exponent = int(match.group(2).translate(str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")))
    return float(match.group(1)) * 10**exponent


def plot(source_dir: Path, output_dir: Path) -> None:
    configure_style()
    source = source_dir / "Source_data_fig3.xlsx"
    pigment = pd.read_excel(source, sheet_name="Fig.3b", header=2)
    growth = pd.read_excel(source, sheet_name="Fig.3c", header=1)
    qpcr = pd.read_excel(source, sheet_name="Fig.3d", header=2)
    toxin = pd.read_excel(source, sheet_name="Fig. 3e", header=2)
    toxin = toxin[toxin.Samples.notna()].copy()

    fig = plt.figure(figsize=(7.2, 7.3), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 1])

    ax = fig.add_subplot(grid[0, 0])
    pigment_cols = ["Bcar", "Chlide a", "MgDVP", "Chl c3", "Chl c2", "Viol", "Diad", "Peri", "Chl a"]
    pigment_colors = ["#D94D7C", "#E99AB3", "#E68C66", "#F0A14A", "#F6C85F", "#FBE6B7", "#93D3D3", "#7FA3C5", "#9BC59B"]
    pigment["day"] = pigment.Samples.str.extract(r"-D(\d+)-")[0].astype(int)
    days = sorted(pigment.day.unique())
    x = np.arange(len(days))
    for direction, strain in ((1, WT), (-1, EDITED)):
        bottom = np.zeros(len(days))
        for column, color in zip(pigment_cols, pigment_colors):
            values = np.array([pigment.loc[(pigment.Strain == strain) & (pigment.day == day), column].mean() for day in days])
            ax.bar(x, direction * values, bottom=direction * bottom, color=color, width=0.78, linewidth=0, label=column if direction == 1 else None)
            bottom += values
    ax.axhline(0, color=COLORS["dark"], lw=0.7)
    ax.set(xticks=x, xticklabels=days, xlabel="Time (days)", ylabel="Culture pigment content (µg L−1)")
    ax.legend(
        ncol=3,
        bbox_to_anchor=(0.5, 1.02),
        loc="lower center",
        fontsize=5.4,
        columnspacing=0.8,
        handletextpad=0.4,
    )
    panel_label(ax, "b")

    ax2 = ax.twinx()
    all_pigments = list(pigment.columns[2:27])
    pigment["total"] = pigment[all_pigments].sum(axis=1)
    pigment["pg_cell"] = pigment.total * 1e6 / pigment["Cell density (cells/L)"]
    for strain, color, marker, sign in ((WT, COLORS["dark"], "o", 1), (EDITED, COLORS["edited"], "o", -1)):
        means = np.array([pigment.loc[(pigment.Strain == strain) & (pigment.day == day), "pg_cell"].mean() for day in days])
        sds = np.array([pigment.loc[(pigment.Strain == strain) & (pigment.day == day), "pg_cell"].std(ddof=1) for day in days])
        ax2.errorbar(x, sign * means, yerr=sds, color=color, marker=marker, ms=2.8, capsize=1.5, lw=0.8)
    ax2.set_ylabel("Cellular pigment content (pg cell−1)")
    lim = max(abs(v) for v in ax2.get_ylim())
    ax2.set_ylim(-lim, lim)

    ax = fig.add_subplot(grid[0, 1])
    growth["density"] = growth["Cell density (cells/L)"].map(parse_density)
    growth["ln_density"] = np.log(growth.density / 1000)
    growth = growth[growth.Days.between(5, 23)]
    daily = growth.groupby(["Strain", "Days"]).ln_density.agg(["mean", "std"]).reset_index()
    for strain, color in ((WT, COLORS["wt"]), (EDITED, COLORS["edited"])):
        sub = daily[daily.Strain == strain]
        ax.errorbar(sub.Days, sub["mean"], yerr=sub["std"], color=color, marker="o", ms=2.8, capsize=1.5, lw=0.8, label=strain)
        fit = stats.linregress(sub.Days, sub["mean"])
        xx = np.linspace(5, 23, 100)
        ax.plot(xx, fit.intercept + fit.slope * xx, color=color, ls="--", lw=0.8)
        ax.text(
            0.98,
            0.92 if strain == EDITED else 0.18,
            f"y = {fit.slope:.2f}x + {fit.intercept:.2f}\n$R^2$ = {fit.rvalue**2:.2f}",
            color=color,
            transform=ax.transAxes,
            ha="right",
            va="top",
        )
    ax.set(xlabel="Time (days)", ylabel="ln(cell density per mL)")
    ax.legend(loc="upper left")
    panel_label(ax, "c")

    ax = fig.add_subplot(grid[1, :])
    qbio = qpcr.groupby(["Samples", "Days", "Cell density (cells/ml)"], as_index=False).Ct.mean()
    qbio["Strain"] = np.where(qbio.Samples.str.contains("91"), WT, EDITED)
    qbio["transcripts"] = 10 ** ((39.67 - qbio.Ct) / 3.59) / qbio["Cell density (cells/ml)"]
    qbio["log10_transcripts"] = np.log10(qbio.transcripts)
    days = [5, 9, 13, 15]
    x = np.arange(len(days))
    width = 0.32
    for offset, strain, color in ((-width / 2, WT, COLORS["wt"]), (width / 2, EDITED, COLORS["edited"])):
        means = [qbio.loc[(qbio.Strain == strain) & (qbio.Days == day), "transcripts"].mean() for day in days]
        sds = [qbio.loc[(qbio.Strain == strain) & (qbio.Days == day), "transcripts"].std(ddof=1) for day in days]
        ax.bar(x + offset, means, width, yerr=sds, capsize=2, color=color, edgecolor=COLORS["dark"], label=strain)
        for i, day in enumerate(days):
            vals = qbio.loc[(qbio.Strain == strain) & (qbio.Days == day), "transcripts"].to_numpy()
            jitter = np.linspace(-0.04, 0.04, len(vals))
            ax.scatter(x[i] + offset + jitter, vals, s=8, color=COLORS["dark"], zorder=3)
    ymax = max(ax.get_ylim()[1], 310)
    ax.set_ylim(0, ymax)
    for i, day in enumerate(days):
        sub = qbio[qbio.Days == day]
        p = stats.ttest_ind(sub.loc[sub.Strain == WT, "log10_transcripts"], sub.loc[sub.Strain == EDITED, "log10_transcripts"], equal_var=True).pvalue
        y = max(sub.transcripts.max() * 1.12, 35)
        add_bracket(ax, x[i] - width / 2, x[i] + width / 2, y, significance_label(p, tiers=2), height=max(ymax * 0.015, 2))
    ax.set(xticks=x, xticklabels=days, xlabel="Time (days)", ylabel="sxtA4 transcripts per cell")
    ax.legend(loc="upper right")
    panel_label(ax, "d")

    ax = fig.add_subplot(grid[2, 0])
    congeners = ["GTX1", "GTX4", "GTX2", "GTX3", "STX", "NEO"]
    colors = ["#F3C64E", "#1EB3D7", "#9AB7CE", "#79221E", "#F0A188", "#C8102E"]
    for xpos, strain in enumerate((WT, EDITED)):
        totals = toxin.loc[toxin.Strain == strain, congeners].sum()
        ax.pie(totals, colors=colors, radius=0.42, center=(xpos, 0), wedgeprops={"edgecolor": "white", "linewidth": 0.4})
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(-0.55, 0.55)
    ax.set_aspect("equal")
    ax.set_xticks([0, 1], ["Wild type", "sxtA4 edited"])
    ax.legend(congeners, ncol=2, bbox_to_anchor=(0.5, 1.02), loc="lower center")
    panel_label(ax, "e1")

    ax = fig.add_subplot(grid[2, 1])
    toxin["day"] = toxin.Samples.str.extract(r"-(?:91|54)-(\d+)-")[0].astype(int)
    toxin["fmol_cell"] = (
        toxin.Total
        * toxin["LC-MS/MS extract volume (mL)"]
        * toxin["Dilution factor"]
        * 1e6
        / (toxin["Cell density (cells/ml)"] * toxin["Culture volume collected (mL)"])
    )
    days = sorted(toxin.day.unique())
    x = np.arange(len(days))
    width = 0.34
    for offset, strain, color in ((-width / 2, WT, COLORS["wt"]), (width / 2, EDITED, COLORS["edited"])):
        means = [toxin.loc[(toxin.Strain == strain) & (toxin.day == day), "fmol_cell"].mean() for day in days]
        sds = [toxin.loc[(toxin.Strain == strain) & (toxin.day == day), "fmol_cell"].std(ddof=1) for day in days]
        ax.bar(x + offset, means, width, yerr=sds, capsize=1.5, color=color, edgecolor=COLORS["dark"], label=strain)
    for i, day in enumerate(days):
        sub = toxin[toxin.day == day]
        p = stats.ttest_ind(sub.loc[sub.Strain == WT, "fmol_cell"], sub.loc[sub.Strain == EDITED, "fmol_cell"], equal_var=False).pvalue
        y = max(sub.fmol_cell.max() * 1.08, 0.04)
        add_bracket(ax, x[i] - width / 2, x[i] + width / 2, y, significance_label(p, tiers=2), height=max(ax.get_ylim()[1] * 0.012, 0.015))
    ax.set(xticks=x, xticklabels=days, xlabel="Time (days)", ylabel="Cellular toxin (fmol cell−1)")
    ax.legend(loc="upper right")
    panel_label(ax, "e2")

    save_figure(fig, output_dir, "Figure3_data_panels")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.source, args.output)


if __name__ == "__main__":
    main()
