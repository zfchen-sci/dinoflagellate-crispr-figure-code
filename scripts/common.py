from __future__ import annotations

import logging
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


logging.getLogger("fontTools.subset").setLevel(logging.ERROR)


COLORS = {
    "wt": "#E99B5B",
    "edited": "#4E8B7B",
    "blue": "#6E97B7",
    "orange": "#D08A61",
    "grey": "#B9B7AF",
    "dark": "#33383E",
    "teal": "#2F9B8F",
    "red": "#C64B46",
    "navy": "#2878A8",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 8,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "lines.linewidth": 1.0,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def save_figure(fig: mpl.figure.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "pdf"):
        fig.savefig(output_dir / f"{stem}.{suffix}", bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def significance_label(p: float | None, tiers: int = 4) -> str:
    if p is None or not np.isfinite(p):
        return ""
    if tiers >= 4 and p < 0.0001:
        return "****"
    if tiers >= 3 and p < 0.001:
        return "***"
    if tiers >= 2 and p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def add_bracket(
    ax: mpl.axes.Axes,
    x1: float,
    x2: float,
    y: float,
    label: str,
    height: float | None = None,
    linewidth: float = 0.8,
) -> None:
    ymin, ymax = ax.get_ylim()
    span = ymax - ymin
    h = height if height is not None else span * 0.025
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=linewidth, clip_on=False)
    ax.text((x1 + x2) / 2, y + h + span * 0.008, label, ha="center", va="bottom", fontsize=6.5)


def mean_sd(values) -> tuple[float, float, int]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.mean(arr)), float(np.std(arr, ddof=1)), int(arr.size)


def safe_log10(values):
    arr = np.asarray(values, dtype=float)
    return np.log10(np.where(arr > 0, arr, np.nan))


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.14, 1.06, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def format_p(p: float | None) -> str:
    if p is None or not math.isfinite(p):
        return ""
    if p < 0.0001:
        return f"{p:.2e}"
    return f"{p:.4f}"
