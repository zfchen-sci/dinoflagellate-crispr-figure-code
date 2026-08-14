from __future__ import annotations

import logging
import math
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


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
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif", "serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 8,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
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


def games_howell(groups: dict[str, object], alpha: float = 0.05) -> list[dict]:
    """Return Games-Howell pairwise comparisons for independent groups."""
    clean = {}
    for label, raw in groups.items():
        values = np.asarray(raw, dtype=float)
        values = values[np.isfinite(values)]
        if values.size < 2:
            raise ValueError("Games-Howell comparisons require at least two observations per group.")
        clean[str(label)] = values
    labels = list(clean)
    k = len(labels)
    results = []
    for i, first in enumerate(labels):
        for second in labels[i + 1 :]:
            a, b = clean[first], clean[second]
            na, nb = len(a), len(b)
            va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
            component_a, component_b = va / na, vb / nb
            se = np.sqrt(component_a + component_b)
            difference = float(np.mean(b) - np.mean(a))
            df = (component_a + component_b) ** 2 / (
                component_a**2 / (na - 1) + component_b**2 / (nb - 1)
            )
            q = np.sqrt(2) * abs(difference) / se
            p_adjusted = float(stats.studentized_range.sf(q, k, df))
            critical = float(stats.studentized_range.ppf(1 - alpha, k, df) / np.sqrt(2))
            half_width = critical * se
            results.append(
                {
                    "group_1": first,
                    "group_2": second,
                    "mean_difference": difference,
                    "df": float(df),
                    "p_adjusted": p_adjusted,
                    "ci_low": float(difference - half_width),
                    "ci_high": float(difference + half_width),
                }
            )
    return results


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
    ax.annotate(
        label,
        xy=((x1 + x2) / 2, y + h),
        xycoords="data",
        xytext=(0, 2),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=6.5,
        annotation_clip=False,
    )


def mean_sd(values) -> tuple[float, float, int]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.mean(arr)), float(np.std(arr, ddof=1)), int(arr.size)


def safe_log10(values):
    arr = np.asarray(values, dtype=float)
    return np.log10(np.where(arr > 0, arr, np.nan))


def infer_fig4_qpcr_strain(samples, wt_label: str, edited_label: str) -> np.ndarray:
    """Map the Fig. 4d culture identifiers to the two analysis groups.

    The current source workbook uses culture identifier 91 for the wild-type
    group and 54 for the edited group.  Match complete hyphen-delimited tokens
    so unrelated digits elsewhere in a sample name cannot silently change the
    assignment.
    """
    labels = np.asarray(samples, dtype=str)
    is_wt = np.array([bool(re.search(r"(?:^|-)91-", value)) for value in labels])
    is_edited = np.array([bool(re.search(r"(?:^|-)54-", value)) for value in labels])
    invalid = is_wt == is_edited
    if invalid.any():
        bad = ", ".join(sorted(set(labels[invalid])))
        raise ValueError(f"Unrecognized or ambiguous Fig. 4d sample identifiers: {bad}")
    return np.where(is_wt, wt_label, edited_label)


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(-0.14, 1.06, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def format_p(p: float | None) -> str:
    if p is None or not math.isfinite(p):
        return ""
    if p < 0.0001:
        return f"{p:.2e}"
    return f"{p:.4f}"
