from __future__ import annotations

import argparse
import json
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

from common import format_p, mean_sd, significance_label


WT = "Wild-type strain"
EDITED = "sxtA4-edited strain"


def t_row(
    figure: str,
    panel: str,
    comparison: str,
    outcome: str,
    a,
    b,
    test_name: str,
    equal_var: bool,
    source_file: str,
    source_sheet: str,
    correction: str = "None",
    p_adjusted: float | None = None,
    tiers: int = 4,
    notes: str = "",
) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    result = stats.ttest_ind(a, b, equal_var=equal_var)
    p_display = float(result.pvalue if p_adjusted is None else p_adjusted)
    ma, sda, na = mean_sd(a)
    mb, sdb, nb = mean_sd(b)
    if equal_var:
        df = na + nb - 2
    else:
        va = np.var(a, ddof=1) / na
        vb = np.var(b, ddof=1) / nb
        df = (va + vb) ** 2 / (va**2 / (na - 1) + vb**2 / (nb - 1))
    return {
        "figure": figure,
        "panel": panel,
        "comparison": comparison,
        "outcome": outcome,
        "group_1_mean": ma,
        "group_1_sd": sda,
        "group_1_n": na,
        "group_2_mean": mb,
        "group_2_sd": sdb,
        "group_2_n": nb,
        "test": test_name,
        "correction": correction,
        "statistic": float(result.statistic),
        "df": float(df),
        "p_raw": float(result.pvalue),
        "p_adjusted": p_adjusted,
        "p_used_for_display": p_display,
        "p_formatted": format_p(p_display),
        "significance": significance_label(p_display, tiers=tiers),
        "source_file": source_file,
        "source_sheet": source_sheet,
        "notes": notes,
    }


def add_holm(rows: list[dict], indices: list[int]) -> None:
    adjusted = multipletests([rows[i]["p_raw"] for i in indices], method="holm")[1]
    for i, p_adj in zip(indices, adjusted):
        rows[i]["p_adjusted"] = float(p_adj)
        rows[i]["p_used_for_display"] = float(p_adj)
        rows[i]["p_formatted"] = format_p(float(p_adj))
        rows[i]["significance"] = significance_label(float(p_adj), tiers=4)
        rows[i]["correction"] = "Holm"


def parse_injection_data(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Injection_data", header=21)
    df = df[df["crRNA"].isin(["Site1", "Site2"])].copy()
    for label, denominator in {
        "per_injected": "Injected cells",
        "per_germinated": "Germinated cells",
        "per_viable": "Viable cells",
    }.items():
        df[label] = 100 * df["Positive cells"] / df[denominator].replace(0, np.nan)
    return df


def parse_growth_density(value: object) -> float:
    text = str(value).replace("×", "x")
    match = re.search(r"([0-9.]+)\s*x\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹]+)", text)
    if not match:
        raise ValueError(f"Unrecognized cell-density value: {value}")
    exponent = int(match.group(2).translate(str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")))
    return float(match.group(1)) * 10**exponent


def calculate_statistics(source_dir: Path) -> dict[str, list[dict]]:
    tests: list[dict] = []
    summaries: list[dict] = []
    checks: list[dict] = []

    f1 = source_dir / "Source_data_fig1.xlsx"
    d1 = pd.read_excel(f1, sheet_name="Fig.1d")
    for outcome in ("germination_rate", "viability_rate"):
        tests.append(
            t_row(
                "Fig. 1",
                "d",
                "500-nm tip vs 50-nm tip",
                outcome,
                d1.loc[d1.tip_size_nm == 500, outcome],
                d1.loc[d1.tip_size_nm == 50, outcome],
                "Two-tailed independent-samples Student t-test",
                True,
                f1.name,
                "Fig.1d",
                tiers=1,
            )
        )
    d1f = pd.read_excel(f1, sheet_name="Fig.1f")
    for outcome in ("germination_rate", "viability_rate"):
        tests.append(
            t_row(
                "Fig. 1",
                "f",
                "0.003 µM vs 0.012 µM UvrD",
                outcome,
                d1f.loc[np.isclose(d1f.uvrd_uM, 0.003), outcome],
                d1f.loc[np.isclose(d1f.uvrd_uM, 0.012), outcome],
                "Two-tailed independent-samples Student t-test",
                True,
                f1.name,
                "Fig.1f",
                tiers=1,
            )
        )

    f2 = source_dir / "Source_data_fig2.xlsx"
    inj = parse_injection_data(f2)
    b = inj[(inj.crRNA == "Site2") & (inj["Injection strategy"] == "UvrD+Nuclei")]
    groups = [b.loc[b["RNP concentration (nM)"] == dose, "per_viable"].dropna() for dose in (5, 30, 1200)]
    omnibus = stats.f_oneway(*groups)
    tests.append(
        {
            "figure": "Fig. 2",
            "panel": "b",
            "comparison": "5 vs 30 vs 1,200 nM RNP",
            "outcome": "editing recovery per viable cell",
            "group_1_mean": None,
            "group_1_sd": None,
            "group_1_n": 3,
            "group_2_mean": None,
            "group_2_sd": None,
            "group_2_n": 3,
            "test": "One-way ANOVA",
            "correction": "Tukey HSD for pairwise comparisons",
            "statistic": float(omnibus.statistic),
            "df": "2, 6",
            "p_raw": float(omnibus.pvalue),
            "p_adjusted": None,
            "p_used_for_display": float(omnibus.pvalue),
            "p_formatted": format_p(float(omnibus.pvalue)),
            "significance": significance_label(float(omnibus.pvalue)),
            "source_file": f2.name,
            "source_sheet": "Injection_data",
            "notes": "Omnibus test; the figure displays Tukey pairwise results.",
        }
    )
    tukey_b = pairwise_tukeyhsd(b["per_viable"], b["RNP concentration (nM)"])
    for (g1, g2), diff, p_adj, ci in zip(
        combinations(tukey_b.groupsunique, 2), tukey_b.meandiffs, tukey_b.pvalues, tukey_b.confint
    ):
        tests.append(
            {
                "figure": "Fig. 2",
                "panel": "b",
                "comparison": f"{g1:g} vs {g2:g} nM RNP",
                "outcome": "editing recovery per viable cell",
                "group_1_mean": float(b.loc[b["RNP concentration (nM)"] == g1, "per_viable"].mean()),
                "group_1_sd": float(b.loc[b["RNP concentration (nM)"] == g1, "per_viable"].std(ddof=1)),
                "group_1_n": 3,
                "group_2_mean": float(b.loc[b["RNP concentration (nM)"] == g2, "per_viable"].mean()),
                "group_2_sd": float(b.loc[b["RNP concentration (nM)"] == g2, "per_viable"].std(ddof=1)),
                "group_2_n": 3,
                "test": "Tukey HSD",
                "correction": "Family-wise error rate",
                "statistic": float(diff),
                "df": None,
                "p_raw": None,
                "p_adjusted": float(p_adj),
                "p_used_for_display": float(p_adj),
                "p_formatted": format_p(float(p_adj)),
                "significance": significance_label(float(p_adj)),
                "source_file": f2.name,
                "source_sheet": "Injection_data",
                "notes": f"95% CI for mean difference: {ci[0]:.4f} to {ci[1]:.4f}.",
            }
        )

    c = inj[(inj["Injection strategy"] == "UvrD+Nuclei") & (inj["RNP concentration (nM)"] == 1200)]
    for outcome in ("per_injected", "per_germinated", "per_viable"):
        tests.append(
            t_row(
                "Fig. 2",
                "c",
                "crRNA-site1 vs crRNA-site2",
                outcome,
                c.loc[c.crRNA == "Site1", outcome],
                c.loc[c.crRNA == "Site2", outcome],
                "Two-tailed independent-samples Student t-test",
                True,
                f2.name,
                "Injection_data",
            )
        )

    edit = pd.read_excel(f2, sheet_name="Fig.2d-f")
    edit["lineage"] = edit["sample"].str.rsplit("-", n=1).str[0]
    metric = "corrected_editing_efficiency_pct"
    lineage = edit.groupby(["Target_site", "lineage"], as_index=False)[metric].mean()
    tests.append(
        t_row(
            "Fig. 2",
            "d",
            "crRNA-site1 vs crRNA-site2",
            "mean corrected editing efficiency per lineage",
            lineage.loc[lineage.Target_site == "crRNA_site1", metric],
            lineage.loc[lineage.Target_site == "crRNA_site2", metric],
            "Two-tailed Welch t-test",
            False,
            f2.name,
            "Fig.2d-f",
            notes="Lineage means are the independent observations; each mean summarizes three sequencing measurements.",
        )
    )
    for site, sub in edit.groupby("Target_site"):
        omnibus = stats.f_oneway(*(g[metric].to_numpy() for _, g in sub.groupby("lineage")))
        tests.append(
            {
                "figure": "Fig. 2",
                "panel": "e",
                "comparison": f"All lineages within {site.replace('_', '-')}",
                "outcome": "corrected editing efficiency",
                "group_1_mean": None,
                "group_1_sd": None,
                "group_1_n": 3,
                "group_2_mean": None,
                "group_2_sd": None,
                "group_2_n": 3,
                "test": "One-way ANOVA",
                "correction": "Tukey HSD for pairwise comparisons",
                "statistic": float(omnibus.statistic),
                "df": f"4, {len(sub) - 5}",
                "p_raw": float(omnibus.pvalue),
                "p_adjusted": None,
                "p_used_for_display": float(omnibus.pvalue),
                "p_formatted": format_p(float(omnibus.pvalue)),
                "significance": significance_label(float(omnibus.pvalue)),
                "source_file": f2.name,
                "source_sheet": "Fig.2d-f",
                "notes": "Omnibus test; the figure displays Tukey pairwise results.",
            }
        )
        tukey = pairwise_tukeyhsd(sub[metric], sub["lineage"])
        for (g1, g2), diff, p_adj, ci in zip(
            combinations(tukey.groupsunique, 2), tukey.meandiffs, tukey.pvalues, tukey.confint
        ):
            tests.append(
                {
                    "figure": "Fig. 2",
                    "panel": "e",
                    "comparison": f"{g1} vs {g2}",
                    "outcome": "corrected editing efficiency",
                    "group_1_mean": float(sub.loc[sub.lineage == g1, metric].mean()),
                    "group_1_sd": float(sub.loc[sub.lineage == g1, metric].std(ddof=1)),
                    "group_1_n": 3,
                    "group_2_mean": float(sub.loc[sub.lineage == g2, metric].mean()),
                    "group_2_sd": float(sub.loc[sub.lineage == g2, metric].std(ddof=1)),
                    "group_2_n": 3,
                    "test": "Tukey HSD",
                    "correction": "Family-wise error rate",
                    "statistic": float(diff),
                    "df": None,
                    "p_raw": None,
                    "p_adjusted": float(p_adj),
                    "p_used_for_display": float(p_adj),
                    "p_formatted": format_p(float(p_adj)),
                    "significance": significance_label(float(p_adj)),
                    "source_file": f2.name,
                    "source_sheet": "Fig.2d-f",
                    "notes": f"95% CI for mean difference: {ci[0]:.4f} to {ci[1]:.4f}.",
                }
            )
    checks.append(
        {
            "severity": "High",
            "location": "Fig. 2b,e statistical annotations and Methods",
            "issue": "The displayed pairwise annotations are reproduced by one-way ANOVA followed by Tukey HSD, whereas the current Methods describes only independent-samples t-tests for these panels.",
            "resolution_in_scripts": "The scripts use one-way ANOVA with Tukey HSD and report the adjusted pairwise P values used by the figures.",
            "submission_action": "Add the ANOVA and Tukey HSD procedure to the statistical analysis subsection and the relevant figure legends.",
        }
    )

    f3 = source_dir / "Source_data_fig3.xlsx"
    pigment = pd.read_excel(f3, sheet_name="Fig.3b", header=2)
    pigment_columns = list(pigment.columns[2:27])
    pigment["day"] = pigment.Samples.str.extract(r"-D(\d+)-")[0].astype(int)
    pigment["total_pigment_ug_L"] = pigment[pigment_columns].sum(axis=1)
    pigment["cellular_pigment_pg_cell"] = (
        pigment["total_pigment_ug_L"] * 1e6 / pigment["Cell density (cells/L)"]
    )
    pigment_indices = []
    for day, sub in pigment.groupby("day"):
        pigment_indices.append(len(tests))
        tests.append(
            t_row(
                "Fig. 3",
                "b",
                f"WT vs edited at day {day}",
                "total cellular pigment content (pg cell−1)",
                sub.loc[sub.Strain == WT, "cellular_pigment_pg_cell"],
                sub.loc[sub.Strain == EDITED, "cellular_pigment_pg_cell"],
                "Two-tailed independent-samples Student t-test",
                True,
                f3.name,
                "Fig.3b",
                tiers=2,
            )
        )
    adjusted = multipletests([tests[i]["p_raw"] for i in pigment_indices], method="holm")[1]
    for i, p_adj in zip(pigment_indices, adjusted):
        tests[i]["notes"] = f"Holm sensitivity value across eight days: {p_adj:.6g}; not used for the figure annotation."

    growth = pd.read_excel(f3, sheet_name="Fig.3c", header=1)
    growth["density_cells_L"] = growth["Cell density (cells/L)"].map(parse_growth_density)
    growth["ln_density_cells_ml"] = np.log(growth["density_cells_L"] / 1000)
    growth = growth[growth.Days.between(5, 23)].copy()
    daily = growth.groupby(["Strain", "Days"], as_index=False)["ln_density_cells_ml"].mean()
    for strain, sub in daily.groupby("Strain"):
        fit = stats.linregress(sub.Days, sub.ln_density_cells_ml)
        tests.append(
            {
                "figure": "Fig. 3",
                "panel": "c",
                "comparison": f"Linear fit for {strain}",
                "outcome": "ln(cell density per mL), days 5–23",
                "group_1_mean": float(fit.slope),
                "group_1_sd": float(fit.stderr),
                "group_1_n": int(len(sub)),
                "group_2_mean": float(fit.intercept),
                "group_2_sd": float(fit.intercept_stderr),
                "group_2_n": int(len(sub)),
                "test": "Ordinary least-squares linear regression",
                "correction": "None",
                "statistic": float(fit.rvalue),
                "df": int(len(sub) - 2),
                "p_raw": float(fit.pvalue),
                "p_adjusted": None,
                "p_used_for_display": float(fit.pvalue),
                "p_formatted": format_p(float(fit.pvalue)),
                "significance": significance_label(float(fit.pvalue), tiers=2),
                "source_file": f3.name,
                "source_sheet": "Fig.3c",
                "notes": f"Slope={fit.slope:.6f}; intercept={fit.intercept:.6f}; R²={fit.rvalue**2:.6f}.",
            }
        )
    growth_indices = []
    for day, sub in growth.groupby("Days"):
        growth_indices.append(len(tests))
        tests.append(
            t_row(
                "Fig. 3",
                "c",
                f"WT vs edited at day {day}",
                "ln(cell density per mL)",
                sub.loc[sub.Strain == WT, "ln_density_cells_ml"],
                sub.loc[sub.Strain == EDITED, "ln_density_cells_ml"],
                "Two-tailed independent-samples Student t-test",
                True,
                f3.name,
                "Fig.3c",
                tiers=2,
                notes="Time-point comparison is tabulated but not annotated in the current panel.",
            )
        )
    adjusted = multipletests([tests[i]["p_raw"] for i in growth_indices], method="holm")[1]
    for i, p_adj in zip(growth_indices, adjusted):
        tests[i]["notes"] += f" Holm sensitivity value across ten plotted days: {p_adj:.6g}."

    q = pd.read_excel(f3, sheet_name="Fig.3d", header=2)
    qbio = q.groupby(["Samples", "Days", "Cell density (cells/ml)"], as_index=False).Ct.mean()
    qbio["Strain"] = np.where(qbio.Samples.str.contains("91"), WT, EDITED)
    qbio["transcripts_per_cell"] = 10 ** ((39.67 - qbio.Ct) / 3.59) / qbio["Cell density (cells/ml)"]
    qbio["log10_transcripts_per_cell"] = np.log10(qbio["transcripts_per_cell"])
    for day, sub in qbio.groupby("Days"):
        tests.append(
            t_row(
                "Fig. 3",
                "d",
                f"WT vs edited at day {day}",
                "log10(sxtA4 transcripts per cell)",
                sub.loc[sub.Strain == WT, "log10_transcripts_per_cell"],
                sub.loc[sub.Strain == EDITED, "log10_transcripts_per_cell"],
                "Two-tailed independent-samples Student t-test",
                True,
                f3.name,
                "Fig.3d",
                tiers=2,
                notes="Technical Ct replicates were averaged within each biological sample before conversion; the test uses log10 abundance.",
            )
        )
    checks.extend(
        [
            {
                "severity": "High",
                "location": "Source_data_fig3.xlsx / Fig.3d",
                "issue": "The Description column assigns both 54 and 91 samples to the same strain at a given day.",
                "resolution_in_scripts": "Strain is recovered from the stable sample prefix: 91 = wild type; 54 = edited.",
                "submission_action": "Correct the source-data Description column before submission.",
            },
            {
                "severity": "Medium",
                "location": "Methods / Fig. 3d",
                "issue": "The displayed significance is reproduced by testing log10-transformed abundance, but the current Methods does not state this transformation.",
                "resolution_in_scripts": "Biological-sample abundances are log10 transformed for the t-test; bars remain on the original scale.",
                "submission_action": "State the transformation explicitly in Methods and the figure legend.",
            },
        ]
    )

    toxin = pd.read_excel(f3, sheet_name="Fig. 3e", header=2)
    toxin = toxin[toxin.Samples.notna()].copy()
    toxin["day"] = toxin.Samples.str.extract(r"-(?:91|54)-(\d+)-")[0].astype(int)
    toxin["cellular_toxin_fmol_cell"] = (
        toxin.Total
        * toxin["LC-MS/MS extract volume (mL)"]
        * toxin["Dilution factor"]
        * 1e6
        / (toxin["Cell density (cells/ml)"] * toxin["Culture volume collected (mL)"])
    )
    for day, sub in toxin.groupby("day"):
        tests.append(
            t_row(
                "Fig. 3",
                "e",
                f"WT vs edited at day {day}",
                "total cellular toxin (fmol cell−1)",
                sub.loc[sub.Strain == WT, "cellular_toxin_fmol_cell"],
                sub.loc[sub.Strain == EDITED, "cellular_toxin_fmol_cell"],
                "Two-tailed Welch t-test",
                False,
                f3.name,
                "Fig. 3e",
                tiers=2,
                notes="Welch test is used because group variances differ substantially.",
            )
        )

    f4 = source_dir / "Source_data_fig4.xlsx"
    sxta = pd.read_excel(f4, sheet_name="Fig.4b")
    fig4b_indices = []
    for day, sub in sxta.groupby("day"):
        fig4b_indices.append(len(tests))
        tests.append(
            t_row(
                "Fig. 4",
                "b",
                f"WT vs edited at day {day}",
                "sxtA expression (FPKM)",
                sub.loc[sub.group == WT, "fpkm"],
                sub.loc[sub.group == EDITED, "fpkm"],
                "Two-tailed Welch t-test",
                False,
                f4.name,
                "Fig.4b",
            )
        )
    add_holm(tests, fig4b_indices)

    pst = pd.read_excel(f4, sheet_name="Fig.4c")
    fig4c_indices = []
    for gene, sub in pst.groupby("pst_gene", sort=False):
        fig4c_indices.append(len(tests))
        tests.append(
            t_row(
                "Fig. 4",
                "c",
                f"WT vs edited for {gene}",
                "FPKM pooled across four sampled days",
                sub.loc[sub.group == WT, "fpkm"],
                sub.loc[sub.group == EDITED, "fpkm"],
                "Two-tailed Welch t-test",
                False,
                f4.name,
                "Fig.4c",
            )
        )
    add_holm(tests, fig4c_indices)

    persistent = pd.read_excel(f4, sheet_name="Fig.4e")
    gene_p = persistent[
        [
            "gene_id",
            "dark_gene",
            "annotation_source_count",
            "response_class",
            "WT_5d_vs_M_5d_log2fc",
            "WT_5d_vs_M_5d_pval",
            "WT_9d_vs_M_9d_log2fc",
            "WT_9d_vs_M_9d_pval",
            "WT_13d_vs_M_13d_log2fc",
            "WT_13d_vs_M_13d_pval",
            "WT_15d_vs_M_15d_log2fc",
            "WT_15d_vs_M_15d_pval",
        ]
    ].replace([np.inf, -np.inf], ["Inf", "-Inf"])
    checks.append(
        {
            "severity": "High",
            "location": "Source_data_fig4.xlsx / Fig.4e and Supplementary Table 3",
            "issue": "Nominal DESeq2 P values are present, but adjusted P values stated as retained in the Methods are not included in the supplied source tables.",
            "resolution_in_scripts": "The workbook reports the available nominal P values without inventing adjusted values.",
            "submission_action": "Add DESeq2 adjusted P values to the source data or revise the Methods statement.",
        }
    )

    fed = source_dir / "Source_data_Extended_data_fig.xlsx"
    ed1b = pd.read_excel(fed, sheet_name="ED_Fig.1b", header=5)
    required_ed1b = {"log10(sxtA4-copies)", "Ct Value"}
    if not required_ed1b.issubset(ed1b.columns):
        raise ValueError("ED_Fig.1b source-data columns are incomplete.")
    ed1b_means = (
        ed1b.dropna(subset=list(required_ed1b))
        .groupby("log10(sxtA4-copies)", as_index=False)["Ct Value"]
        .mean()
        .sort_values("log10(sxtA4-copies)")
    )
    ed1b_regression = stats.linregress(ed1b_means["log10(sxtA4-copies)"], ed1b_means["Ct Value"])
    ed1b_n = len(ed1b_means)
    tests.append(
        {
            "figure": "Extended Data Fig. 1",
            "panel": "b",
            "comparison": "qPCR standard curve",
            "outcome": "Ct vs log10(copy number)",
            "group_1_mean": float(ed1b_regression.slope),
            "group_1_sd": None,
            "group_1_n": ed1b_n,
            "group_2_mean": float(ed1b_regression.intercept),
            "group_2_sd": None,
            "group_2_n": None,
            "test": "Linear regression",
            "correction": "None",
            "statistic": float(ed1b_regression.rvalue**2),
            "df": float(ed1b_n - 2),
            "p_raw": float(ed1b_regression.pvalue),
            "p_adjusted": None,
            "p_used_for_display": float(ed1b_regression.pvalue),
            "p_formatted": format_p(ed1b_regression.pvalue),
            "significance": significance_label(ed1b_regression.pvalue),
            "source_file": fed.name,
            "source_sheet": "ED_Fig.1b",
            "notes": "Regression used the mean Ct value at each log10 copy-number level; group_1_mean is the slope and group_2_mean is the intercept.",
        }
    )

    ed1c = pd.read_excel(fed, sheet_name="ED_Fig.1c", header=2)
    ed1c_bio = ed1c.groupby(["Number of cells", "Samples"], as_index=False)["Ct Value"].mean()
    ed1c_bio["copy_number_per_cell"] = (
        10 ** ((39.67 - ed1c_bio["Ct Value"]) / 3.59) * 80 / ed1c_bio["Number of cells"]
    )
    for cells, sub in ed1c_bio.groupby("Number of cells"):
        mean, sd, n = mean_sd(sub.copy_number_per_cell)
        summaries.append(
            {
                "figure": "Extended Data Fig. 1",
                "panel": "c",
                "group": f"{cells:g} cells",
                "outcome": "sxtA4 copy number per cell",
                "mean": mean,
                "sd": sd,
                "n": n,
                "source_file": fed.name,
                "source_sheet": "ED_Fig.1c",
                "notes": "Technical Ct replicates averaged per biological sample; ×80 corrects fourfold dilution and 1/20 template fraction.",
            }
        )

    ed7 = pd.read_excel(fed, sheet_name="ED_Fig.7")
    ed7_rows = ed7[
        ["module", "GO.ID", "term", "timepoint", "up", "down", "DEG", "FDR", "FDR_bin"]
    ].to_dict(orient="records")

    panel_map = [
        ("Fig. 1", "d", "Fig.1d", "plot_figure1.py", "Student t-tests"),
        ("Fig. 1", "f", "Fig.1f", "plot_figure1.py", "Student t-tests"),
        ("Fig. 1", "g,h", "Fig.1g_base_comp; Fig.1g_h_SNPs", "plot_figure1.py", "Descriptive"),
        ("Fig. 2", "a–c", "Injection_data; Fig.2a–c", "plot_figure2.py", "ANOVA/Tukey and t-tests"),
        ("Fig. 2", "d–f", "Fig.2d-f", "plot_figure2.py", "Welch t-test; ANOVA/Tukey"),
        ("Fig. 2", "g–i", "Fig.2g–i", "plot_figure2.py", "Descriptive"),
        ("Fig. 3", "b", "Fig.3b", "plot_figure3.py", "Time-point t-tests"),
        ("Fig. 3", "c", "Fig.3c", "plot_figure3.py", "Regression and time-point t-tests"),
        ("Fig. 3", "d", "Fig.3d", "plot_figure3.py", "Log10-scale t-tests"),
        ("Fig. 3", "e", "Fig. 3e", "plot_figure3.py", "Welch t-tests"),
        ("Fig. 4", "a", "Fig.4a", "plot_figure4.py", "Descriptive"),
        ("Fig. 4", "b", "Fig.4b", "plot_figure4.py", "Welch tests with Holm correction"),
        ("Fig. 4", "c", "Fig.4c", "plot_figure4.py", "Welch tests with Holm correction"),
        ("Fig. 4", "d,e", "Fig.4d; Fig.4e", "plot_figure4.py", "DESeq2 gene-level tests"),
        ("Extended Data Fig. 1", "b,c", "ED_Fig.1b; ED_Fig.1c", "plot_extended_data.py", "Regression/descriptive"),
        ("Extended Data Fig. 2", "a,b", "ED_Fig.2a; ED_Fig.2b", "plot_extended_data.py", "Descriptive heat maps"),
        ("Extended Data Fig. 6", "a–e", "ED_Fig.6a–e", "plot_extended_data.py", "PCA/Pearson/descriptive"),
        ("Extended Data Fig. 7", "", "ED_Fig.7", "plot_extended_data.py", "GO enrichment FDR"),
    ]
    panel_rows = [
        {
            "figure": fig,
            "panel": panel,
            "source_file": (
                "Source_data_Extended_data_fig.xlsx" if fig.startswith("Extended") else f"Source_data_fig{fig.split('.')[1].strip()}.xlsx"
            ),
            "source_sheet": sheet,
            "script": script,
            "statistics": stat,
        }
        for fig, panel, sheet, script, stat in panel_map
    ]

    return {
        "panel_map": panel_rows,
        "tests": tests,
        "summaries": summaries,
        "fig4e_gene_p": gene_p.to_dict(orient="records"),
        "ed7_go_fdr": ed7_rows,
        "checks": checks,
    }


def write_excel_results(result: dict, path: Path) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name in ("panel_map", "tests"):
            rows = result[sheet_name]
            frame = pd.DataFrame(rows)
            frame.to_excel(writer, sheet_name=sheet_name, index=False)
            sheet = writer.sheets[sheet_name]
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            sheet.sheet_view.showGridLines = False
            for cell in sheet[1]:
                cell.font = Font(bold=True, color="FFFFFF")
                cell.fill = PatternFill("solid", fgColor="315F78")
                cell.alignment = Alignment(horizontal="center", vertical="center")
            for column_index, column_name in enumerate(frame.columns, 1):
                values = [str(column_name)] + [str(value) for value in frame[column_name].dropna()]
                width = min(max(max(len(value) for value in values) + 2, 10), 45)
                sheet.column_dimensions[get_column_letter(column_index)].width = width
            for row in sheet.iter_rows(min_row=2):
                for cell in row:
                    cell.alignment = Alignment(vertical="top", wrap_text=isinstance(cell.value, str) and len(cell.value) > 40)


def write_results(source_dir: Path, output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    result = calculate_statistics(source_dir)
    json_path = output_dir / "statistics_results.json"
    excel_path = output_dir / "statistics_results.xlsx"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    write_excel_results(result, excel_path)
    for name in ("checks.csv", "ed7_go_fdr.csv", "fig4e_gene_p.csv", "panel_map.csv", "summaries.csv", "tests.csv"):
        path = output_dir / name
        if path.exists():
            path.unlink()
    return json_path, excel_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Recalculate figure statistics from the source-data workbooks.")
    parser.add_argument("--source", type=Path, required=True, help="Folder containing the source-data workbooks.")
    parser.add_argument("--output", type=Path, required=True, help="Folder for tabulated results.")
    args = parser.parse_args()
    for path in write_results(args.source, args.output):
        print(path)


if __name__ == "__main__":
    main()
