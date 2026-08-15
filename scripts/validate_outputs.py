from __future__ import annotations

import argparse
import json
import math
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


REQUIRED_SOURCE_FILES = (
    "Source_Data_Figure_1.xlsx",
    "Source_Data_Figure_2.xlsx",
    "Source_Data_Figure_3.xlsx",
    "Source_Data_Figure_4.xlsx",
    "Source_Data_Supplementary_Figure_S1.xlsx",
    "Source_Data_Supplementary_Figure_S2.xlsx",
)

EXPECTED_FIGURE_STEMS = (
    "Figure1_data_panels",
    "Figure2_data_panels",
    "Figure3_data_panels",
    "Figure4_data_panels",
    "Supplementary_Figure_S1b",
    "Supplementary_Figure_S1c",
    "Supplementary_Figure_S1d_f",
    "Supplementary_Figure_S2",
)

RECALCULATED_PANELS = {
    2: {"d", "f"},
    3: {"f", "g", "h", "i"},
    4: {"d", "e"},
}


def _as_float(value: object) -> float | None:
    if value is None or pd.isna(value):
        return None
    return float(value)


def _assert_close(actual: object, expected: object, label: str) -> None:
    actual_number = _as_float(actual)
    expected_number = _as_float(expected)
    if actual_number is None and expected_number is None:
        return
    if actual_number is None or expected_number is None:
        raise AssertionError(f"{label}: {actual!r} != {expected!r}")
    if not math.isclose(actual_number, expected_number, rel_tol=1e-9, abs_tol=1e-12):
        raise AssertionError(f"{label}: {actual_number:.12g} != {expected_number:.12g}")


def _load_expected_statistics(source_dir: Path, figure_number: int) -> pd.DataFrame:
    source = source_dir / f"Source_Data_Figure_{figure_number}.xlsx"
    expected = pd.read_excel(source, sheet_name="Statistics", header=1)
    expected = expected[expected["Panel"].notna()].reset_index(drop=True)
    return expected


def validate_statistics(source_dir: Path, results_dir: Path) -> int:
    json_path = results_dir / "statistics_results.json"
    excel_path = results_dir / "statistics_results.xlsx"
    result = json.loads(json_path.read_text(encoding="utf-8"))
    tests = result.get("tests", [])

    expected_total = 0
    for figure_number in range(2, 5):
        figure = f"Fig. {figure_number}"
        expected = _load_expected_statistics(source_dir, figure_number)
        actual = [row for row in tests if row["figure"] == figure]
        recalculated_panels = {str(row["panel"]) for row in actual}
        if recalculated_panels != RECALCULATED_PANELS[figure_number]:
            raise AssertionError(
                f"{figure}: recalculated panels {sorted(recalculated_panels)}; "
                f"expected {sorted(RECALCULATED_PANELS[figure_number])}"
            )
        expected = expected[
            expected["Panel"].astype(str).isin(recalculated_panels)
        ].reset_index(drop=True)
        expected_total += len(expected)
        if len(actual) != len(expected):
            raise AssertionError(
                f"{figure}: generated {len(actual)} statistical rows; expected {len(expected)}"
            )
        for index, (generated, recorded) in enumerate(zip(actual, expected.to_dict("records")), start=1):
            prefix = f"{figure} statistics row {index}"
            if str(generated["panel"]) != str(recorded["Panel"]):
                raise AssertionError(
                    f"{prefix}: panel {generated['panel']!r} != {recorded['Panel']!r}"
                )
            _assert_close(generated["group_1_n"], recorded["n1"], f"{prefix} n1")
            _assert_close(generated["group_2_n"], recorded["n2"], f"{prefix} n2")
            _assert_close(generated["statistic"], recorded["Statistic"], f"{prefix} statistic")
            _assert_close(generated["p_raw"], recorded["P"], f"{prefix} raw P")
            _assert_close(
                generated["p_used_for_display"],
                recorded["Adjusted P"],
                f"{prefix} reported P",
            )

    if len(tests) != expected_total:
        raise AssertionError(
            f"Generated {len(tests)} total statistical rows; "
            f"source Statistics sheets contain {expected_total}"
        )
    if result.get("checks"):
        raise AssertionError(f"Unresolved analysis checks remain: {result['checks']}")

    workbook = pd.read_excel(excel_path, sheet_name="statistical_tests")
    if len(workbook) != len(tests):
        raise AssertionError(
            f"Statistics workbook contains {len(workbook)} rows; JSON contains {len(tests)}"
        )
    for index, (generated, exported) in enumerate(
        zip(tests, workbook.to_dict("records")), start=1
    ):
        prefix = f"statistics export row {index}"
        _assert_close(generated["p_raw"], exported["p_raw"], f"{prefix} raw P")
        _assert_close(
            generated["p_used_for_display"], exported["p_reported"], f"{prefix} reported P"
        )
        if generated["significance"] != exported["figure_label"]:
            raise AssertionError(
                f"{prefix}: figure label {generated['significance']!r} != {exported['figure_label']!r}"
            )
    return len(tests)


def _validate_png(path: Path) -> None:
    data = path.read_bytes()
    if not data.startswith(b"\x89PNG\r\n\x1a\n") or len(data) < 24:
        raise AssertionError(f"Invalid PNG: {path.name}")
    width, height = struct.unpack(">II", data[16:24])
    if width < 100 or height < 100:
        raise AssertionError(f"Implausible PNG dimensions for {path.name}: {width} x {height}")


def _validate_pdf(path: Path) -> None:
    data = path.read_bytes()
    if not data.startswith(b"%PDF-") or b"%%EOF" not in data[-2048:]:
        raise AssertionError(f"Invalid PDF: {path.name}")


def _validate_svg(path: Path) -> None:
    root = ET.parse(path).getroot()
    if root.tag.rsplit("}", 1)[-1] != "svg":
        raise AssertionError(f"Invalid SVG root element: {path.name}")


def validate_figures(figures_dir: Path) -> int:
    expected = {
        f"{stem}.{suffix}"
        for stem in EXPECTED_FIGURE_STEMS
        for suffix in ("svg", "pdf", "png")
    }
    actual = {
        path.name
        for path in figures_dir.iterdir()
        if path.is_file() and path.suffix in {".svg", ".pdf", ".png"}
    }
    missing = sorted(expected - actual)
    unexpected = sorted(actual - expected)
    if missing or unexpected:
        raise AssertionError(f"Figure artifact mismatch; missing={missing}, unexpected={unexpected}")

    for name in sorted(expected):
        path = figures_dir / name
        if path.stat().st_size < 512:
            raise AssertionError(f"Figure artifact is unexpectedly small: {name}")
        validators = {".png": _validate_png, ".pdf": _validate_pdf, ".svg": _validate_svg}
        validators[path.suffix](path)
    return len(expected)


def validate_run(source_dir: Path, results_dir: Path, figures_dir: Path) -> dict[str, int]:
    missing = [name for name in REQUIRED_SOURCE_FILES if not (source_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing source-data workbooks: {', '.join(missing)}")
    return {
        "statistical_rows": validate_statistics(source_dir, results_dir),
        "figure_artifacts": validate_figures(figures_dir),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate statistics and figure artifacts.")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()
    summary = validate_run(args.source, args.results, args.figures)
    print(
        "Validation passed: "
        f"{summary['statistical_rows']} statistical rows and "
        f"{summary['figure_artifacts']} figure artifacts."
    )


if __name__ == "__main__":
    main()
