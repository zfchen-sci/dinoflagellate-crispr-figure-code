from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


REQUIRED_SOURCE_FILES = (
    "Source_data_fig1.xlsx",
    "Source_data_fig2.xlsx",
    "Source_data_fig3.xlsx",
    "Source_data_fig4.xlsx",
    "Source_data_Extended_data_fig.xlsx",
)


def run(script: Path, source: Path, output: Path) -> None:
    command = [sys.executable, str(script), "--source", str(source), "--output", str(output)]
    subprocess.run(command, check=True)


def resolve_source(root: Path, requested: Path | None) -> Path:
    candidates = [requested.resolve()] if requested else [root, root.parent]
    reports = []
    for folder in candidates:
        missing = [name for name in REQUIRED_SOURCE_FILES if not (folder / name).is_file()]
        if not missing:
            return folder
        reports.append(f"{folder}: missing {', '.join(missing)}")
    detail = "\n".join(reports)
    raise FileNotFoundError(
        f"Required source-data workbooks were not found.\n{detail}\n"
        'Use: python run_all.py --source "path\\to\\source-data"'
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Recalculate statistics and recreate Excel-based figure panels.")
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help="Folder containing the Source_data_*.xlsx files.",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    source = resolve_source(root, args.source)
    scripts = root / "scripts"
    results = root / "results"
    figures = root / "figures"
    results.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)

    run(scripts / "statistics.py", source, results)
    for name in ("plot_figure1.py", "plot_figure2.py", "plot_figure3.py", "plot_figure4.py", "plot_extended_data.py"):
        run(scripts / name, source, figures)

    print(f"Source data: {source}")
    print(f"Statistics: {results}")
    print(f"Figures: {figures}")


if __name__ == "__main__":
    main()
