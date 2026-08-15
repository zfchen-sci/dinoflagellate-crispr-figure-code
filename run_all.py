from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

from scripts.validate_outputs import (
    EXPECTED_FIGURE_STEMS,
    REQUIRED_SOURCE_FILES,
    validate_run,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clear_generated_outputs(results: Path, figures: Path) -> None:
    for name in ("statistics_results.json", "statistics_results.xlsx"):
        (results / name).unlink(missing_ok=True)
    for stem in EXPECTED_FIGURE_STEMS:
        for suffix in ("svg", "pdf", "png"):
            (figures / f"{stem}.{suffix}").unlink(missing_ok=True)


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
    parser = argparse.ArgumentParser(
        description="Recalculate statistics and recreate Excel-based figure panels."
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help="Folder containing the six Source_Data_*.xlsx files.",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    source = resolve_source(root, args.source)
    scripts = root / "scripts"
    results = root / "results"
    figures = root / "figures"
    results.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    source_hashes = {name: sha256(source / name) for name in REQUIRED_SOURCE_FILES}
    clear_generated_outputs(results, figures)

    run(scripts / "statistics.py", source, results)
    for name in (
        "plot_figure1.py",
        "plot_figure2.py",
        "plot_figure3.py",
        "plot_figure4.py",
        "plot_extended_data.py",
    ):
        run(scripts / name, source, figures)

    final_hashes = {name: sha256(source / name) for name in REQUIRED_SOURCE_FILES}
    if final_hashes != source_hashes:
        changed = sorted(
            name for name in source_hashes if source_hashes[name] != final_hashes[name]
        )
        raise RuntimeError(f"Source-data workbooks changed during analysis: {', '.join(changed)}")
    validation = validate_run(source, results, figures)

    print(f"Source data: {source}")
    print(f"Statistics: {results}")
    print(f"Figures: {figures}")
    print(
        "Validation passed: "
        f"{validation['statistical_rows']} statistical rows and "
        f"{validation['figure_artifacts']} figure artifacts; source SHA-256 hashes unchanged."
    )


if __name__ == "__main__":
    main()
