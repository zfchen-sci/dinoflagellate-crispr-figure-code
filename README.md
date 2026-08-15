# UvrD-assisted CRISPR editing of a dinoflagellate giant genome illuminates saxitoxin biosynthesis

This repository contains the Python scripts used to recreate the data-driven figure panels and statistical outputs for the associated article.

## Contents

- `run_all.py`: runs all analyses, figures and validation checks.
- `scripts/plot_figure1.py` to `scripts/plot_figure4.py`: main-figure plotting scripts.
- `scripts/plot_extended_data.py`: Supplementary Figures S1 and S2 plotting script.
- `scripts/statistics.py`: statistical analyses and tabulated results.
- `scripts/validate_outputs.py`: regression and output-integrity checks.
- `requirements.txt`: Python dependencies.

## Installation

Python 3.10 or later is required.

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Source data

Download the six source-data workbooks distributed with the article:

- `Source_Data_Figure_1.xlsx`
- `Source_Data_Figure_2.xlsx`
- `Source_Data_Figure_3.xlsx`
- `Source_Data_Figure_4.xlsx`
- `Source_Data_Supplementary_Figure_S1.xlsx`
- `Source_Data_Supplementary_Figure_S2.xlsx`

Place them in this repository, its parent directory, or another directory supplied with `--source`.

## Usage

```bash
python run_all.py
```

For a separate source-data directory:

```bash
python run_all.py --source "/path/to/source-data"
```

The command creates editable SVG/PDF figures, PNG previews, and JSON/XLSX statistical results. It also verifies the 46 statistical rows against the workbooks, checks all 24 figure artifacts, and confirms that the six input files retain their original SHA-256 hashes. Any mismatch stops the run.

TIFF files are not generated.

## Citation

Chen, Z. *et al.* *UvrD-assisted CRISPR editing of a dinoflagellate giant genome illuminates saxitoxin biosynthesis*. Citation details and DOI will be added after publication.

## License

The code is released under the MIT License. Source-data workbooks distributed with the article are not covered by this repository license.
