# UvrD-assisted CRISPR editing of a dinoflagellate giant genome illuminates saxitoxin biosynthesis

This repository contains the Python scripts used to generate the data-driven
figure panels and statistical outputs for the article *UvrD-assisted CRISPR
editing of a dinoflagellate giant genome illuminates saxitoxin biosynthesis*.

## Repository contents

- `run_all.py`: main entry point for all analyses and figures.
- `scripts/plot_figure1.py` to `scripts/plot_figure5.py`: current main-figure plotting scripts.
- `scripts/plot_extended_data.py`: Extended Data figure plotting script.
- `scripts/statistics.py`: statistical analyses and tabulated results.
- `scripts/validate_outputs.py`: regression checks against the source-data
  Statistics sheets and structural checks for every generated artifact.
- `requirements.txt`: required Python packages.

The nine source-data workbooks are distributed separately with the article and
are read at run time.

## Requirements and installation

Python 3.10 or later is required. Create and activate a virtual environment,
then install the dependencies:

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Required source-data files

Download the source-data files associated with the article. The scripts expect
the following workbook names:

- `Source_data_fig1.xlsx`
- `Source_data_fig2.xlsx`
- `Source_data_fig3.xlsx`
- `Source_data_fig4.xlsx`
- `Source_data_fig5.xlsx`
- `Source_data_Extended_Data_Fig1.xlsx`
- `Source_data_Extended_Data_Fig2.xlsx`
- `Source_data_Extended_Data_Fig6.xlsx`
- `Source_data_Extended_Data_Fig7.xlsx`

Place all nine workbooks either in this repository directory, in its parent
directory, or in a separate directory supplied with `--source`.

## Usage

If the workbooks are in this directory or its parent directory, run:

```bash
python run_all.py
```

To use another source-data directory, run:

```bash
python run_all.py --source "/path/to/source-data"
```

The command creates:

- `figures/`, containing SVG, PDF and PNG figure outputs;
- `results/statistics_results.json`, containing machine-readable results; and
- `results/statistics_results.xlsx`, containing one concise `statistical_tests`
  worksheet with the tests, sample sizes, multiplicity adjustments, exact P
  values and figure labels needed for review.

Existing files with the same output names are replaced.

After generation, `run_all.py` verifies that all 55 statistical rows reproduce
the values recorded in the source-data Statistics sheets, validates all 33
SVG/PDF/PNG artifacts and confirms that the nine input workbooks retain their
original SHA-256 hashes. Any mismatch stops the run with a non-zero exit code.

## Reproducibility notes

The scripts generate figures and statistical outputs from the nine source-data
workbooks without modifying the input files. SVG and PDF are the preferred
editable formats, and PNG files are convenient previews; TIFF is not generated.

The analyses reported in the manuscript used Python 3.11.3, SciPy 1.10.0 and
statsmodels 0.13.5. The dependency ranges in `requirements.txt` support current
compatible environments; the built-in regression checks protect the reported
statistics from numerical or workbook-schema drift.

## Citation

If you use or adapt these scripts, please cite the associated article:

Chen, Z. *et al.* *UvrD-assisted CRISPR editing of a dinoflagellate giant
genome illuminates saxitoxin biosynthesis*. Citation details and DOI will be
added after publication.

## Questions and issues

Please use the GitHub issue tracker to report reproducibility problems or ask
questions about the scripts.

## License

The code in this repository is released under the MIT License. See `LICENSE`
for details. The source-data workbooks distributed with the article are not
covered by this repository license.
