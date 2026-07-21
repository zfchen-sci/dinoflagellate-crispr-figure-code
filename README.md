# Figure and statistical analysis scripts

Python scripts used to generate the data-driven figure panels and statistical
outputs for the article *Editing a dinoflagellate giant genome illuminates
saxitoxin biosynthesis*.

## Repository contents

- `run_all.py`: main entry point for all analyses and figures.
- `scripts/plot_figure1.py` to `scripts/plot_figure4.py`: main-figure plotting scripts.
- `scripts/plot_extended_data.py`: Extended Data figure plotting script.
- `scripts/statistics.py`: statistical analyses and tabulated results.
- `requirements.txt`: required Python packages.

The source-data Excel workbooks are distributed with the article and are not
duplicated in this repository.

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
- `Source_data_Extended_data_fig.xlsx`

Place all five workbooks either in this repository directory, in its parent
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
- `results/statistics_results.xlsx`, containing tabulated statistical results.

Existing files with the same output names are replaced.

## Reproducibility notes

The scripts validate required workbook names, worksheet structures and several
expected data dimensions before generating outputs. They do not modify the
source-data workbooks. SVG and PDF are the preferred editable and publication
formats; PNG files are provided for convenient previewing.

## Citation

If you use or adapt these scripts, please cite the associated article:

Chen, Z. *et al.* *Editing a dinoflagellate giant genome illuminates saxitoxin
biosynthesis*. Citation details and DOI will be added after publication.

## Questions and issues

Please use the GitHub issue tracker to report reproducibility problems or ask
questions about the scripts.

## License

The code in this repository is released under the MIT License. See `LICENSE`
for details. The source-data workbooks distributed with the article are not
covered by this repository license.
