# Figure and statistical analysis scripts

Python scripts used to generate the data-driven figure panels and statistical
outputs for the article *Editing a dinoflagellate giant genome illuminates
saxitoxin biosynthesis*.

## Repository contents

- `run_all.py`: main entry point for all analyses and figures.
- `scripts/plot_figure1.py` to `scripts/plot_figure5.py`: current main-figure plotting scripts.
- `scripts/plot_extended_data.py`: Extended Data figure plotting script.
- `scripts/statistics.py`: statistical analyses and tabulated results.
- `Statistics_analysis_by_figure_V3.2.xlsx`: by-figure results for figures with statistical tests.
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
- `Source_data_fig5.xlsx`
- `Source_data_Extended_data_fig.xlsx`

Place all six workbooks either in this repository directory, in its parent
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

## Reproducibility notes

The scripts validate required workbook names, worksheet structures and expected
data dimensions before generating outputs. They do not modify the source-data
workbooks. SVG and PDF are the preferred editable and publication formats; PNG
files are provided only as convenient previews. No TIFF files are generated.
All figure text uses Times New Roman, including math-text fallbacks.

The current mapping is Figure 1 (copy-aware target and free-crRNA predictions),
Figure 2 (microinjection and UvrD optimization), Figure 3 (delivery, editing and
repair), Figure 4 (molecular and toxin phenotypes) and Figure 5
(transcriptomics). Cell recovery in Figure 3b–e is descriptive and displays the
three independent cell sets. Figure 3f uses one-way ANOVA with Tukey HSD;
Figure 3g and 3h use two-sided Welch t-tests; Figure 3i uses Welch one-way ANOVA
with Games–Howell pairwise tests. Figure 4d and 4e and Figure 5b and 5c use
two-sided Welch tests with the multiplicity adjustments recorded in
`statistics_results.xlsx`. Figure 5e is descriptive and is not subjected to an
additional hypothesis test.

Extended Data Figure 1d–f is recreated from the free-RNA prediction tables.
Extended Data Figure 1g is a qualitative gel assay: its workbook sheet records
the displayed lane map and splice disclosure, and the scripts deliberately do
not invent densitometry.

## Panel-specific drawing boundaries

`plot_figure1.py` intentionally redraws only Fig. 1b and Fig. 1f. Fig. 1a is
excluded, and Fig. 1e is excluded because the formal panel contains RNA
secondary-structure diagrams; a horizontal MFE bar chart would not be a
faithful replacement. Fig. 3b-e uses the same raw-point plus mean +/- s.d.
geometry as the formal main figure. ED Fig. 1f uses the same teal probability
scale, mature-crRNA boundaries and segment labels as the formal Extended Data
figure.

The source-data `Description` field is treated as original metadata and is not
rewritten by these scripts. Analysis grouping derived from sample identifiers
is explicit in the code and should be checked against the authors' final sample
map before submission.

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
