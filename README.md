# DMWA propagation medium: code and source data

Code and derived source data accompanying **Directed multi-order networks reveal resource–spectral rebalancing of the seizure propagation medium**.

This version follows the manuscript and Supplementary Materials dated **17 September 2026**: **Figs. 1–7, Figs. S1–S8 and Tables S1–S3**. Document hashes are recorded in [the snapshot manifest](docs/manuscript_snapshot.json). The repository does not distribute the manuscript files.

## Scientific scope

Directed multi-order weighted adjacency (DMWA), denoted by `W`, represents structural coupling among overlapping regional groups. The propagation medium parameter is `eta = R / rho(W)`, where `R` is mean outgoing DMWA weight and `rho(W)` is the spectral radius. Across 24 seizures from 14 patients, mean coupling increased proportionally more than spectral radius, producing an ictal and postictal increase in eta.

- **Fig. 6:** retaining group order and exact overlap degree recovers the projected response under the specified linear model. Typical node-only approximation errors were small; exact recovery is not necessarily dimensional compression.
- **Fig. 7:** with uniform leakage and mean-weight normalization, eta sets the critical normalized feedback gain and asymptotic recovery scale. These are model relationships, not independent physiological predictions.

The current term is **propagation medium parameter**, replacing the earlier effective-refractive-index analogy. Structural directionality does not imply causal neural transmission. Groups are pairs and closed triangles extracted from a thresholded pairwise PLV scaffold.

## Contents

| Path | Contents |
|---|---|
| [Source_Data.xlsx](Source_Data.xlsx) | Panel data for Figs. 2–7 and S1–S8, updated Tables S1–S3, Data_Map and Data_Dictionary. Fig. 1 is a construction schematic. |
| [docs/figure_map.md](docs/figure_map.md) | Current questions, source sheets, code entry points and numbering. |
| [mechanism/](mechanism/) | Portable Python analyses, design records, CSVs and 120 de-identified derived DMWA matrices for Figs. 6–7. |
| [r_figures/](r_figures/) | R scripts arranged by current figure number. |
| [matlab/](matlab/) | Original empirical, control, ablation, robustness and simulation scripts. Historical directory names are mapped in the figure index. |

## Reproduce the added mechanism analyses

Use Python 3.11 or later:

```bash
python -m pip install -r mechanism/requirements.txt
python mechanism/code/w_mechanism_analysis.py
python mechanism/code/eta_response_analysis.py
```

The W script regenerates the 100 formal scaffolds and three illustrative examples for Fig. 6. The eta script regenerates the synthetic examples and all 360 responses from the 120 supplied matrices. No raw SEEG or private MATLAB cache is required. For a short synthetic demonstration:

```bash
python mechanism/code/eta_response_analysis.py --examples-only
```

By default, scripts replace derived results within `mechanism/data/` and write diagnostics to `mechanism/reports/`. Set `DMWA_OUTPUT_ROOT` to a separate directory to preserve distributed results. `DMWA_INPUT_DIR` can select another folder containing the Fig. 7 selection table and `matrices/`. See [mechanism/README.md](mechanism/README.md) for equations, settings and units.

## Draw the current Figs. 6 and 7

Plotting reads the included CSVs; rerunning simulations is optional.

```r
install.packages(c("ggplot2", "patchwork", "dplyr", "tidyr", "readr",
                   "ggtext", "svglite", "ragg", "jsonlite"))
```

```bash
Rscript r_figures/Fig_6/make_fig6_group_context.R
Rscript r_figures/Fig_7/make_fig7_feedback_response.R
```

Use patchwork 1.3.0 or later, including `free()`. Outputs are SVG, PDF, 600-dpi PNG/TIFF and a 300-dpi preview under `outputs/Fig_6/` and `outputs/Fig_7/`. Arial matches the manuscript; font substitution can alter spacing. Scripts preserve the latest Fig. 6D and Fig. 7 layouts. Input overrides are `DMWA_W_DATA_DIR` and `DMWA_ETA_DATA_DIR`; output overrides are `DMWA_FIG6_OUT` and `DMWA_FIG7_OUT`.

Older R scripts are retained as analysis/plotting provenance, with current folder numbering. Some require analysis-specific CSV or MATLAB intermediates that are not distributed. Their historical filenames and environment-variable prefixes remain unchanged. **Source_Data.xlsx supplies the public panel-level values; it is not a drop-in replacement for every legacy input file.** Install openpyxl and run `python tools/extract_source_data.py` to extract its source blocks. Fig. 1 has no numerical plotting script.

## Analysis units and interpretation

- Primary inference for the eta stage effect uses **14 patient-level averages** (Fig. S2), each formed by equally averaging that patient's seizure-stage estimates. The 24 seizure-level summaries in Figs. 2–5 are exploratory.
- Fig. 3 uses arithmetic stage means. Component bars describe log changes in the ratio of stage-mean R to stage-mean spectral radius; the eta points are separately averaged. These summaries need not obey the window-level log identity. Panels B–D use 100 times log changes.
- Fig. S6 AUC uses seizure-level stage summaries; Fig. S7 AUC uses windows. Both describe within-cohort discrimination, without patient-held-out validation.
- Fig. S5 contains 50 descriptive sensitivity settings; P values are unadjusted across settings.
- Fig. 5 burden is an eigen-sensitivity descriptor, not a deletion effect or causal propagation contribution.
- Fig. S8 does not establish a monotonic eta response to increasing simulated epileptiform drive.
- Eta is invariant to uniform scaling and transposition, with no universal upper bound of one. Its threshold interpretation requires the specified leakage and normalization.

## Data availability

The repository contains de-identified derived values and group-coupling matrices. Mechanism identifiers `S01`–`S24` and `P01`–`P14` correspond to Table S1 with zero padding. Patient labels were harmonized without changing seizure grouping or numerical values.

Raw/preprocessed SEEG, recording dates, electrode labels, names, private clinical identifiers and private linkage keys are excluded. Access to raw recordings, where permissible, requires review by the responsible clinical institution and appropriate ethics and data-use arrangements.

The journal Source Data file should use the same version. This workbook is a versioned public mirror, not a substitute for the journal upload.

See [CHANGELOG.md](CHANGELOG.md) and [validation.md](docs/validation.md) for this update and the checks performed.
