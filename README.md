# DMWA Propagation Medium Code

This repository contains code used to generate and analyse Fig. 2-Fig. 9 for the manuscript on directed higher-order propagation-medium reorganization in epileptic seizures.

The repository is intentionally code-only. Raw SEEG recordings, patient-level data, intermediate result tables, MATLAB data files, rendered figures, and other derived outputs are not included.

## Contents

- `matlab/`: MATLAB scripts copied from the Fig. 2-Fig. 9 code folders under `2_Fig*`.
- `r_figures/`: the latest R script for each manuscript figure from `nature_fig/Fig_2` to `nature_fig/Fig_9`.

## Data and Results

No clinical SEEG data or result datasets are distributed in this repository. The scripts assume that required data files are available locally under the user's approved data-access environment.

The following file types are excluded from version control by design: MATLAB data files, spreadsheets, CSV/TSV result tables, rendered figures, PDFs, TIFFs, images, electrophysiology recordings, and other large or sensitive data formats.

## Notes

Some scripts may contain local file paths from the analysis workstation. These paths should be adapted to the local data-access environment before use.

For the main empirical DMWA construction, the MATLAB scripts select one PLV threshold per seizure from the elbow of the density-threshold curve and apply that threshold to all windows from the same seizure. The Fig. 6 robustness scripts intentionally override this with a 50-setting grid spanning window lengths from 1 to 5 s and retained PLV edge fractions from 0.45 to 0.90.

All paths in the public scripts are example paths. Raw data, intermediate result tables, rendered figures and manuscript files are intentionally excluded from the repository.
