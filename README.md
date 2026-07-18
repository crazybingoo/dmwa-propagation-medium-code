# DMWA Propagation Medium Code and Source Data

This repository contains the current analysis and figure code for Fig. 2–Fig. 9 and Supplementary Fig. 1–Supplementary Fig. 5 of the manuscript on directed higher-order propagation-medium reorganization in epileptic seizures.

## Contents

- `Source_Data.xlsx`: de-identified numerical source data underlying every figure panel. Each panel has an independent worksheet; Fig. 9a uses separate trace and PLV worksheets.
- `matlab/`: MATLAB scripts for the DMWA, control, robustness, representation-comparison and simulation analyses.
- `r_figures/`: R scripts for the final manuscript and Supplementary Information figures.

The workbook contains a `README`, a panel-to-sheet `Data_Map`, and a `Data_Dictionary`. Analysis units are stated explicitly because they differ across panels.

## Final R figure scripts

The final figure versions use 10.5-pt main text and 12.5-pt panel labels.

- Fig. 2: `r_figures/Fig_2/make_fig2_effective_refractive_index_redraw_font10p5_label12p5_review.R`
- Fig. 3: run `r_figures/Fig_3/render_fig3_delta_log_notation.R`
- Fig. 4: run `r_figures/Fig_4/render_fig4_control_label_nouns.R`
- Fig. 5: `r_figures/Fig_5/make_fig5_eta_burden_font10p5_label12p5_legends_below.R`
- Fig. 6: run `r_figures/Fig_6/render_fig6_retained_fraction_d_title_bold.R`
- Fig. 7: run `r_figures/Fig_7/render_fig7c_degree_bias_spacing.R`
- Fig. 8: run `r_figures/Fig_8/render_fig8c_order_labels.R`
- Fig. 9: `r_figures/Fig_9/draw_fig9_hr_mechanism_final.R`

The supplementary-figure scripts are located in the corresponding `r_figures/Supplementary_Fig_*` directories.

Scripts resolve their own directory by default. Approved users can override input/output locations with the environment variables documented in each script, such as `FIG2_DATA_DIR`, `FIG5_DATA_DIR`, `FIG8_SOURCE_DIR`, or `FIG9_SOURCE_DIR`. No personal workstation paths are embedded in the public code.

## Data availability and privacy boundary

`Source_Data.xlsx` contains only de-identified, derived numerical values used for the figures and statistical summaries. Seizures are labelled `Seizure 01`–`Seizure 24`, and patients are labelled `P01`–`P14`. The private linkage keys are not included.

The repository intentionally excludes:

- raw or preprocessed SEEG recordings;
- participant names, clinical identifiers, recording dates and electrode labels;
- internal case codes and patient-linkage tables;
- local absolute file paths and workstation-specific information;
- intermediate files that could reconnect public labels to private clinical records.

The raw SEEG recordings are human-participant clinical neurophysiology data and are not publicly redistributable through GitHub because of participant privacy, ethics approval, consent terms and clinical data-use restrictions. Access, when permissible, requires review through the responsible clinical institution, appropriate ethics approval and a data-use agreement.

For journal submission, `Source_Data.xlsx` should also be uploaded as the manuscript's formal Source Data file; the GitHub copy is a versioned public mirror and does not replace the journal upload.

## Analysis notes

For the main empirical DMWA construction, one PLV threshold is selected per seizure from the elbow of its density–threshold curve and applied to all windows from the same seizure. The Fig. 6 robustness analysis uses 50 settings spanning window lengths from 1 to 5 s and retained PLV edge fractions from 0.45 to 0.90.
