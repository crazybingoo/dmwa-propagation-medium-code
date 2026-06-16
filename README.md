# DMWA Propagation Medium Code

This repository contains code used to generate and analyse Fig. 2-Fig. 9 and Supplementary Fig. 1-Supplementary Fig. 4 for the manuscript on directed higher-order propagation-medium reorganization in epileptic seizures.

This GitHub release is a code release. Clinical data access is governed by the Data Availability route described below.

## Contents

- `matlab/`: MATLAB scripts for the final analyses reported in Fig. 2-Fig. 9.
- `r_figures/`: R scripts for the final manuscript and Supplementary Information figures.

## Data Availability

The raw SEEG recordings analysed in this study are human-participant clinical neurophysiology data. They are not publicly redistributable through GitHub because they may contain sensitive participant information and were collected under ethics approval, consent terms and clinical data-use conditions that restrict open redistribution. Qualified researchers seeking access should apply through the responsible clinical institution or its data-access/ethics review process. Access would require appropriate ethics approval, a data-use agreement and approval by the responsible data controller.

The code expects approved users to place local input files under the example paths shown in each script, such as `data/` and `figures/`. These paths are placeholders and should be adapted to the user's approved data-access environment. Figure source tables and Supplementary Tables supporting the manuscript should be provided with the manuscript or a controlled repository record where permitted by ethics and data-use conditions.

## Notes

For the main empirical DMWA construction, the MATLAB scripts select one PLV threshold per seizure from the elbow of the density-threshold curve and apply that threshold to all windows from the same seizure. The Fig. 6 robustness scripts use a 50-setting grid spanning window lengths from 1 to 5 s and retained PLV edge fractions from 0.45 to 0.90.
