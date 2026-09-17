# Validation of the September 2026 update

The checks below concern the current repository update. The original raw-SEEG analyses were not rerun.

## Numerical reproduction

- Ran the portable Fig. 6 analysis on all 100 formal scaffolds and three fixed examples. No failed, unstable or omitted full-model cases were reported.
- Maximum structural residual in the exact controls: approximately 4.95e-15. Maximum relative response error in the exact controls: approximately 4.21e-16.
- Ran the portable Fig. 7 analysis from the 120 supplied sparse matrices, producing all 360 responses. Both fixed fitting intervals passed the 0.0005 criterion in every response.
- Largest fitted-versus-analytic growth-rate discrepancy across both fitting intervals in this rerun: approximately 1.70e-5 per model-time unit. Largest relative trajectory difference between the two solvers: approximately 2.99e-9. The archived manuscript calculation reported 3.35e-9; both satisfy the 1e-7 solver criterion. Small numerical differences are retained in the validation record, not substituted into manuscript results.
- Compared complete closure, trajectory-error, separation and rate tables with the distributed results. The manuscript result tables retain their original values.

## Figures and workbook

- Executed both new R plotting scripts from the public package; inspected the resulting Fig. 6 and Fig. 7 previews. Data counts are 100 formal structures for Fig. 6 and 120 matrices/360 responses for Fig. 7.
- Updated workbook numbering and checked all retained empirical/simulation numeric cells against the previous public workbook.
- Checked new mechanism workbook blocks against their CSV values, preserving small residuals as numerical values rather than rounded zeros.
- Checked Data_Map ranges against data-row counts. Tables S2 and S3 were transcribed from the current SM.
- Matched the mechanism anonymous patient labels to Table S1 while verifying that the seizure-to-patient partition is unchanged.
- Exported all 55 workbook source blocks and parsed all 20 R scripts. Updated metadata to identify patient-level primary inference and the Fig. 3 stage-mean log-change convention.

The machine-readable record is [validation.json](validation.json). The two numerical entry points use NumPy, SciPy, pandas and threadpoolctl. Plotting uses R and the packages listed in the top-level README.

## Reproduction boundary

The Fig. 6 synthetic analyses and Fig. 7 derived-matrix responses can be rerun from this repository alone. Regenerating the patient DMWA matrices from original SEEG still requires controlled data access. Some original MATLAB/R pipelines also require undistributed intermediate tables; the public workbook records the panel-level numbers, not all such intermediates. Analytical/numerical agreement does not establish independent physiological validation.
