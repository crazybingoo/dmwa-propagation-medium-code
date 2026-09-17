# Current figure and code map

Numbering follows the 17 September 2026 manuscript and SM. Current panel letters are uppercase. Historical script basenames and intermediate filenames remain traceable to the original analysis.

| Current figure | Question / evidence | Workbook | R entry point | Numerical analysis |
|---|---|---|---|---|
| Fig. 1 | How are overlapping groups represented by DMWA? | Schematic; no numerical sheet | No numerical plotting script | `matlab/2_Fig2/` construction |
| Fig. 2 | Does eta change across seizure stages? | `Fig2` | `r_figures/Fig_2/make_fig2_effective_refractive_index_redraw_font10p5_label12p5_review.R` | `matlab/2_Fig2/` |
| Fig. 3 | Does mean coupling increase proportionally more than spectral radius? | `Fig3` | `r_figures/Fig_3/render_fig3_delta_log_notation.R` | `matlab/2_Fig3/` |
| Fig. 4 | Does support organization contribute beyond density and weight assignment? | `Fig4` | `r_figures/Fig_4/render_fig4_control_label_nouns.R` | `matlab/2_Fig4/` |
| Fig. 5 | How is resource–spectral burden distributed? | `Fig5` | `r_figures/Fig_5/make_fig5_eta_burden_font10p5_label12p5_legends_below.R` | `matlab/2_Fig5_size_adjusted/` |
| Fig. 6 | Which group context gives exact aggregated responses? | `Fig6` | `r_figures/Fig_6/make_fig6_group_context.R` | `mechanism/code/w_mechanism_analysis.py` |
| Fig. 7 | What feedback threshold and recovery scale does eta describe? | `Fig7` | `r_figures/Fig_7/make_fig7_feedback_response.R` | `mechanism/code/eta_response_analysis.py` |
| Fig. S1 | What are the individual seizure trajectories? | `SuppFig1` | `r_figures/Supplementary_Fig_1/make_Supplementary_Fig_2_eta_trajectories.R` | Derived eta time series |
| Fig. S2 | Does eta change at patient level, and is the effect sensitive to leaving one patient out? | `SuppFig2` | `r_figures/Supplementary_Fig_2/make_Supplementary_Fig_4_patient_sensitivity.R` | Primary eta stage inference on 14 patients |
| Fig. S3 | How do control effects vary across seizures? | `SuppFig3` | `r_figures/Supplementary_Fig_3/make_Supplementary_Fig_3_null_control_delta_eta.R` | `matlab/2_Fig4/` |
| Fig. S4 | How was the PLV threshold selected? | `SuppFig4` | `r_figures/Supplementary_Fig_4/make_Supplementary_Fig_1_threshold_elbow.R` | `select_seizure_plv_elbow_threshold.m` |
| Fig. S5 | How sensitive is eta to window length and edge fraction? | `SuppFig5` | `r_figures/Supplementary_Fig_5/render_fig6_retained_fraction_d_title_bold.R` | `matlab/2_Fig6/` |
| Fig. S6 | What changes when weighting components are removed? | `SuppFig6` | `r_figures/Supplementary_Fig_6/render_fig7c_degree_bias_spacing.R` | `matlab/2_Fig7/` |
| Fig. S7 | What changes when group order is restricted? | `SuppFig7` | `r_figures/Supplementary_Fig_7/render_fig8c_order_labels.R` | `matlab/2_Fig8/` |
| Fig. S8 | How do eta, synchrony and density behave in simulations? | `SuppFig8` | `r_figures/Supplementary_Fig_8/draw_fig9_hr_mechanism_final.R` | `matlab/2_Fig9_HR_eta_simulation/` and `matlab/2_Fig9_eta_mechanism/` |

## Numbering changes and retained provenance

Old main Figs. 6, 7, 8 and 9 are now Figs. S5, S6, S7 and S8. Current Figs. 6 and 7 are new mechanism analyses.

The old supplementary R folders for threshold selection, trajectories and patient sensitivity were numbered 1, 2 and 4. They are now in folders S4, S1 and S2, matching the workbook and SM. The control figure remains S3.

Older R scripts can emit historical filenames, panel typography and draft narrative text. The current manuscript/SM takes precedence over automatically generated draft legends. Figs. 6–7 have portable, directly tested final-layout scripts. Not every legacy analysis can be rerun without restricted or undistributed intermediates.
