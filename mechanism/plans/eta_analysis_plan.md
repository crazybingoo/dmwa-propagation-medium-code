# Frozen analysis plan: conditional response interpretation of eta

Frozen on 2026-09-16 before computing new responses. Original manuscripts, figures, cached matrices and earlier reconstruction files remain read-only.

## Question and claim

Question: which response difference does eta quantify under an explicit, limited dynamical model?
Claim: for a nonnegative DMWA matrix W with R = sum(W)/m > 0 and rho(W) > 0, eta = R/rho(W) is the critical normalized gain in a homogeneous-leak linear feedback model. This is an operational mathematical interpretation, not independent physiological validation or proof of higher-order necessity.

## Model and endpoints

A = W/R; dx/dt = delta[-I + g A^T]x, delta = 1; W_ij is source i to target j. Initial condition x(0) = ones(m)/m. Mean outgoing normalized weight is exactly 1 for every matrix. Units of time are model units, not seconds.

Gains are fixed at 0.5, 0.8 and 0.9. No gain will be selected between example thresholds after seeing responses. An additional analytic gain curve on [0, 1.2] with step 0.005 describes the threshold continuously. Primary outputs: spectral growth rate s = g/eta - 1, critical gain g_c = eta, and asymptotic recovery time -1/s only if s < 0. The trajectory endpoint is total model activity sum(x). No p-value is assigned to a mathematical identity.

Numerical time grid: 0 to 80 inclusive, step 0.05. Fit natural-log total activity against time in the fixed intervals [20,40] and [40,80]. Preserve both slope estimates. A trajectory is called asymptotically resolved only when both slopes differ from the analytic rate by at most 0.0005 absolute model-time^-1 (one order smaller than the approximately 1% eta-scale effect near g = 0.5-0.9). Also report the slope difference between intervals and the relative error. Nonconvergence is retained and counted, without outcome-based extensions or omissions.

Use sparse expm_multiply for trajectories; avoid full dense matrix exponentials at every time point. For an independent solver check, DOP853 integrates y' = g(A^T - rho(A)I)y with x = exp(s t)y, rtol 1e-11 and atol 1e-13. This algebraically equivalent shift prevents absolute-tolerance domination of small decayed signals. Check relative trajectory error <= 1e-7 and slope error between numerical solvers <= 1e-7. For examples additionally compare with direct, unshifted ODE integration while trajectory amplitude is above 1e-8. Eigenvalue residual <= 1e-9 and cached R/rho/eta discrepancies are reported. Numerical failures are reported rather than silently suppressed.

## Fixed legal DMWA examples

Two 5-region scaffolds are fixed before responses: chain-attached edges {(1,2),(1,3),(2,3),(3,4),(4,5)}; hub-attached edges {(1,2),(1,3),(2,3),(1,4),(1,5)}. Retain every pair edge and every closed triangle (the single triangle {1,2,3}); no arbitrary positive matrix is presented as DMWA. Build W by the manuscript Coverage/degree-bias rule. Raw W may have different R; response operators use W/R, fixing effective R to 1 while preserving eta. Preserve incidence, scaffold edges, raw W and normalized W. These examples share 5 regions, 5 pairs, 1 triangle, and 6 hyperedge states. They are explanatory cases, not empirical observations or response-selected samples.

## Patient matrix substitution

Reuse the previously cached 24 seizure DMWA series and their original five-stage labels. For each seizure and stage choose the window index nearest (minimum eligible index + maximum eligible index)/2; break ties toward the earlier index. With uniformly spaced original windows this is the midpoint among available window centres. Maximum 120 matrices; missing stages stay missing. Selection uses only seizure/stage/index, not eta or response.

Do not reload or reprocess raw SEEG. Use anonymous patient/seizure identifiers. Recompute R and rho from each selected cached W and check against the saved window metrics. Retain all selected W, including a nonpositive-R/rho or numerical failure record. For valid W simulate every fixed gain and report all failures. Patient matrices are model substitutions, not 120 independent patients, interventions, or out-of-sample physiology. Descriptive summaries identify 24 seizures and 14 patients; no stage-inference claim is added from this deterministic transformation of eta.

## Boundary audits

For both legal examples and selected patient matrices verify eta(cW)=eta(W) with c in {0.5,2}, eta(W^T)=eta(W), and model threshold s(g=eta)=0. Explain that W/R is scale-invariant while raw fixed-coupling dynamics is not. Transposition preserves eta and the spectral rate, not every localized source-target response.

For the two legal examples audit heterogeneous leakage D = diag(linspace(0.5,1.5,m)) and its reversed assignment, keeping mean leakage 1. Its threshold is 1/rho(D^-1 A^T), not generally eta. Preserve departures even if small. Report finite-time activity/peak descriptively; do not claim eta determines nonnormal transient amplification, local transfer, nonlinear spread, or a conserved Laplacian diffusion process.

## Figure and outputs

Four panels: A explicit model and conditional identity; B the two legal same-normalized-R examples; C response curves for all three frozen gains; D patient model substitution and finite-time agreement/failures. R only for figures; source root's manuscript_style.R. Arial main text 10.5 pt; bold uppercase panel labels 12.5 pt; manuscript blue/orange/green/grey palette; SVG, Cairo PDF, 600 dpi PNG/TIFF and preview. Preserve source CSVs, provenance hashes, numerical QA, alignment/collision QA and English Results/legend/Methods. Analyses use Python without plotting. No claim of physiological validation, predictive superiority or higher-order necessity may be introduced.
