# Group context and feedback response analyses

## Matrix and model conventions

`W[i,j]` couples source group `i` to target group `j`. Incidence matrix `H` has groups in rows and regions in columns. Every retained pair and closed triangle is included. Off-diagonal weights are overlap divided by source-group size. Same-order overlaps also receive `degree[j] / (degree[i] + degree[j])`. Diagonal weights are zero.

Group activity is a column vector:

```text
R = sum(W) / number_of_groups
eta = R / spectral_radius(W)
J = -I + (g/R) W.T
```

Time is in model units. The group-to-node projection shares each group's activity equally among its members. Initial unit mass normalizes activity; it is not a conserved physical mass.

## Fig. 6

The formal set has 100 connected 12-node scaffolds: 34 random, 33 core-periphery and 33 modular. Seeds are 930000–930099. Each includes a ring. Extra-edge probabilities are 0.20 for random graphs, 0.70/0.25/0.08 for core–core/core–periphery/periphery–periphery, and 0.65/0.05 within/between modules. Three fixed illustrative examples are separate.

Paired nonnegative group states have identical node projections. Their hidden direction uses incidence geometry and seed `1000000 + scaffold_seed`, independently of response outcomes. Models use gain 0.5, unit leakage and 101 times over 0–10. Pseudoinverse relative cutoff is 1e-12; exact structural residual tolerance is 1e-10.

Coverage-only closes at node level; constant same-order bias closes with order retained; full DMWA closes with order and exact overlap degree retained. These are sufficient representations, not claims of a minimal or compressed representation. Finite-time errors are reported separately from algebraic exactness. Independent group states are a model assumption.

`data/w/` includes scaffold, group, closure and response records. Full DMWA formal-scaffold rows underlie panels C and D. Toy examples are excluded from formal medians.

## Fig. 7

Two five-region structures each have five edges and one triangle. Coupling is normalized by mean outgoing weight. Fixed gains are 0.5, 0.8 and 0.9. With unit uniform leakage, dominant growth rate is `s = g/eta - 1`, critical gain is `eta`, and recovery time is `-1/s` only when `s < 0`.

One window is selected per seizure/stage, nearest the midpoint of the available window-index range, with an earlier-index tie break. The release supplies **120 derived matrices from 24 seizures in 14 patients**, selected-window records and cached resource/spectral values. Raw SEEG is not included.

Sparse matrix-exponential multiplication and DOP853 solve the same shifted system at 1,601 times over 0–80. DOP853 tolerances are 1e-11 relative and 1e-13 absolute. Log-total-activity slopes in intervals 20–40 and 40–80 must both agree with the analytic rate within 0.0005 per model-time unit. Panel D uses 40–80. Solver agreement tolerance is 1e-7. Failures are recorded and cause a nonzero exit.

The 360 responses are model evaluations, not independent patients or observed neural responses. Scaling and transposition preserve eta. Heterogeneous leakage changes the critical gain; that audit is included.

## Portability and provenance

`plans/` preserves original dated internal design records, not external preregistration. The public eta script loads saved anonymous matrices instead of a private MATLAB cache and identifier lookup. Equations and numerical settings are unchanged.

Patient labels now match Table S1. Grouping and numerical values are unchanged. `input_manifest.json` records matrix hashes. Local clinical identifiers and workstation paths are excluded.

R scripts read the distributed CSVs. If numerical output is redirected using `DMWA_OUTPUT_ROOT`, point the R input variables to that root's `data/w` and `data/eta` to plot recomputed results.
