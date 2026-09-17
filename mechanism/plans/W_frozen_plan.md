# Frozen W mechanism analysis plan

Status: frozen before running this new structural analysis. This is a conditional mathematical analysis of the existing W construction, not an independent neurophysiological validation. Earlier simulation outcomes are known. No old figures, reconstruction outputs or manuscript files will be changed.

## Scientific question and representation

Does forgetting hyperedge order and overlap-neighbourhood context prevent node-level closure of a specified W-defined group response model? H has hyperedges in rows and original nodes in columns. Every scaffold contributes all its edges and all its closed triangles, without arbitrary independent hyperedge sampling. W_ij is the source-i to target-j weight: overlap/source size, multiplied by D_j/(D_i+D_j) only for same-size groups. Group states are columns. J = -I + (g/R) W^T, with g = 0.5, local decay = 1 and R = mean outgoing W weight. Time is fixed from 0 to 10 in 101 equal steps. Stability is checked and reported, not imposed by changing g.

C = H^T D_size^-1 splits group-state mass equally among members. Independent nonnegative group states are admissible model states. They are not claimed to be reachable through a single node-input lift, and they are not interchangeable with group states constrained to be instantaneous member averages. In particular, a pair differing in ker(C) is generally unavailable if all group states are restricted to range(C^T).

## Frozen scaffold set

Three toys: (1) triangle ABC; (2) triangle ABC with branches AD and AE; (3) triangles ABC and BCD sharing BC. There are no other scaffold edges; all implied closed triangles are included. The complete H/W stays identical across the two initial states. Toy initial state 1 puts mass 1 on ABC; state 2 puts mass 1/3 on each of AB, AC and BC. Their node projections and total initial masses agree exactly.

Formal set: 100 total connected, 12-node scaffolds, seeds 930000–930099. Families cycle through random, core-periphery, modular, giving 34/33/33 instances. Each includes a fixed 12-node ring before adding independently sampled undirected edges. Random extra-edge probability is 0.20. Core-periphery uses the first four nodes as core and probabilities 0.70/0.25/0.08 for core-core/core-periphery/periphery-periphery. Modular uses three blocks of four nodes with probabilities 0.65 within and 0.05 between blocks. Edge presence, not new artificial weighted hyperedges, defines the scaffold. No scaffold is regenerated or discarded based on a result.

For each formal scaffold, z_bar = 1/m. A fixed independent normal vector is projected into ker(C) using geometry only, without W or response optimization. The direction is scaled so its largest absolute entry is 0.8/m, and z_plus/minus = z_bar ± d. Both remain nonnegative, with the same C projection and total mass. If the projection is injective or the projected vector is numerically zero, that instance is explicitly marked as having no hidden-state pair; no replacement graph is drawn. Closure statistics remain included.

## Comparisons and analytical acceptance lines

Three weighting rules: Coverage-only; constant-bias (same order 0.5, cross order 1); full DMWA. Three retained-state projections: node only, node × order, and node × exact(order, overlap degree). Structural closure error is ||P W^T (I-P^+P)||_F / ||P W^T||_F. Zero denominator is recorded as an identically silent projected operator. P^+ is the SVD pseudoinverse. The unconstrained best instantaneous linear closure is B = P J P^+. It is not claimed to be the best finite-time predictor.

Mandatory analytical checks, tolerance 1e-10: Coverage-only closes at node level; constant-bias closes at node × order; full DMWA closes at node × exact(order,degree). These identities are tested on every toy and formal scaffold. Full DMWA need not fail the simpler projections; all outcomes are retained. Exact context recovery may use more coordinates than the original group space. Coordinate count, active-coordinate count and rank are reported; recovery is not automatically compression.

## Fixed response endpoints

The full group model is integrated exactly as a matrix exponential, as are the candidate projected closures. For the same-C initial pair, report separately: (a) time-integrated absolute total-mass separation divided by integrated mean total mass; (b) mean total-variation distance between normalized node profiles over 0–10. Both include all 101 times. This distinguishes total amplification from spatial redistribution. At fixed t=1 and t=10, values are also retained; the maximum is diagnostic only.

For each projected closure, reconstruct the node readout by summing contexts. The primary approximation error is sqrt(integrated squared node-response error / integrated squared true node response), pooled equally over the two initial states. Exact context restoration must agree below 1e-8 relative error (apart from explicitly diagnosed ill-conditioned cases). Negative approximating states, unstable operators, condition numbers and failures are recorded rather than silently clipped. Node-profile TV is computed only from the two nonnegative full-model trajectories, not from possibly signed least-squares approximations.

## Decision rules and figure contract

If full W closes at node level, there is no evidence of loss in that case. If order alone restores closure, interpretation is limited to order context. Coverage-only failure or full order-degree recovery failure is an implementation/theorem problem, not new physiology. If finite-response discrepancies are small, report their magnitude without declaring practical necessity from a nonzero algebraic residual. If only total response changes, do not claim changed spatial propagation. No p-value threshold is used to decide which graphs or panels to show.

The four-panel figure will explain: A fixed admissible same-node initial states on a legal toy; B closure across weighting and retained context; C total versus spatial response loss in all formal instances; D restoration of finite-time node responses by added context. Plotting uses R only, the supplied manuscript style, uppercase labels, editable SVG/PDF and 600-dpi PNG/TIFF. Figure text and the Results/legend will state the W-defined model boundary. An exact context representation remains an explicit recovery control; full W is not claimed to be the only sufficient representation.
