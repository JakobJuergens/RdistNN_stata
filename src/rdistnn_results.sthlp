{smcl}
{* *! version 0.1.0  2026-10-18}{...}
{title:rdistnn_results — Extract fitted results without drawing a graph}
{pstd}
After {help np_dnn}, {help np_tdnn}, {help cate_dnn}, {help clate_dnn} or {help plqr_dnn}:
{p 8 12 2}{cmd:rdistnn_results} [{cmd:, type(}{it:type}{cmd:) targets(}{it:numlist}{cmd:)}
{cmd:level(}{it:#}{cmd:) saving(}{it:filename}{cmd:) replace}]{p_end}
{pstd}
The default type is estimates. Other types are vcov, diagnostics, weights,
components and status (CLATE), and influence (CATE/CLATE). targets() selects
original target IDs in the requested order. level() is a percentage.
The active data and e() are preserved, including when saving fails.
{pstd}
Weights and weight-based diagnostics use the final mean rank weights within
exact represented selected-metric distance groups, separately for each target.
These are reconstructed from validated saved ranks, distances and groups.
Realized support and concentration can differ from the rank table. Weight
reconstruction warns that inference with persistent distance ties is not
established. Schema 3 and exact_distance_block_averaging_v1 are required;
refit older models. See {help weights_dnn} for survivor redistribution.
{title:Returned results}
{pstd}
r(fit) and r(se_fit) are q by m matrices: one row per component, one column
per selected target (q=1 outside PLQR). r(vcov) retains joint target covariance,
in target-major component order. r(index) identifies the target and component
for each flattened column. r(table) has Stata's nine reporting rows;
r(summary) has estimate, std_error, z, pvalue, ci_lower and ci_upper.
{pstd}
r(data) is the selected numeric export. Estimates include target, component,
x1,...,xp, and the six summary columns. vcov is a long table with row_target,
row_component, column_target, column_component and covariance. Diagnostics,
components and status prepend target to the corresponding e() matrix columns.
Weights contain observation_id, target, rank, distance and weight.
Influence contains target, observation_id, rank, distance, weight, score and
contribution; CLATE appends phi_y and phi_d. Observation IDs refer to fitted rows
in the original dataset. saving() writes r(data) as a double-precision .dta file.
No graph is created. Copy r() results before another command replaces them.
{pstd}
r(quantity) is moment for PLQR, with r(theta_eval) separately retaining the
supplied coefficient vector. Its tests and intervals concern score moments,
not coefficients. Other fits return quantity=estimate. r(level) is a percentage.
Undefined CLATE targets remain in canonical tables; unavailable uncertainty
remains missing. Confidence limits are pointwise, not simultaneous.
{pstd}
Reporting requires current fitted state. Weights and influence use its numeric
model snapshot. Refit objects from older development versions.
Current snapshots are structurally validated before observation-level
extraction; refit an incomplete or inconsistent snapshot. Validation errors
preserve the active data and e(). Replaying stored numeric estimates does not
reconstruct an older model with current weights.
For grid extraction, first use {help rdistnn_candidate}; the existing
estimator-specific diagnostics commands retain all-candidate table exports.
{title:See also}
{pstd}{help rdistnn}, {help rdistnn_candidate}, {help rdistnn_predict}{p_end}
