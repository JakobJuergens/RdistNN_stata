{smcl}
{* *! version 0.2.1  2026-10-17}{...}
{viewerjumpto "Syntax" "np_dnn_diagnostics##syntax"}{...}
{viewerjumpto "Interpretation" "np_dnn_diagnostics##interpretation"}{...}
{viewerjumpto "Options" "np_dnn_diagnostics##options"}{...}
{viewerjumpto "Stored results and schemas" "np_dnn_diagnostics##results"}{...}
{viewerjumpto "Examples" "np_dnn_diagnostics##examples"}{...}

{title:Title}

{phang}{bf:np_dnn_diagnostics} {hline 2} Inspect nonparametric regression localization and scientific diagnostics

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}{cmd:np_dnn_diagnostics} [{cmd:,} {opt targets(numlist)}
{opt candidates(numlist)} {opt allcandidates} {opt detail}
{opt saving(filename)} {opt replace}]

{pstd}
Requires active {help np_dnn} estimates. Returns stored diagnostic tables
without refitting or changing the active estimates or data. See {help rdistnn}
for the package overview and {help np_dnn_plot} for figures.

{marker interpretation}{...}
{title:Interpretation and availability}

{pstd}
ess=1/sum(w^2) measures concentration of normalized weights. max_weight and
top-1/top-5/top-10 shares measure the largest raw weights. For exact top-k DNN,
mathematical positive support has n-s+k observations (n-s+1 at k=1).
n_nonzero_weight counts abs(w)>sqrt(machine epsilon), so extremely small
positive weights can be excluded. Larger s tightens localization.

{pstd}
radius_rank_k and radius_rank_s are distances to those nearest ranks under
the fitted L2/Mahalanobis metric, not positive-support radii. Weighted distance
mean/median describe where the mass lies. Current Stata fits retain these
diagnostics automatically; there is no W/dx retention option to specify.
A high ESS or favorable radius does not establish Gaussian inference,
negligible target bias, or credible nuisance/causal assumptions.
{pstd}
NP summaries describe the weighted outcome distribution. weighted_mean_y equals
 the estimate for normalized weights. weighted_sd_y is the square root of the
 weighted second central moment divided by sum(w), without a sample correction.
 Weighted q10/q50/q90 select the first observed outcome whose cumulative
 normalized weight reaches the probability; they are not interpolated.
 Outcome quantile ranges are not confidence bounds for the mean.

{marker options}{...}
{title:Options}

{phang}{opt targets(numlist)} selects original e(targets) row numbers,
positive integers in 1..e(m), in supplied order. Default is all targets.
It cannot accompany candidate selection.{p_end}
{phang}{opt candidates(numlist)} selects original tuning-candidate rows,
in supplied order, after tuningmode(grid). Infeasible rows remain reportable.{p_end}
{phang}{opt allcandidates} selects every stored candidate, including
infeasible rows. It cannot accompany candidates().{p_end}
{phang}{opt detail} displays all columns; the default compact display does
not reduce r(diagnostics) or the saved table.{p_end}
{phang}{opt saving(filename)} saves the full selected table as a Stata
dataset with double variables named after matrix columns. Explicit IDs and
coordinates are variables; row labels are not saved. The saved table has the same columns as r(diagnostics).
Data are restored even if saving fails; a report dataset is not a saved fit.{p_end}
{phang}{opt replace} permits overwriting the saved dataset.{p_end}

{marker results}{...}
{title:Stored results and column schemas}

{synoptset 24}{...}
{synopt:{cmd:r(diagnostics)}}full selected scientific/candidate matrix{p_end}
{synopt:{cmd:r(compact)}}compact displayed selection{p_end}
{synopt:{cmd:r(localization)}}common localization block; target scope only{p_end}
{synopt:{cmd:r(n_rows)}}number of selected rows{p_end}
{synopt:{cmd:r(scope)}}{cmd:targets} or {cmd:candidates}{p_end}
{synopt:{cmd:r(selected)}}original IDs in returned order{p_end}
{synopt:{cmd:r(saved)}}filename if saving() supplied{p_end}

{pstd}
Target r(localization) columns, in order:
{cmd:target_id x0_1 ... x0_p ess max_weight top1_weight_share
 top5_weight_share top10_weight_share n_nonzero_weight distance_available
 radius_rank_k radius_rank_s weighted_mean_distance weighted_median_distance}.
Target r(diagnostics) starts with this block and appends:
{cmd:outcome_available weighted_mean_y weighted_sd_y weighted_y_q10
 weighted_y_q50 weighted_y_q90}.
The x0_j coordinates follow the fitted xcols() order.
Row identities and target IDs survive selection.

{pstd}
Target r(compact) columns:
{cmd:target_id ess max_weight top5_weight_share weighted_mean_y weighted_sd_y}.

{pstd}
Candidate r(diagnostics), and the candidate saved dataset, have 13 columns:
{cmd:candidate_id computed_id is_default feasible reason status lambda s
 ess max_weight top_n top5_weight_share degenerate}.
computed_id identifies duplicate resolved orders sharing computation;
is_default marks the ordinary fit. status is 1 available or 0 infeasible.
reason is 0 feasible, 1 s<k, 2 s>=n with inference, or 3 s>n with nose.
lambda is missing for sgrid(); top_n=min(5,n).
degenerate means ess>=.9*n, a near-uniform-weight guide, not a test.
Skipped numerical fields remain missing. Candidate compact columns are
{cmd:candidate_id is_default feasible status lambda s ess max_weight
 top5_weight_share degenerate}.
Candidate scope has no localization, component, or status matrix.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. np_dnn_diagnostics}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, targets(3 1) detail}{p_end}
{phang2}{cmd:. matrix report = r(diagnostics)}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, allcandidates}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, saving("np-diagnostics.dta") replace}{p_end}

{pstd}
Copy r() results before another rclass command. Candidate plots use only
available candidates, while allcandidates here deliberately retains skipped rows.
For the fitted estimate table, replay {cmd:np_dnn, level(90)} and inspect
r(summary) or r(table); diagnostics do not replace coefficient reporting.

{title:Also see}

{psee}{help rdistnn}, {help np_dnn}, {help np_dnn_plot}{p_end}
