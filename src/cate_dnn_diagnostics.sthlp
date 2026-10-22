{smcl}
{* *! version 0.2.1  2026-10-17}{...}
{viewerjumpto "Syntax" "cate_dnn_diagnostics##syntax"}{...}
{viewerjumpto "Interpretation" "cate_dnn_diagnostics##interpretation"}{...}
{viewerjumpto "Options" "cate_dnn_diagnostics##options"}{...}
{viewerjumpto "Stored results and schemas" "cate_dnn_diagnostics##results"}{...}
{viewerjumpto "Examples" "cate_dnn_diagnostics##examples"}{...}

{title:Title}

{phang}{bf:cate_dnn_diagnostics} {hline 2} Inspect CATE localization and scientific diagnostics

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}{cmd:cate_dnn_diagnostics} [{cmd:,} {opt targets(numlist)}
{opt candidates(numlist)} {opt allcandidates} {opt detail}
{opt saving(filename)} {opt replace}]

{pstd}
Requires active {help cate_dnn} estimates. Returns stored diagnostic tables
without refitting or changing the active estimates or data. See {help rdistnn}
for the package overview and {help cate_dnn_plot} for figures.

{marker interpretation}{...}
{title:Interpretation and availability}

{pstd}
ess=1/sum(w^2) measures concentration of normalized weights. max_weight and
top-1/top-5/top-10 shares measure the largest raw weights. For exact top-k DNN,
the unaveraged rank table has positive support through n-s+k (n-s+1 at k=1).
n_nonzero_weight counts abs(w)>sqrt(machine epsilon), so extremely small
positive weights can be excluded. Larger s tightens localization.

{pstd}
The support formula describes the unaveraged rank table. Averaging across
identical covariate rows can extend nonzero observation weights to the end of a
block crossing that edge. Use the actual returned weights for realized support
and concentration; rank-only tables do not contain the duplicate groups.

{pstd}
radius_rank_k and radius_rank_s are distances to those nearest ranks under
the fitted localization metric (L1, L2, scaled L1/L2 or Mahalanobis), not positive-support radii. Weighted distance
mean/median describe where the mass lies. Current Stata fits retain these
diagnostics automatically; there is no W/dx retention option to specify.
A high ESS or favorable radius does not establish Gaussian inference,
negligible target bias, or credible nuisance/causal assumptions.
{pstd}
Arm mass and normalized within-arm ESS concern treatment D=1/D=0.
Within-arm ESS is (sum_arm w)^2/sum_arm(w^2), missing when arm mass is zero.
Overlap, arm-ESS and propensity flags are zero for precomputed pseudo-outcome fits,
with unavailable values missing. AIPW routes supply these fields.
Propensity summaries use clipped probabilities, which does not establish overlap.

{pstd}
physical_radius_90 is the ordinary Euclidean radius enclosing 90% of fitted
weight in original X coordinates, even with Mahalanobis ranking. It is distinct
from radius_rank_s. top5_score_share is the top-five share of absolute
w*(rho-tau) contributions, distinct from raw-weight share and
jackknife pseudo-values. The share is missing when all contributions are zero,
even if score_contribution_available is one. These diagnostics are also available with nose.

{marker options}{...}
{pstd}
Physical radii use scaled Euclidean norms to avoid intermediate square overflow
or underflow. Weighted quantiles select the first observation whose compensated
prefix divided by compensated total reaches the probability, in stable value
order, without a boundary tolerance.

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
{cmd:overlap_available treated_mass control_mass treated_share
 propensity_available weighted_mean_propensity weighted_propensity_q10
 weighted_propensity_q90 physical_radius_90 arm_ess_available
 treated_ess control_ess score_contribution_available top5_score_share}.
The x0_j coordinates follow the fitted xcols() order.
Row identities and target IDs survive selection.

{pstd}
Target r(compact) columns:
{cmd:target_id ess top5_weight_share physical_radius_90 treated_ess control_ess top5_score_share}.

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

{phang2}{cmd:. cate_dnn_diagnostics}{p_end}
{phang2}{cmd:. cate_dnn_diagnostics, targets(3 1) detail}{p_end}
{phang2}{cmd:. matrix report = r(diagnostics)}{p_end}
{phang2}{cmd:. cate_dnn_diagnostics, allcandidates}{p_end}
{phang2}{cmd:. cate_dnn_diagnostics, saving("cate-diagnostics.dta") replace}{p_end}

{pstd}
Copy r() results before another rclass command. Candidate plots use only
available candidates, while allcandidates here deliberately retains skipped rows.
For the fitted estimate table, replay {cmd:cate_dnn, level(90)} and inspect
r(summary) or r(table); diagnostics do not replace coefficient reporting.

{title:Also see}

{psee}{help rdistnn}, {help cate_dnn}, {help cate_dnn_plot}{p_end}
