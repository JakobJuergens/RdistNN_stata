{smcl}
{* *! version 0.2.1  2026-10-17}{...}
{viewerjumpto "Syntax" "clate_dnn_diagnostics##syntax"}{...}
{viewerjumpto "Interpretation" "clate_dnn_diagnostics##interpretation"}{...}
{viewerjumpto "Options" "clate_dnn_diagnostics##options"}{...}
{viewerjumpto "Stored results and schemas" "clate_dnn_diagnostics##results"}{...}
{viewerjumpto "Examples" "clate_dnn_diagnostics##examples"}{...}

{title:Title}

{phang}{bf:clate_dnn_diagnostics} {hline 2} Inspect CLATE localization and scientific diagnostics

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}{cmd:clate_dnn_diagnostics} [{cmd:,} {opt targets(numlist)}
{opt candidates(numlist)} {opt allcandidates} {opt detail}
{opt saving(filename)} {opt replace}]

{pstd}
Requires active {help clate_dnn} estimates. Returns stored diagnostic tables
without refitting or changing the active estimates or data. See {help rdistnn}
for the package overview and {help clate_dnn_plot} for figures.

{marker interpretation}{...}
{title:Interpretation and availability}

{pstd}
ess=1/sum(w^2) measures concentration of normalized weights. max_weight and
top-1/top-5/top-10 shares measure the largest raw weights. For exact top-k DNN,
the unaveraged rank table has positive support through n-s+k (n-s+1 at k=1).
n_nonzero_weight counts every represented nonzero weight, including tiny
positive weights. Larger s tightens localization.

{pstd}
The support formula describes the unaveraged rank table. Averaging assigned
weights across exact-distance groups can extend nonzero observation weights to
the end of a group crossing that edge. Use actual returned weights for realized
support and concentration; rank-only tables do not contain target-specific groups.

{pstd}
radius_rank_k and radius_rank_s are distances to those nearest ranks under
the fitted localization metric (L1, L2, scaled L1/L2 or Mahalanobis), not positive-support radii. Weighted distance
mean/median describe where the mass lies. Current Stata fits retain these
diagnostics automatically; there is no W/dx retention option to specify.
A high ESS or favorable radius does not establish Gaussian inference,
negligible target bias, or credible nuisance/causal assumptions.
{pstd}
Arm mass and normalized within-arm ESS concern instrument A=1/A=0.
Within-arm ESS is (sum_arm w)^2/sum_arm(w^2), missing when arm mass is zero.
One positive-weight observation gives arm ESS one even with tiny arm mass.
Inspect arm mass alongside ESS when judging local support.
Overlap, arm-ESS and propensity flags are zero for paired-signal fits,
with unavailable values missing. Supplied IV-nuisance routes supply these fields.
Propensity summaries use clipped probabilities, which does not establish overlap.

{pstd}
physical_radius_90 is the ordinary Euclidean radius enclosing 90% of fitted
weight in original X coordinates, even with Mahalanobis ranking. It is distinct
from radius_rank_s. top5_score_share is the top-five share of absolute
w*(phi_y-theta*phi_d) contributions, distinct from raw-weight share and
jackknife pseudo-values. The share is missing when all contributions are zero,
even if score_contribution_available is one. Undefined ratios or nonfinite contributions have score_contribution_available=0.

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
coordinates are variables; row labels are not saved. Target scope appends the eight component columns after target_id and all five status columns; candidate scope saves the candidate table.
Data are restored even if saving fails; a report dataset is not a saved fit.{p_end}
{phang}{opt replace} permits overwriting the saved dataset.{p_end}

{marker results}{...}
{title:Stored results and column schemas}

{synoptset 24}{...}
{synopt:{cmd:r(diagnostics)}}full selected scientific/candidate matrix{p_end}
{synopt:{cmd:r(compact)}}compact displayed selection{p_end}
{synopt:{cmd:r(localization)}}common localization block; target scope only{p_end}
{synopt:{cmd:r(components)}}selected components; target scope only{p_end}
{synopt:{cmd:r(status)}}selected five status columns; target scope only{p_end}
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
{cmd:overlap_available instrument_1_mass instrument_0_mass instrument_1_share
 propensity_available weighted_mean_propensity weighted_propensity_q10
 weighted_propensity_q90 physical_radius_90 arm_ess_available
 instrument_1_ess instrument_0_ess score_contribution_available top5_score_share}.
The x0_j coordinates follow the fitted xcols() order.
Row identities and target IDs survive selection.

{pstd}
Target r(compact) columns:
{cmd:target_id ess top5_weight_share physical_radius_90 instrument_1_ess instrument_0_ess top5_score_share}.

{pstd}
Candidate r(diagnostics), and the candidate saved dataset, have 13 columns:
{cmd:candidate_id computed_id is_default feasible reason status lambda s
 ess max_weight top_n top5_weight_share degenerate}.
computed_id identifies duplicate resolved orders sharing computation;
is_default marks the ordinary fit. status is 1 available or 0 infeasible.
reason is 0 feasible, 1 s<k, 2 s>=n with inference, or 3 s>n with nose.
lambda is missing for sgrid(); top_n=min(5,n).
ess is the minimum across targets; max_weight and top5_weight_share are
the corresponding maxima across targets. degenerate is one if any target
has ESS>=.9*n, a near-uniform-weight guide, not a test.
Skipped numerical fields remain missing. Candidate compact columns are
{cmd:candidate_id is_default feasible status lambda s ess max_weight
 top5_weight_share degenerate}.
Candidate scope has no localization, component, or status matrix.

{marker components}{...}
{title:CLATE components and computation status}

{pstd}
For target scope, r(components) has columns:
{cmd:target_id reduced_form first_stage numerator_se denominator_se
component_covariance component_correlation denominator_stat estimate}.
The covariance is on the estimator scale. Correlation requires two positive
component SEs; denominator_stat is signed first_stage/denominator_se and
requires a positive denominator SE. Otherwise these standardized fields are missing.
Component and status rows follow the same selected target order, including
undefined ratios. e(component_vcov) retains interleaved Y1,D1,Y2,D2 covariance.

{pstd}
r(status) has five columns: ratio_status (0 computed, 1 zero denominator,
2 numerical failure); negative_denominator (0/1); near_zero_denominator (0/1);
covariance_status (0 computed, 1 not computed, 2 zero score variance,
3 numerical failure); and wald_available (0/1).
The near-zero flag compares the nonzero first stage with the maximum absolute
input denominator signal at denominatortol(); it never regularizes a ratio.
Negative first stages retain a computed signed ratio. Zero denominators retain
components but have an undefined ratio. A zero ratio SE is retained but has no
Wald bounds when wald_available=0. These are pointwise diagnostics with no
weak-IV-uniform guarantee or automatic weak-instrument cutoff.

{pstd}
Underlying e(diagnostics) retains internal treated/control arm slots;
this command names them instrument_1/instrument_0 in displayed, returned,
and saved tables. Effective sample size is not a complier count.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. clate_dnn_diagnostics}{p_end}
{phang2}{cmd:. clate_dnn_diagnostics, targets(3 1) detail}{p_end}
{phang2}{cmd:. matrix report = r(diagnostics)}{p_end}
{phang2}{cmd:. clate_dnn_diagnostics, allcandidates}{p_end}
{phang2}{cmd:. clate_dnn_diagnostics, saving("clate-diagnostics.dta") replace}{p_end}

{pstd}
Copy r() results before another rclass command. Candidate plots use only
available candidates, while allcandidates here deliberately retains skipped rows.
For the fitted estimate table, replay {cmd:clate_dnn, level(90)} and inspect
r(summary) or r(table); diagnostics do not replace coefficient reporting.

{title:Also see}

{psee}{help rdistnn}, {help clate_dnn}, {help clate_dnn_plot}{p_end}
