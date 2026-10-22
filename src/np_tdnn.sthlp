{smcl}
{* *! version 0.1.0 18oct2026}{...}
{title:np_tdnn — Two-scale nonparametric regression}

{p 4 4 2}
{cmd:np_tdnn} {it:y} {ifin}, {cmd:xcols(}{it:varlist}{cmd:)}
{cmd:x0(}{it:numlist}{cmd:)} [{cmd:s1(}{it:#}{cmd:)} {cmd:s2(}{it:#}{cmd:)}
{cmd:dim(}{it:#}{cmd:)} {cmd:method(exact|asymp)} {cmd:nonormalize}
{cmd:distance(}{it:metric}{cmd:)} {cmd:scalerule(}{it:rule}{cmd:)}
{cmd:scalevec(}{it:numlist}{cmd:)} {cmd:scaleeps(}{it:#}{cmd:)}
{cmd:nose} {cmd:level(}{it:#}{cmd:)}]

{p 4 4 2}
Use {cmd:x0matrix(}{it:M}{cmd:)} instead of {cmd:x0()} for multiple target
rows. Defaults are {cmd:s1(1)},
{cmd:s2(2)}, and {cmd:dim()} equal to the number of covariates.
Require 1 <= s1 < s2 <= n; jackknife inference requires s2 < n.

{title:Signed two-scale weights}

{p 4 4 2}
The estimate is the weighted outcome sum using w1*a(s1)+w2*a(s2), where
w1=1/(1-(s1/s2)^(-2/dim)) and w2=1-w1. Individual weights may be negative.
Exact tables use the k=1 DNN formula and have support through n-s1+1.
Asymptotic tables use Beta(1,s) CDF increments; {cmd:nonormalize} disables
their numerical normalization. Exact tables always retain the R reference's
normalization safeguards. Outcomes and distance geometry remain frozen under
common delete-one resamples. Joint covariance retains cross-target dependence.

{pstd}
The support formula describes the unaveraged rank table. Averaging across
identical covariate rows can extend nonzero observation weights to the end of a
block crossing that edge. Use the actual returned weights for realized support
and concentration; rank-only tables do not contain the duplicate groups.

{p 4 4 2}
Distances are {cmd:l2}, {cmd:l1}, {cmd:l2_scaled}, {cmd:l1_scaled}, and
{cmd:mahalanobis}. For scaled metrics, {cmd:scalerule()} is {cmd:sd},
{cmd:mad}, {cmd:iqr}, {cmd:range}, or {cmd:custom}; default is {cmd:sd}.
{cmd:scalevec()} supplies positive scales with {cmd:scalerule(custom)}.
The default positive scale floor is 1e-12. Fitted geometry is saved for reuse.

{title:Inference}

{pstd}
The documented sampling framework uses independent observations; cluster
and survey-design inference are not supplied. Exactly identical raw covariate
rows receive their block's average signed weight. Each deletion removes the
actual observation before averaging the n-1 table within the remaining block.
The resulting estimates and covariance are invariant to aligned permutations
up to numerical precision. A warning identifies this modification as outside
the paper's theory. Different vectors at equal distance are not grouped.


{p 4 4 2}
Intervals are pointwise normal intervals based on delete-one jackknife
covariance. Signed bias correction alone does not establish negligible bias,
Gaussianity, or valid coverage. {cmd:nose} leaves uncertainty unavailable.
Replay with {cmd:np_tdnn, level(90)}. Signed weights do not define probability
distributions, so TDNN does not report weighted outcome quantiles or variances.

{title:Saved results}

{pstd}
{cmd:e(depvar)} records the outcome. {cmd:e(vce)} is signed_rank_shift or
none; {cmd:e(vcetype)} is Jackknife, or empty with {cmd:nose}.

{p 4 4 2}
{cmd:e(estimates)}, {cmd:e(se_targets)}, {cmd:e(estimate_vcov)},
{cmd:e(targets)}, and {cmd:e(diagnostics)} retain all targets.
{cmd:e(s1)}, {cmd:e(s2)}, {cmd:e(dim)}, {cmd:e(normalize)}, and the
{cmd:e(model_*)} numeric state support saved-fit reuse.
Compatibility fields {cmd:e(neff_targets)}, {cmd:e(neff_min)},
{cmd:e(neff_max)} and single-target {cmd:e(neff)} retain inverse squared-weight
concentration. With signed weights these fields do not measure a
probability-based effective sample size.
{cmd:r(summary)} has estimate, std_error, z, pvalue, ci_lower, ci_upper;
{cmd:r(table)} is the coefficient table. Copy r() results before another r-class
command. See {help np_tdnn_diagnostics}, {help np_tdnn_plot},
{help weights_tdnn}, {help rdistnn_results}, and {help rdistnn_predict}.

{title:Example}

{phang2}{cmd:. np_tdnn y, xcols(x1 x2) x0(0 0) s1(3) s2(6)}
{phang2}{cmd:. np_tdnn_diagnostics, detail}
{phang2}{cmd:. np_tdnn_plot, type(weights) replace}
