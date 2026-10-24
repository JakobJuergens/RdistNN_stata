{smcl}
{* *! version 0.1.3  2026-10-17}{...}
{viewerjumpto "Syntax" "np_dnn##syntax"}{...}
{viewerjumpto "Description" "np_dnn##description"}{...}
{viewerjumpto "Sample and targets" "np_dnn##sample"}{...}
{viewerjumpto "Inference" "np_dnn##inference"}{...}
{viewerjumpto "Options" "np_dnn##options"}{...}
{viewerjumpto "Stored results" "np_dnn##results"}{...}
{viewerjumpto "Examples" "np_dnn##examples"}{...}
{viewerjumpto "Postestimation" "np_dnn##postestimation"}{...}

{title:Title}

{phang}
{bf:np_dnn} {hline 2} Nonparametric regression via distributional nearest neighbors

{pstd}
Stata 17 or later is required. See {help rdistnn} for the package version
and overview.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:np_dnn} {it:outcome} {ifin}{cmd:,}
{opth xcols(varlist)}
{c -(}{opt x0(numlist)} | {opt x0matrix(matrix_name)}{c )-}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Localization (required)}
{synopt:{opth xcols(varlist)}}covariates used to define the neighborhood{p_end}
{synopt:{opt x0(numlist)}}single evaluation point; one value per covariate{p_end}
{synopt:{opt x0matrix(matrix_name)}}named m x kx target matrix corresponding to R's {cmd:X0}{p_end}

{syntab:Tuning}
{synopt:{opt k(#)}}nearest observations selected per size-s subsample; default {cmd:k(1)}{p_end}
{synopt:{opt s(#)}}kernel order; default uses floor(0.4*n^((kx+1)/(kx+2))){p_end}
{synopt:{opt tuningmode(string)}}{cmd:default} (default) or integrated {cmd:grid}{p_end}
{synopt:{opt lambda(#)}}default paper-rule multiplier; default 0.4{p_end}
{synopt:{opt lambdagrid(numlist)}}candidate multipliers for grid mode{p_end}
{synopt:{opt sgrid(numlist)}}direct integer candidate orders for grid mode{p_end}
{synopt:{opt sdefault(#)}}default order for {cmd:sgrid()}{p_end}
{synopt:{opt distance(string)}}l1, l2 (default), l1_scaled, l2_scaled, or mahalanobis{p_end}
{synopt:{opt scalerule(string)}}sd (default), mad, iqr, range, or custom{p_end}
{synopt:{opt scalevec(numlist)}}custom scales, one per covariate{p_end}
{synopt:{opt scaleeps(#)}}positive scale floor; default 1e-12{p_end}
{synopt:{opt method(string)}}{cmd:exact} (default) or {cmd:asymp}{p_end}

{syntab:Reporting}
{synopt:{opt level(#)}}pointwise confidence level; Stata default, normally 95{p_end}
{synopt:{opt nose}}skip rank-shift jackknife inference; permits s=n{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
The documented sampling framework uses independent observations. Cluster and
survey-design inference are not supplied. At each target, observations with
exactly equal represented selected-metric distance share the mean of their
assigned rank weights, including each deletion's surviving group. No tolerance
groups nearby distances. Aligned row permutations preserve estimates and
covariance up to numerical precision. The paper assumes a continuous distance
distribution, under which old and new rules agree almost surely; inference with
persistent distance ties is not established. Reported order fractions
are descriptive and do not establish inference conditions.

{pstd}
{cmd:np_dnn} estimates E[Y | X=x0] at one or more target covariate profiles.
It applies deterministic DNN rank weights to the outcome and returns a
localized weighted mean for every target.  Supply exactly one of {opt x0()}
or {opt x0matrix()}.  Rows of the target matrix are targets and columns follow
the order in {opt xcols()}.

{pstd}
The canonical covariate option {opt xcols()} corresponds to R's {cmd:x_cols}.
R's {cmd:X0}, {cmd:tuning_mode}, and
{cmd:lambda_grid} correspond to {opt x0matrix()}, {opt tuningmode()}, and
{opt lambdagrid()}.

{pstd}
Standard errors use the fast rank-shift delete-1 jackknife.  For multiple
targets, {cmd:e(V)} is the full joint covariance matrix, retaining cross-target
covariance.  {opt nose} omits {cmd:e(V)} and permits an order equal to the
estimation-sample size.  Initial display and replay explicitly report that
inference was not requested; {cmd:e(se_method)} and {cmd:e(jackknife)} are
both {cmd:none}.

{pstd}
Every fit stores R-aligned per-target weight, distance, and outcome diagnostics.
Use {help np_dnn_diagnostics} to display or save them and {help np_dnn_plot}
for estimate, tuning, localization, outcome, weight, and dashboard plots.

{pstd}
The outcome may also be a caller-defined scalar signal, with the resulting
interpretation determined by that signal. To estimate treatment effects from
supplied nuisances or CATE pseudo-outcomes, see {help cate_dnn}. For two-scale
nonparametric regression, see {help np_tdnn}.

{marker sample}{...}
{title:Sample and targets}

{pstd}
{cmd:if} and {cmd:in} restrict candidate training observations, not target
rows. Observations missing the outcome or any {opt xcols()} variable are
excluded. {cmd:e(N)} records the resulting complete-case count;
{cmd:e(sample)} identifies the training rows while the matching data are
active. The default order and covariance metric use this complete sample.
An empty sample raises an error.

{pstd}
The target matrix must have at least one row, exactly one column per
localization covariate, and no missing entries. Its columns follow
{opt xcols()} order; column names do not reorder target values. Resolved
targets are retained in {cmd:e(targets)} with rows {cmd:target1},
{cmd:target2}, etc. Targets need not be observed training rows. Check
distance and concentration diagnostics when assessing whether a profile
has suitable nearby observations.

{marker inference}{...}
{title:Inference}

{pstd}
The point estimate is sum_i w_i(x0)*Y_i. With inference requested, the
delete-1 rank-shift jackknife holds the supplied outcomes fixed. L2 deletion
ranks coincide with recomputing distances on the remaining sample. For
Mahalanobis distance, the full-sample covariance metric is held fixed;
the jackknife does not re-estimate that metric under deletion.

{pstd}
Displayed confidence intervals are normal pointwise intervals for local
conditional means. Valid inference needs the regression assumptions and
localization-bias control; the jackknife does not remove smoothing bias.
Batch {cmd:e(V)} retains cross-target covariance and can be used for
prespecified linear contrasts. Pointwise intervals over many targets do not
provide simultaneous confidence bands.

{pstd}
NP diagnostics describe the weighted outcome distribution as well as the
weights. {cmd:weighted_sd_y} is outcome spread, not a standard error of the
mean. Weighted outcome q10/q50/q90 describe outcome distribution rather than
uncertainty for the mean. Use the estimate plot for mean intervals
and the outcome plot for these distributional summaries.

{marker options}{...}
{title:Options}

{phang}{opth xcols(varlist)} specifies the covariates defining distance and
rank localization.

{phang}{opt x0(numlist)} supplies one target with one value per covariate.

{phang}{opt x0matrix(matrix_name)} names an existing numeric matrix with at
least one row, exactly kx columns, and no missing entries.  It is mutually
exclusive with {opt x0()}.

{phang}{opt k(#)} sets how many nearest observations are selected within each
size-s subsample; it is not the full-sample positive-weight support size.  The
default is 1, matching the paper theory and R package.
Larger fixed k is an implementation-supported extension; the stated paper
theory uses k=1.

{phang}{opt s(#)} sets the subset/kernel order, a localization or bandwidth
analogue whose direction is opposite to the usual bandwidth convention.  For
exact top-k weights, positive support is ranks 1 through n - s + k and has
size n - s + k; for the paper-default {cmd:k(1)}, this becomes n - s + 1.
Holding n and k fixed, larger s shrinks support and tightens localization,
while smaller s spreads positive weight across more observations and raises
the effective local sample size.  The order requires k <= s.  With jackknife
inference, s must be less than n; with {opt nose}, s may equal n.  Omitting
{opt s()} uses floor(0.4*n^((kx+1)/(kx+2))).
{cmd:s(0)} is the sentinel for the default rule. The uniform cases s=1 or
k=s give no localization. An infeasible default order raises an error.

{pstd}
The support formula describes the unaveraged rank table. Averaging assigned
weights across exact-distance groups can extend nonzero observation weights to
the end of a group crossing that edge. Use actual returned weights for realized
support and concentration; rank-only tables do not contain target-specific groups.

{phang}{opt tuningmode(grid)} evaluates a tuning grid after computing target
ranks once.  The primary {cmd:e(b)}, {cmd:e(V)}, {cmd:e(s)}, and ESS results
remain anchored to the marked default candidate.  Grid mode cannot be combined
with explicit {opt s()}.
It is a sensitivity calculation, not automatic selection. The grid-specific
controls require {cmd:tuningmode(grid)}.

{phang}{opt lambda(#)} multiplies the paper order rule.  The finite-sample
default is 0.4.

{phang}{opt lambdagrid(numlist)} supplies positive multipliers.  The default
sensitivity grid is .3, .4, .5, and .75.  Rounded duplicate orders remain separate
candidate rows but are computed once.

{phang}{opt sgrid(numlist)} directly supplies positive integer orders.
{opt sdefault(#)} selects the default order; otherwise the first entry is used.
A default absent from {opt sgrid()} is prepended. {opt lambdagrid()} and
{opt sgrid()} are mutually exclusive; {opt sdefault()} requires {opt sgrid()}.
Infeasible candidates remain in the table with missing estimates. The default
candidate must be feasible.

{phang}{opt distance(string)} selects {cmd:l1}, Euclidean {cmd:l2} (default),
{cmd:l1_scaled}, {cmd:l2_scaled}, or {cmd:mahalanobis} distance.
The full-sample covariance metric is held fixed
through jackknife deletions. L2 distances depend on covariate scales.
Mahalanobis requires at least two complete observations, including with
{opt nose}; one localization covariate is permitted. Rank-deficient covariance
matrices use centered-design SVD whitening, retaining singular values greater
than max(n,p)*2^-52 times the largest singular value. Directions in its null
space do not contribute to distance.

{phang}{opt method(string)} selects exact finite-sample DNN weights or the
asymptotic beta-order-statistic approximation, written {cmd:method(asymp)}.
The exact support formula n-s+k does not apply to the approximation, which
is defined over the full rank grid. Weights are normalized; there is no
Stata {cmd:normalize()} option.

{phang}{opt level(#)} sets the normal pointwise confidence level. It uses
Stata's configured level, normally 95. Replay accepts {opt level()} without
refitting or changing stored estimates and covariance.

{phang}{opt nose} skips the jackknife, omits {cmd:e(V)}, and permits s=n.
Estimates and diagnostics remain available; common SEs, tests and interval
bounds are missing. {cmd:e(se_method)} and {cmd:e(jackknife)} are {cmd:none}.

{title:Distance scaling}

{pstd}
{cmd:distance(l1)} uses Manhattan distance. {cmd:distance(l2)} uses Euclidean
distance (default). {cmd:l1_scaled} and {cmd:l2_scaled} divide each coordinate
by a scale fitted on the complete estimation sample before computing the
corresponding distance. {cmd:scalerule(sd)} is the default; alternatives are
{cmd:mad} (median absolute deviation, constant 1.4826), {cmd:iqr} (R type-7
quartiles), {cmd:range}, and {cmd:custom}. Sample SD uses the n-1 divisor and
requires two observations. Supply {cmd:scalevec(numlist)} only with
{cmd:scalerule(custom)}, with one finite scale per covariate. Custom values
are absolute scales, not inverse scales. Each scale is floored at the positive
finite {cmd:scaleeps()} (default 1e-12), including constant covariates and zero
custom scales. Scaling controls do not change unscaled or Mahalanobis metrics.
The fitted scales or Mahalanobis whitening factor remain fixed across
jackknife deletions and saved-model predictions. Ranks use represented selected-metric distance then original row index.
At each target, observations with exactly equal distance share the mean of their
assigned rank weights, including different covariate profiles and zero-weight
ranks beyond nominal support. No tolerance or rounding groups nearby distances.
Each deletion removes the actual observation before averaging the n-1 table
over each remaining group using its survivor count. Signed and unnormalized
mass is preserved; individual weights, support and ESS can change. Aligned row
permutations preserve estimates and covariance up to numerical precision.
The paper assumes a continuous distance distribution, under which old and new
rules agree almost surely. Inference with persistent distance ties and a uniform
MSE improvement are not established; outcomes and nuisances do not define groups.

{marker results}{...}
{title:Stored results}

{pstd}
CATE, CLATE and NP regression share a reporting contract. Fit and replay
return a genuine 9-by-m {cmd:r(table)} with rows {cmd:b se z pvalue ll ul df crit eform},
and an m-by-6 {cmd:r(summary)} with columns
{cmd:estimate std_error z pvalue ci_lower ci_upper}. {cmd:r(level)} records
the confidence level. These are temporary r-class results; copy them before
running another r-class command. Replay with {cmd:level()} changes reported
bounds without changing the fit. Full target order is retained, including
undefined CLATE ratios. {cmd:nose} leaves uncertainty missing.

{pstd}
Common stored matrices are {cmd:e(estimates)} and {cmd:e(se_targets)}
(1-by-m) and {cmd:e(estimate_vcov)} (m-by-m); existing names remain available.
{cmd:e(signature_vars)} and {cmd:e(ordered_signature)} support value and
observation-order checks for observation-level plots. Stored-result plots
and summaries do not require the original data. CLATE retains its status
mask: zero score variance can have SE zero with no Wald statistic or interval.
For CATE/NP, zero variance gives point confidence bounds; Stata represents
an infinite z statistic as missing.

{pstd}{cmd:np_dnn} stores the following in {cmd:e()}:{p_end}

{synoptset 25}{...}
{p2col 5 25 29 2:Matrices}{p_end}
{synopt:{cmd:e(b)}}one conditional-mean estimate per target{p_end}
{synopt:{cmd:e(estimates)}}1 x m alias of {cmd:e(b)}{p_end}
{synopt:{cmd:e(se_targets)}}1 x m marginal standard errors; missing with {cmd:nose}{p_end}
{synopt:{cmd:e(V)}}joint rank-shift jackknife covariance; omitted with {cmd:nose}{p_end}
{synopt:{cmd:e(estimate_vcov)}}joint covariance; all missing with {cmd:nose}{p_end}
{synopt:{cmd:e(targets)}}resolved m x kx target matrix{p_end}
{synopt:{cmd:e(x0_matrix)}}alias of {cmd:e(targets)}{p_end}
{synopt:{cmd:e(ess_targets)}}one effective sample size per target{p_end}
{synopt:{cmd:e(neff_targets)}}alias of {cmd:e(ess_targets)}{p_end}
{synopt:{cmd:e(diagnostics)}}17 per-target weight, distance, and outcome fields{p_end}
{synopt:{cmd:e(target_diagnostics)}}alias of {cmd:e(diagnostics)}{p_end}
{synopt:{cmd:e(tuning_candidates)}}candidate metadata in grid mode{p_end}
{synopt:{cmd:e(tuning_estimates)}}candidate-by-target estimates in grid mode{p_end}
{synopt:{cmd:e(tuning_se)}}candidate-by-target standard errors in grid mode{p_end}
{synopt:{cmd:e(tuning_diagnostics)}}candidate ESS and concentration diagnostics{p_end}

{p2col 5 25 29 2:Scalars}{p_end}
{synopt:{cmd:e(N)}}estimation-sample size{p_end}
{synopt:{cmd:e(m)}}number of targets{p_end}
{synopt:{cmd:e(k)}}neighbor index{p_end}
{synopt:{cmd:e(s)}}default kernel order{p_end}
{synopt:{cmd:e(lambda)}}default multiplier, or missing for direct s grids{p_end}
{synopt:{cmd:e(kx)}, {cmd:e(p)}}covariate dimension{p_end}
{synopt:{cmd:e(ess)}, {cmd:e(neff)}}single-target ESS{p_end}
{synopt:{cmd:e(ess_min)}, {cmd:e(ess_max)}}ESS range across targets{p_end}
{synopt:{cmd:e(neff_min)}, {cmd:e(neff_max)}}aliases of the ESS range{p_end}
{synopt:{cmd:e(tuning_default_candidate)}}marked default candidate row, in grid mode{p_end}
{synopt:{cmd:e(tuning_n_candidates)}}number of candidate rows, in grid mode{p_end}
{synopt:{cmd:e(tuning_n_feasible)}, {cmd:e(tuning_n_skipped)}}candidate feasibility counts, in grid mode{p_end}
{synopt:{cmd:e(tuning_n_failed)}, {cmd:e(tuning_n_unique)}}failed and unique feasible computation counts, in grid mode{p_end}

{p2col 5 25 29 2:Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:np_dnn}{p_end}
{synopt:{cmd:e(depvar)}}positional outcome variable{p_end}
{synopt:{cmd:e(vce)}}{cmd:rank_shift} or {cmd:none} with {cmd:nose}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Jackknife}; empty when uncertainty is unavailable{p_end}
{synopt:{cmd:e(xvars)}, {cmd:e(x_cols)}}covariate variable names{p_end}
{synopt:{cmd:e(x0)}, {cmd:e(x0matrix)}}supplied single-target values or source target-matrix name{p_end}
{synopt:{cmd:e(target_source)}}{cmd:x0} or R-aligned {cmd:X0}{p_end}
{synopt:{cmd:e(distance)}}distance metric{p_end}
{synopt:{cmd:e(method)}}weight method{p_end}
{synopt:{cmd:e(tuning_mode)}}{cmd:default} or {cmd:grid}{p_end}
{synopt:{cmd:e(tuning_source)}}paper rule, explicit s, lambda grid, or s grid{p_end}
{synopt:{cmd:e(se_method)}, {cmd:e(jackknife)}}{cmd:rank_shift}/{cmd:fast} or {cmd:none}{p_end}
{synopt:{cmd:e(outcome)}, {cmd:e(signal)}}positional outcome variable{p_end}
{synopt:{cmd:e(tuning_status_codes)}, {cmd:e(tuning_reason_codes)}}candidate status/reason maps{p_end}
{synopt:{cmd:e(cmdline)}}estimation command{p_end}

{pstd}
{cmd:distance(mahalanobis)} requires at least two complete estimation
observations, including with {cmd:nose}. Rank-deficient covariance matrices for
larger samples use the SVD Moore-Penrose pseudoinverse.

{pstd}
Columns of {cmd:e(diagnostics)} are ESS, maximum weight, top-1/top-5/top-10
weight shares, nonzero support count, distance availability, rank-k and rank-s
radii, weighted mean and median distance, outcome availability, weighted
outcome mean and standard deviation, and weighted outcome 10th, 50th, and
90th percentiles.

{pstd}
Candidate metadata columns are {cmd:candidate_id computed_id is_default
feasible reason status lambda s}. Status/feasible is 0 for infeasible and
1 for available. Reason is 0 for OK, 1 for k>s, 2 for a jackknife order
s>=n, and 3 for s>n in an estimate-only fit. Rounded duplicates retain
their rows and share {cmd:computed_id}. Candidate covariance matrices are
not retained; {cmd:e(tuning_se)} contains marginal SEs.

{pstd}
Effective sample size is 1/sum_i w_i^2. It measures concentration, not
physical localization, support count or the assumptions for Gaussian
inference. {cmd:e(diagnostics)} contains the default candidate's target
diagnostics even after a grid fit.

{marker examples}{...}
{title:Examples}

{pstd}Generate a reproducible regression example; this clears the dataset.{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 180}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. gen double x1 = rnormal()}{p_end}
{phang2}{cmd:. gen double x2 = rnormal()}{p_end}
{phang2}{cmd:. gen double y = sin(x1) - .5*x2 + .3*rnormal()}{p_end}
{phang2}{cmd:. np_dnn y, xcols(x1 x2) x0(0 0)}{p_end}
{phang2}{cmd:. matrix one_target_summary = r(summary)}{p_end}

{pstd}Fit a target matrix with joint covariance.{p_end}
{phang2}{cmd:. matrix X0 = (0,0 \ .5,-.25 \ -.5,.75)}{p_end}
{phang2}{cmd:. np_dnn y, xcols(x1 x2) x0matrix(X0) k(2) s(20)}{p_end}
{phang2}{cmd:. matrix list e(V)}{p_end}

{pstd}Inspect the grid while retaining the ordinary default fit.{p_end}
{phang2}{cmd:. np_dnn y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)}{p_end}
{phang2}{cmd:. matrix list e(tuning_estimates)}{p_end}
{phang2}{cmd:. estimates store np_fit}{p_end}
{phang2}{cmd:. np_dnn, level(90)}{p_end}
{phang2}{cmd:. matrix summary90 = r(summary)}{p_end}
{phang2}{cmd:. np_dnn_diagnostics}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, allcandidates}{p_end}
{phang2}{cmd:. np_dnn_plot, type(estimates) replace}{p_end}
{phang2}{cmd:. np_dnn_plot, type(outcome) replace}{p_end}
{phang2}{cmd:. np_dnn_plot, type(dashboard)}{p_end}

{pstd}Request an estimate without inference.{p_end}
{phang2}{cmd:. np_dnn y if !missing(y,x1,x2), xcols(x1 x2) x0(0 0) nose}{p_end}
{phang2}{cmd:. estimates restore np_fit}{p_end}
{phang2}{cmd:. np_dnn}{p_end}

{marker postestimation}{...}
{title:Postestimation}

{pstd}
Each fit replaces active {cmd:e()} results. Use {cmd:estimates store} and
{cmd:estimates restore} to keep fits and reactivate the matching
postestimation commands. Single-target coefficients are named {cmd:mu};
batch coefficients are {cmd:mu1}, {cmd:mu2}, etc., in target order.
For a prespecified first-minus-second target contrast after a batch fit with
inference, use {cmd:lincom mu1-mu2}.

{pstd}
Copy {cmd:r(summary)} or {cmd:r(table)} immediately after fitting or replay:
other r-class commands replace them. {help np_dnn_diagnostics} selects
target or candidate tables and can save a {cmd:.dta} file. {help np_dnn_plot}
provides estimates, diagnostics, tuning, weights, outcome and dashboard
views. Current effective-kernel views use saved numeric observations and
geometry after active data changes.

{pstd}
{help rdistnn_predict} evaluates saved frozen signals at new targets with
optional covariance. {help rdistnn_candidate} extracts standalone grid fits.
{help rdistnn_results} exports numeric reports, including long weight tables,
without graphs. {help weights_dnn} exports one target's weights as a variable.

{pstd}
See also {help rdistnn}, {help cate_dnn}, {help weights_dnn},
{help np_dnn_diagnostics}, and {help np_dnn_plot}.

{title:Saved model reporting and reuse}
{pstd}
Current fits retain versioned numeric observations, frozen signals, ranks and
geometry in e(model_*), including in estimates save/use files.
{help rdistnn_results} exports numeric reports without a graph;
{help rdistnn_candidate} extracts reusable grid fits;
{help rdistnn_predict} evaluates saved signals at new targets with optional
frozen-score covariance. Active data and e() are preserved by default.
post explicitly activates a derived fit; store() retains a named one.
Fits without the current required numeric snapshot need refitting.
