{smcl}
{* *! version 0.1.2  2026-07-12}{...}
{viewerjumpto "Syntax" "np_dnn##syntax"}{...}
{viewerjumpto "Description" "np_dnn##description"}{...}
{viewerjumpto "Options" "np_dnn##options"}{...}
{viewerjumpto "Stored results" "np_dnn##results"}{...}
{viewerjumpto "Examples" "np_dnn##examples"}{...}

{title:Title}

{phang}
{bf:np_dnn} {hline 2} Nonparametric regression via distributional nearest neighbors

{pstd}Stata 17 or later is required.{p_end}

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
{synopt:{opt distance(string)}}{cmd:l2} (default) or {cmd:mahalanobis}{p_end}
{synopt:{opt method(string)}}{cmd:exact} (default) or {cmd:asymp}{p_end}

{syntab:Reporting}
{synopt:{opt level(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt nose}}skip rank-shift jackknife inference; permits s=n{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:np_dnn} estimates E[Y | X=x0] at one or more target covariate profiles.
It applies deterministic DNN rank weights to the outcome and returns a
localized weighted mean for every target.  Supply exactly one of {opt x0()}
or {opt x0matrix()}.  Rows of the target matrix are targets and columns follow
the order in {opt xcols()}.

{pstd}
The canonical covariate option {opt xcols()} corresponds to R's {cmd:x_cols}.
The earlier Stata option {opt x()} remains a compatibility alias but cannot be
combined with {opt xcols()}.  R's {cmd:X0}, {cmd:tuning_mode}, and
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

{marker options}{...}
{title:Options}

{phang}{opth xcols(varlist)} specifies the covariates defining distance and
rank localization.  Alias: {opt x()}.

{phang}{opt x0(numlist)} supplies one target with one value per covariate.

{phang}{opt x0matrix(matrix_name)} names an existing numeric matrix with at
least one row, exactly kx columns, and no missing entries.  It is mutually
exclusive with {opt x0()}.

{phang}{opt k(#)} sets how many nearest observations are selected within each
size-s subsample; it is not the full-sample positive-weight support size.  The
default is 1, matching the paper theory and R package.

{phang}{opt s(#)} sets the subset/kernel order, a localization or bandwidth
analogue whose direction is opposite to the usual bandwidth convention.  For
exact top-k weights, positive support is ranks 1 through n - s + k and has
size n - s + k; for the paper-default {cmd:k(1)}, this becomes n - s + 1.
Holding n and k fixed, larger s shrinks support and tightens localization,
while smaller s spreads positive weight across more observations and raises
the effective local sample size.  The order requires k <= s.  With jackknife
inference, s must be less than n; with {opt nose}, s may equal n.  Omitting
{opt s()} uses floor(0.4*n^((kx+1)/(kx+2))).

{phang}{opt tuningmode(grid)} evaluates a tuning grid after computing target
ranks once.  The primary {cmd:e(b)}, {cmd:e(V)}, {cmd:e(s)}, and ESS results
remain anchored to the marked default candidate.  Grid mode cannot be combined
with explicit {opt s()}.

{phang}{opt lambda(#)} multiplies the paper order rule.  The finite-sample
default is 0.4.

{phang}{opt lambdagrid(numlist)} supplies positive multipliers.  The default
sensitivity grid is .3, .4, .5, and .75.  Rounded duplicate orders remain separate
candidate rows but are computed once.

{phang}{opt sgrid(numlist)} directly supplies positive integer orders.
{opt sdefault(#)} selects the default order; otherwise the first entry is used.

{phang}{opt distance(string)} selects Euclidean {cmd:l2} distance or
{cmd:mahalanobis} distance.  The full-sample covariance metric is held fixed
through jackknife deletions.

{phang}{opt method(string)} selects exact finite-sample DNN weights or the
asymptotic beta-order-statistic approximation.

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
{synopt:{cmd:e(V)}}joint rank-shift jackknife covariance; omitted with {cmd:nose}{p_end}
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

{p2col 5 25 29 2:Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:np_dnn}{p_end}
{synopt:{cmd:e(xvars)}, {cmd:e(x_cols)}}covariate variable names{p_end}
{synopt:{cmd:e(target_source)}}{cmd:x0} or R-aligned {cmd:X0}{p_end}
{synopt:{cmd:e(distance)}}distance metric{p_end}
{synopt:{cmd:e(method)}}weight method{p_end}
{synopt:{cmd:e(tuning_mode)}}{cmd:default} or {cmd:grid}{p_end}
{synopt:{cmd:e(tuning_source)}}paper rule, explicit s, lambda grid, or s grid{p_end}
{synopt:{cmd:e(se_method)}, {cmd:e(jackknife)}}{cmd:rank_shift}/{cmd:fast} or {cmd:none}{p_end}

{pstd}
{cmd:distance(mahalanobis)} requires at least two complete estimation
observations, including with {cmd:nose}. Rank-deficient covariance matrices for
larger samples continue to use the generalized inverse.

{pstd}
Columns of {cmd:e(diagnostics)} are ESS, maximum weight, top-1/top-5/top-10
weight shares, nonzero support count, distance availability, rank-k and rank-s
radii, weighted mean and median distance, outcome availability, weighted
outcome mean and standard deviation, and weighted outcome 10th, 50th, and
90th percentiles.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. np_dnn y, xcols(x1 x2) x0(0 0)}{p_end}
{phang2}{cmd:. matrix X0 = (0,0 \ .5,-.25 \ -.5,.75)}{p_end}
{phang2}{cmd:. np_dnn y, xcols(x1 x2) x0matrix(X0) k(2) s(20)}{p_end}
{phang2}{cmd:. matrix list e(V)}{p_end}
{phang2}{cmd:. np_dnn y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)}{p_end}
{phang2}{cmd:. matrix list e(tuning_estimates)}{p_end}
{phang2}{cmd:. np_dnn_diagnostics}{p_end}
{phang2}{cmd:. np_dnn_plot, type(dashboard)}{p_end}
