{smcl}
{* *! version 0.0.8  2026-07-12}{...}
{viewerjumpto "Syntax"          "cate_dnn##syntax"}{...}
{viewerjumpto "Description"     "cate_dnn##description"}{...}
{viewerjumpto "Inference contract" "cate_dnn##inference_contract"}{...}
{viewerjumpto "Options"         "cate_dnn##options"}{...}
{viewerjumpto "Stored results"  "cate_dnn##results"}{...}
{viewerjumpto "Examples"        "cate_dnn##examples"}{...}
{viewerjumpto "Postestimation"  "cate_dnn##postestimation"}{...}
{viewerjumpto "References"      "cate_dnn##references"}{...}

{title:Title}

{phang}
{bf:cate_dnn} {hline 2} Conditional average treatment effect via
distributional nearest neighbors

{pstd}
Stata 17 or later is required.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:cate_dnn}
{it:outcome_or_signal}
{ifin}{cmd:,}
{opth xcols(varlist)}
{c -(}{opt x0(numlist)} | {opt x0matrix(matrix_name)}{c )-}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Localization (required)}
{synopt:{opth xcols(varlist)}}covariates used to define the neighborhood{p_end}
{synopt:{opt x0(numlist)}}single evaluation point; one value per covariate{p_end}
{synopt:{opt x0matrix(matrix_name)}}named m x kx matrix corresponding to R's {cmd:X0}{p_end}

{syntab:Score input (choose one path)}
{synopt:{opt pseudooutcome}}positional variable is R's {cmd:pseudo_outcome}{p_end}
{synopt:{opth tcol(varname)}}binary treatment, corresponding to R's {cmd:t_col}{p_end}
{synopt:{opth ecol(varname)}}propensity, corresponding to R's {cmd:e_col}{p_end}
{synopt:{opth m0col(varname)}}control outcome mean, corresponding to R's {cmd:m0_col}{p_end}
{synopt:{opth m1col(varname)}}treated outcome mean, corresponding to R's {cmd:m1_col}{p_end}
{synopt:{opth mcol(varname)}}pooled outcome mean, corresponding to R's {cmd:m_col}{p_end}
{synopt:{opt clip(#)}}propensity clipping constant in (0, 0.5); default {cmd:1e-8}{p_end}

{syntab:Tuning}
{synopt:{opt k(#)}}nearest observations selected per size-s subsample; default {cmd:k(1)}{p_end}
{synopt:{opt s(#)}}kernel order; default {cmd:s(0)} uses
    floor(0.4*n^((kx+1)/(kx+2))){p_end}
{synopt:{opt tuningmode(string)}}{cmd:default} (default) or integrated {cmd:grid}{p_end}
{synopt:{opt lambda(#)}}default multiplier in the paper tuning rule; default 0.4{p_end}
{synopt:{opt lambdagrid(numlist)}}candidate multipliers for {cmd:tuningmode(grid)}{p_end}
{synopt:{opt sgrid(numlist)}}direct integer candidate orders for grid mode{p_end}
{synopt:{opt sdefault(#)}}default order for {cmd:sgrid()}; otherwise its first entry{p_end}
{synopt:{opt distance(string)}}{cmd:l2} (default) or {cmd:mahalanobis}{p_end}
{synopt:{opt method(string)}}{cmd:exact} (default) or {cmd:asymp}{p_end}

{syntab:Reporting}
{synopt:{opt level(#)}}set confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt nose}}skip the rank-shift jackknife standard error; permits s=n{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:cate_dnn} estimates the conditional average treatment effect
tau(x0) = E[Y(1) - Y(0) | X = x0] at one or more covariate points using the
distributional nearest-neighbor (DNN) estimator of Juergens (2025).  Supply
exactly one of {opt x0()} for a single target or {opt x0matrix()} for a named Stata
matrix whose rows are targets and whose columns correspond to {opt xcols()}.

{pstd}
The estimator applies DNN rank weights w_i(x0) to per-observation
pseudo-outcomes, yielding a localized weighted mean.  Three score-input paths
mirror the R package: a precomputed pseudo-outcome, a full AIPW score, or a
Robinson/R-learner score.

{pstd}
With {opt pseudooutcome}, the positional variable is used directly as R's
{cmd:pseudo_outcome}; no outcome or nuisance variables are supplied.  The full
AIPW pseudo-outcome is

{phang2}phi_i = D_i*(Y_i-m1_i)/eclip_i - (1-D_i)*(Y_i-m0_i)/(1-eclip_i) + (m1_i-m0_i){p_end}

{pstd}
The Robinson pseudo-outcome is

{phang2}phi_i = (D_i-eclip_i)*(Y_i-mhat_i)/(eclip_i*(1-eclip_i)).{p_end}

{pstd}
For nuisance-based paths, eclip_i = min(max(e_i, clip), 1-clip), where
{opt clip()} defaults to 1e-8. Raw propensity values must be finite and lie in
[0,1]; values outside that range are rejected rather than mapped back into it.
Supply {opt ecol()} and either {opt m0col()} with {opt m1col()} (AIPW),
or {opt mcol()} (Robinson).

{pstd}
{bf:Inference warning:} Neyman orthogonality alone does not validate arbitrary
same-sample nuisance estimates.  See {it:Inference contract} below before
using full-sample predictions for inference.

{pstd}
Standard errors are computed via the rank-shift jackknife.  For
{cmd:distance(l2)}, this is the exact delete-1 recomputation because the
remaining distances do not change.  For {cmd:distance(mahalanobis)}, the
sample covariance metric is computed from the full estimation sample and
held fixed across jackknife deletions.  With multiple targets, {cmd:e(V)} is
the full joint covariance of the leave-one-out target-estimate vector.  Its
diagonal standard errors match separate single-target calls and the R package;
its off-diagonal entries retain cross-target covariance for Stata inference.

{marker inference_contract}{...}
{title:Inference contract}

{pstd}
Orthogonality controls first-order population nuisance drift.  It does not, by
itself, justify treating arbitrary nuisance predictions estimated on the same
observations as valid inputs for the reported inference.

{pstd}
Valid full-sample reuse additionally requires the nuisance-rate and
full-state stability conditions in the accompanying paper.  {cmd:cate_dnn}
computes the localized estimator and rank-shift variance procedure, but it
cannot diagnose whether an arbitrary learner satisfies those theoretical
conditions.

{pstd}
Users may supply oracle quantities, cross-fitted or out-of-fold predictions,
or full-sample predictions from learners for which the paper's rate and
stability conditions are credible.  The design is deliberately modular:
Stata does not automatically fit or select the nuisance models.

{marker options}{...}
{title:Options}

{dlgtab:Localization}

{phang}
{opth xcols(varlist)} specifies the covariates X used to define the local
neighborhood around x0.

{phang}
{opt x0(numlist)} specifies the evaluation point.  Provide one number per
variable in {opt xcols()}, separated by spaces.  Example: {cmd:x0(0 0)} for
the two-dimensional origin.

{phang}
{opt x0matrix(matrix_name)} specifies a named numeric Stata matrix with one
target per row and one column per variable in {opt xcols()}.  It corresponds
to the R argument {cmd:X0}.  The matrix must exist before
calling {cmd:cate_dnn}, have no missing entries, and have exactly kx columns.
It is mutually exclusive with {opt x0()}.  Resolved values are stored in
{cmd:e(targets)} with canonical covariate column names.

{dlgtab:Score input}

{phang}
{opt pseudooutcome} specifies that the positional variable is an already
constructed CATE pseudo-outcome and corresponds to R's {cmd:pseudo_outcome}.
Do not combine it with nuisance options.

{phang}
{opth tcol(varname)} specifies the binary treatment indicator D (0 or 1),
corresponding to R's {cmd:t_col}.

{phang}
{opth ecol(varname)} provides e(X) = P(D=1|X), corresponding to R's
{cmd:e_col}.  Values must lie in [0, 1].  Values at or very near the endpoints
are replaced by the clipped propensity used consistently throughout the AIPW
or Robinson score and propensity diagnostics.

{phang}
{opt clip(#)} specifies the finite propensity clipping constant. It must lie
strictly between 0 and 0.5 and defaults to {cmd:1e-8}. It has no effect on a
precomputed pseudo-outcome.

{phang}
{opth m0col(varname)} provides E[Y|X, D=0], corresponding to R's
{cmd:m0_col}.  Use it with {opt m1col()} for full AIPW.

{phang}
{opth m1col(varname)} provides E[Y|X, D=1], corresponding to R's
{cmd:m1_col}.  Use it with {opt m0col()} for full AIPW.

{phang}
{opth mcol(varname)} provides E[Y|X] for the Robinson/R-learner path and
corresponds to R's {cmd:m_col}.  It is mutually exclusive with {opt m0col()}
and {opt m1col()}.

{dlgtab:Compatibility aliases}

{pstd}
Stata option names cannot contain underscores.  The canonical mapping removes
only R's underscores: {cmd:x_cols -> xcols()}, {cmd:t_col -> tcol()},
{cmd:e_col -> ecol()}, {cmd:m0_col -> m0col()}, {cmd:m1_col -> m1col()},
{cmd:m_col -> mcol()}, and {cmd:pseudo_outcome -> pseudooutcome}.  R's
{cmd:X0} matrix is {cmd:x0matrix()}, {cmd:tuning_mode} is
{cmd:tuningmode()}, and {cmd:lambda_grid} is {cmd:lambdagrid()}; {cmd:x0()}
and {cmd:lambda()} are identical in both APIs.
The earlier Stata names {cmd:x()}, {cmd:treat()}, {cmd:propensity()},
{cmd:m0()}, {cmd:m1()}, and {cmd:mhat()} remain accepted but cannot be combined
with their canonical counterparts.

{dlgtab:Tuning}

{phang}
{opt k(#)} sets how many nearest observations are selected within each
size-s subsample.  It is not the number of observations receiving positive
weight in the full-sample DNN rank-weight vector.  Default is {cmd:k(1)},
matching the paper theory and the R package.

{phang}
{opt s(#)} sets the subset/kernel order, a localization or bandwidth analogue
whose direction is opposite to the usual bandwidth convention.  For exact
top-k DNN weights, positive support ends at global rank n - s + k, so its size
is n - s + k.  With the paper-default {cmd:k(1)}, the edge and support size are
n - s + 1.  Holding n and k fixed, larger s shrinks the positive-weight support
and produces tighter localization; smaller s spreads weight across more
observations and increases the effective local sample size.  Default
{cmd:s(0)} uses floor(0.4*n^((kx+1)/(kx+2))).

{phang}
{opt tuningmode(string)} selects {cmd:default}, which evaluates only the
resolved default order and preserves the ordinary fit, or {cmd:grid}, which
evaluates a candidate grid while reusing one target-rank calculation.  As in
R, explicit {opt s()} cannot be combined with grid mode.  The main
{cmd:e(b)}, {cmd:e(V)}, {cmd:e(s)}, and effective-sample-size results remain
anchored to the candidate marked as default.

{phang}
{opt lambda(#)} is the positive multiplier in
floor(lambda*n^((kx+1)/(kx+2))).  Its finite-sample default is 0.4.

{phang}
{opt lambdagrid(numlist)} provides paper-rule multipliers in grid mode.
The default sensitivity grid is {cmd:.3 .4 .5 .75}.  Candidate rows are preserved even
when rounding maps multiple multipliers to the same order; the common order
is computed once.

{phang}
{opt sgrid(numlist)} is the Stata-native alternative for supplying candidate
orders directly.  Values must be positive integers.  {opt sdefault(#)} marks
the default order; if omitted, the first entry is the default.  A requested
default absent from {opt sgrid()} is prepended, matching R's synthetic-default
behavior.

{phang}
{opt distance(string)} selects the distance metric.  {cmd:l2} (default)
uses Euclidean distance.  {cmd:mahalanobis} uses the sample covariance
matrix (via {cmd:invsym()}) to form an ellipsoidal neighborhood;
recommended when covariates are strongly correlated
({it:cf.} DGP 5 of Juergens 2025).  Requires at least 2 covariates and at
least two complete estimation observations. Rank-deficient covariance matrices
with two or more observations continue to use {cmd:invsym()}.

{phang}
{opt method(string)} selects {cmd:exact} finite-sample DNN weights (default)
or the {cmd:asymp} approximation.

{phang}
{opt nose} skips the delete-1 rank-shift jackknife and posts only {cmd:e(b)}.
Without {opt nose}, {opt s()} must be strictly less than the estimation-sample
size.  With {opt nose}, {opt s()} may equal that size.

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

{pstd}
{cmd:cate_dnn} stores the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(m)}}number of evaluation targets{p_end}
{synopt:{cmd:e(neff)}}effective sample size for a one-target fit{p_end}
{synopt:{cmd:e(ess)}}R-aligned alias of {cmd:e(neff)} for a one-target fit{p_end}
{synopt:{cmd:e(neff_min)}}minimum target-specific effective sample size{p_end}
{synopt:{cmd:e(neff_max)}}maximum target-specific effective sample size{p_end}
{synopt:{cmd:e(ess_min)}, {cmd:e(ess_max)}}R-aligned aliases of the effective-sample-size range{p_end}
{synopt:{cmd:e(k)}}neighbor count k{p_end}
{synopt:{cmd:e(s)}}kernel order s{p_end}
{synopt:{cmd:e(lambda)}}default paper-rule multiplier; missing for {cmd:sgrid()}{p_end}
{synopt:{cmd:e(kx)}}covariate dimension{p_end}
{synopt:{cmd:e(p)}}R-aligned alias of {cmd:e(kx)}{p_end}
{synopt:{cmd:e(clip)}}propensity clipping constant used for nuisance-based scores{p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}1 x m vector of target-specific CATE estimates{p_end}
{synopt:{cmd:e(V)}}m x m joint jackknife covariance; not posted with {cmd:nose}{p_end}
{synopt:{cmd:e(targets)}}m x kx resolved target matrix{p_end}
{synopt:{cmd:e(x0_matrix)}}R-aligned alias of {cmd:e(targets)}{p_end}
{synopt:{cmd:e(neff_targets)}}1 x m target-specific effective sample sizes{p_end}
{synopt:{cmd:e(ess_targets)}}R-aligned alias of {cmd:e(neff_targets)}{p_end}
{synopt:{cmd:e(diagnostics)}}per-target weight, distance, and CATE-overlap diagnostics{p_end}
{synopt:{cmd:e(target_diagnostics)}}explicit alias of {cmd:e(diagnostics)}{p_end}
{synopt:{cmd:e(tuning_candidates)}}candidate schema in grid mode: ID, computed ID, default/feasible flags, reason/status codes, lambda, and s{p_end}
{synopt:{cmd:e(tuning_estimates)}}candidate-by-target estimate matrix in grid mode{p_end}
{synopt:{cmd:e(tuning_se)}}candidate-by-target marginal-SE matrix in grid mode{p_end}
{synopt:{cmd:e(tuning_diagnostics)}}candidate diagnostics: ESS, maximum weight, top-n, top-five share, and degeneracy flag{p_end}

{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cate_dnn}{p_end}
{synopt:{cmd:e(x0)}}single evaluation point supplied through {cmd:x0()}{p_end}
{synopt:{cmd:e(x0matrix)}}source matrix name supplied through {cmd:x0matrix()}{p_end}
{synopt:{cmd:e(target_source)}}{cmd:x0} or {cmd:X0}{p_end}
{synopt:{cmd:e(xvars)}}covariate variable names{p_end}
{synopt:{cmd:e(x_cols)}}R-aligned alias of {cmd:e(xvars)}{p_end}
{synopt:{cmd:e(distance)}}distance metric used{p_end}
{synopt:{cmd:e(method)}}weight method used{p_end}
{synopt:{cmd:e(score_method)}}{cmd:aipw}, {cmd:robinson}, or {cmd:pseudo_outcome}{p_end}
{synopt:{cmd:e(propensity_version)}}{cmd:clipped} for nuisance scores or {cmd:not_used}{p_end}
{synopt:{cmd:e(se_method)}}{cmd:rank_shift} or {cmd:none}{p_end}
{synopt:{cmd:e(jackknife)}}R-aligned {cmd:fast} or {cmd:none}{p_end}
{synopt:{cmd:e(tuning_mode)}}{cmd:default} or {cmd:grid}{p_end}
{synopt:{cmd:e(tuning_source)}}{cmd:paper}, {cmd:lambda_grid}, {cmd:sgrid}, or explicit {cmd:s}{p_end}
{synopt:{cmd:e(tuning_status_codes)}}numeric status-code map used by {cmd:e(tuning_candidates)}{p_end}
{synopt:{cmd:e(tuning_reason_codes)}}numeric reason-code map used by {cmd:e(tuning_candidates)}{p_end}
{synopt:{cmd:e(signal)}}positional outcome or pseudo-outcome variable{p_end}
{synopt:{cmd:e(t_col)}, {cmd:e(e_col)}}resolved treatment and propensity variables{p_end}
{synopt:{cmd:e(m0_col)}, {cmd:e(m1_col)}, {cmd:e(m_col)}}resolved nuisance variables{p_end}

{marker examples}{...}
{title:Examples}

{pstd}
{bf:Inference note:} The examples below use oracle quantities.  Orthogonality
alone is not sufficient for arbitrary same-sample nuisance estimates; such
reuse also requires the paper's nuisance-rate and full-state stability
conditions.{p_end}

{pstd}Setup: simulate data with oracle nuisances.{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 1000}{p_end}
{phang2}{cmd:. gen x1 = rnormal()}{p_end}
{phang2}{cmd:. gen x2 = rnormal()}{p_end}
{phang2}{cmd:. gen e_oracle = invlogit(0.5*x1 - 0.5*x2)}{p_end}
{phang2}{cmd:. gen D = (runiform() < e_oracle)}{p_end}
{phang2}{cmd:. gen tau_true = 0.5 + 0.5*sin(x1)}{p_end}
{phang2}{cmd:. gen m0_true = 1 + x1 - x2}{p_end}
{phang2}{cmd:. gen m1_true = m0_true + tau_true}{p_end}
{phang2}{cmd:. gen Y = m0_true + D*(m1_true-m0_true) + rnormal()*0.5}{p_end}

{pstd}Estimate CATE at the origin using oracle nuisances.{p_end}
{phang2}{cmd:. cate_dnn Y, tcol(D) xcols(x1 x2) x0(0 0)}{...}
{phang2}{cmd:    ecol(e_oracle) m0col(m0_true) m1col(m1_true)}{p_end}

{pstd}With Mahalanobis distance (for correlated covariates).{p_end}
{phang2}{cmd:. cate_dnn Y, tcol(D) xcols(x1 x2) x0(0.75 0)}{...}
{phang2}{cmd:    ecol(e_oracle) m0col(m0_true) m1col(m1_true)}{...}
{phang2}{cmd:    distance(mahalanobis)}{p_end}

{pstd}Use a precomputed pseudo-outcome and skip standard errors.{p_end}
{phang2}{cmd:. gen rho = D*(Y-m1_true)/e_oracle - (1-D)*(Y-m0_true)/(1-e_oracle) + m1_true-m0_true}{p_end}
{phang2}{cmd:. cate_dnn rho, xcols(x1 x2) x0(0 0) pseudooutcome nose}{p_end}

{pstd}Use the Robinson score with a pooled outcome regression.{p_end}
{phang2}{cmd:. gen mhat_true = e_oracle*m1_true + (1-e_oracle)*m0_true}{p_end}
{phang2}{cmd:. cate_dnn Y, tcol(D) xcols(x1 x2) x0(0 0)}{...}
{phang2}{cmd:    ecol(e_oracle) mcol(mhat_true)}{p_end}

{pstd}Evaluate the AIPW score at three targets in one call.{p_end}
{phang2}{cmd:. matrix X0 = (0,0 \ 0.5,-0.25 \ -0.5,0.75)}{p_end}
{phang2}{cmd:. cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0)}{...}
{phang2}{cmd:    ecol(e_oracle) m0col(m0_true) m1col(m1_true)}{p_end}
{phang2}{cmd:. matrix list e(targets)}{p_end}
{phang2}{cmd:. matrix list e(V)}{p_end}

{pstd}Evaluate the R paper tuning grid and inspect candidate results.{p_end}
{phang2}{cmd:. cate_dnn rho, xcols(x1 x2) x0matrix(X0) pseudooutcome tuningmode(grid)}{p_end}
{phang2}{cmd:. matrix list e(tuning_candidates)}{p_end}
{phang2}{cmd:. matrix list e(tuning_estimates)}{p_end}
{phang2}{cmd:. matrix list e(tuning_diagnostics)}{p_end}

{pstd}Inspect R-aligned diagnostics for the default fit.{p_end}
{phang2}{cmd:. matrix list e(diagnostics)}{p_end}

{pstd}
{cmd:e(diagnostics)} has one row per target.  Its columns are ESS, maximum
weight, fixed top-1/top-5/top-10 weight shares, nonzero-weight count,
distance-availability flag, rank-k and rank-s radii, weighted mean and median
distance, overlap-availability flag, treated/control mass and treated share,
propensity-availability flag, weighted mean propensity, and weighted
propensity 10th/90th quantiles.  Paper diagnostics then report the ordinary
Euclidean physical radius containing 90 percent of weight mass, an arm-ESS
availability flag, treated/control ESS computed after normalizing weights
within each arm, a score-contribution availability flag, and the share of
total absolute {it:w_i(score_i-tau)} accounted for by the five largest
contributions.  These are distinct from the rank-s radius, overall ESS, and
top-five raw-weight share.  Arm-specific fields are missing for the
precomputed-score path; physical-radius and score-contribution diagnostics
remain available.

{marker postestimation}{...}
{title:Postestimation}

{pstd}
Use {help cate_dnn_diagnostics} to inspect or save target and tuning-candidate
diagnostic tables.  Use {help cate_dnn_plot} for Stata-native analogues of the
R package's CATE estimate, diagnostics, tuning, overlap, weight, influence,
and dashboard plots.

{marker references}{...}
{title:References}

{phang}
Juergens, J. (2026).  Companion conditional Z-estimation working paper;
see the package README for current title and repository metadata.
