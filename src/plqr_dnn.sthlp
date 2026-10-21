{smcl}
{* *! version 0.1.1  17oct2026}{...}
{viewerjumpto "Syntax" "plqr_dnn##syntax"}{...}
{viewerjumpto "Description and score" "plqr_dnn##description"}{...}
{viewerjumpto "Options" "plqr_dnn##options"}{...}
{viewerjumpto "Moment covariance" "plqr_dnn##covariance"}{...}
{viewerjumpto "Stored results" "plqr_dnn##results"}{...}
{viewerjumpto "Examples" "plqr_dnn##examples"}{...}
{viewerjumpto "R companion and scope" "plqr_dnn##scope"}{...}

{title:Title}

{phang}
{cmd:plqr_dnn} {hline 2} Experimental localized PLQR moment evaluation at supplied coefficients

{pstd}
Stata 17 or later is required. See {help rdistnn} for package navigation.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:plqr_dnn} {it:outcome} {ifin}{cmd:,}
{opth x(varlist)} {opt x0(numlist)} {opth d(varlist)}
{opth g(varname)} {opth md(varlist)} {opt thetaeval(numlist)}
[{it:options}]

{p 8 17 2}{cmd:plqr_dnn} [{cmd:, level(}{it:#}{cmd:)}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required inputs}
{synopt:{opth x(varlist)}}numeric reporting covariates defining localization{p_end}
{synopt:{opt x0(numlist)}}one target, coordinates in x() order{p_end}
{synopt:{opth d(varlist)}}numeric partially linear regressors D{p_end}
{synopt:{opth g(varname)}}supplied additive shape nuisance g_0(U){p_end}
{synopt:{opth md(varlist)}}supplied E[D|U] variables, in d() order{p_end}
{synopt:{opt thetaeval(numlist)}}supplied coefficient vector, in d() order{p_end}
{syntab:Localization}
{synopt:{opt tau(#)}}quantile index in (0,1); default .5{p_end}
{synopt:{opt k(#)}}neighbors averaged within each size-s subsample; default 1{p_end}
{synopt:{opt s(#)}}integer subsample order, k<=s<n; default paper rule{p_end}
{synopt:{opt distance(string)}}l1, l2 (default), l1_scaled, l2_scaled, or mahalanobis{p_end}
{synopt:{opt scalerule(string)}}sd (default), mad, iqr, range, or custom{p_end}
{synopt:{opt scalevec(numlist)}}custom scales, one per covariate{p_end}
{synopt:{opt scaleeps(#)}}positive scale floor; default 1e-12{p_end}
{synopt:{opt method(string)}}exact (default) or asymp{p_end}
{synopt:{opt level(#)}}accepted on fit/replay; produces no confidence intervals or tests{p_end}
{synoptline}

{marker description}{...}
{title:Description and score}

{pstd}
{bf:Experimental status.} {cmd:plqr_dnn} evaluates a localized vector score
moment at the coefficients supplied in {cmd:thetaeval()}. It is not a coefficient estimator.
It does not solve
for a coefficient, estimate a conditional density, or report coefficient
standard errors, confidence intervals, p-values or significance tests.
Its API may change. The
{browse "https://jakobjuergens.com/research/#orthogonal-inference-conditional-z":paper}
gives statistical background.

{pstd}
For observation i and score component j, the score is

{phang2}
psi_ij(theta) = (tau - indicator[Y_i <= D_i'theta + g_i]) * (D_ij - mD_ij).

{pstd}
The indicator includes equality: an observation exactly at the threshold
receives indicator one. The command averages these scores using common DNN
weights at {cmd:x0()}. Its covariance is the covariance of that localized
vector moment, including cross-component entries. No local root solving is
performed and the supplied coefficient vector remains unchanged.

{pstd}
The shape nuisance {cmd:g()} is the additive g_0(U) in the threshold
D'theta+g_0(U); g_0 alone is not a complete conditional quantile.
{cmd:md()} supplies m_D(U)=E[D|U] to residualize each regressor.
U is the full adjustment state; these predictions can depend on variables
outside the reporting covariates {cmd:x()}. The package accepts evaluated
predictions and fits no nuisance models.

{pstd}
Y and g use outcome units. Each supplied coefficient has outcome units per
unit of its regressor. Because the indicator residual is dimensionless,
score/moment component j has regressor-j units; covariance entry (j,h) has
the product of regressor-j and regressor-h units. It is not coefficient
covariance. A small sample moment does not turn supplied theta into an
estimated coefficient, and oracle examples do not establish coverage for
arbitrary fitted nuisances.

{marker options}{...}
{title:Options and sample}

{phang}
{cmd:x()} selects the reporting covariates. {cmd:x0()} supplies exactly one
coordinate per variable, in that order. Initial estimation uses a single target;
this command has neither {cmd:x0matrix()} nor {cmd:xcols()}.
After fitting, {help rdistnn_predict} accepts a batch of new targets and returns
joint moment covariance while retaining supplied coefficients separately.

{phang}
{cmd:d()} lists q partially linear regressors, with q>=1. They need not be
binary. {cmd:md()} must have q prediction variables in exactly the same
order. {cmd:thetaeval()} must have q coefficients in that order.
Changing the order of one without the others changes the score.

{phang}
{cmd:g()} supplies the additive shape variable. {cmd:thetaeval()} supplies
the coefficient vector at which the moment is evaluated.

{phang}
The estimation sample consists of complete outcome, x, d, g and md rows in
the requested if/in range. Prepare aligned external nuisance predictions
before the call. {cmd:e(N)} and {cmd:e(sample)} identify the retained sample.

{phang}
{cmd:tau()} sets the quantile index, strictly between zero and one.
The default is .5, the median score.

{phang}
{cmd:k()} is the number of observations averaged within a size-s subsample,
not the full-sample positive-weight support. {cmd:s()} is the integer
subsample order. With complete-sample size n and p reporting covariates,
the default is floor(.4*n^((p+1)/(p+2))). The covariance calculation requires
1<=k<=s<n. There is no PLQR {cmd:nose} or tuning-grid option.

{phang}
For exact top-k DNN, positive support ends at global rank n - s + k, or n - s + 1
under k(1). Larger s tightens localization,
opposite to the usual bandwidth convention;
exact k=s gives uniform weights and no target-specific localization.
{cmd:method(asymp)} uses normalized Beta-rank weights. There is no separate
normalization option in this command.

{phang}
{cmd:distance(l2)} uses Euclidean geometry. {cmd:distance(mahalanobis)}
requires at least two reporting covariates and uses full-estimation-sample
covariance geometry, with generalized-inverse behavior for singular covariance.
Input observation position breaks distance ties. The metric stays fixed
during deletions. See {help weights_dnn} for the shared rank-weight rules.

{marker covariance}{...}
{title:Moment covariance and inferential scale}

{pstd}
The documented sampling framework uses independent observations; cluster
and survey-design inference are not supplied. Input row order breaks distance
ties, so permuting tied observations can change estimates.


{pstd}
Each deletion removes one observation, shifts the retained ranks, and
recomputes weights at the same target, s, k and method. Supplied coefficients,
shape predictions, residualizer predictions and observation-level scores
remain frozen. Full-sample Mahalanobis covariance is not re-estimated after
deletion. If L_i is the row-vector moment after deleting i, the covariance is
(n-1)/n times sum_i (L_i-mean(L))'*(L_i-mean(L)).

{pstd}
{cmd:e(moment_vcov)} already contains the estimator-scale covariance of the
weighted moment. Do not divide it by {cmd:e(neff)} again.
{cmd:e(neff)}=1/sum_i w_i^2 is an effective weight-concentration measure;
it is not a certificate of a CLT, nuisance accuracy or point-target bias control.

{pstd}
The square root of a diagonal entry describes a moment component's jackknife
SD. It is not a coefficient SE. Coefficient inference needs a valid Jacobian
and sandwich construction, together with its identification and asymptotic
conditions. This Stata command supplies neither a density input nor a
coefficient Jacobian/sandwich. Generic coefficient postestimation that treats
{cmd:e(b)} and {cmd:e(V)} as coefficient estimates/covariance uses the wrong
statistical object. The custom display and replay print moments/covariance
without coefficient tests; {cmd:level()} changes no interval because none is
computed. There are no Stata PLQR diagnostics or plot commands.

{marker results}{...}

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
jackknife deletions and saved-model predictions. Exact distance ties use the
original covariates lexicographically in their supplied column order, then
observation order for identical covariate rows. Outcomes and nuisances do not
enter this ordering; rank weights and effective sample size are preserved.

{title:Stored results}

{pstd}
For q variables in {cmd:d()}, {cmd:plqr_dnn} stores:

{synoptset 25 tabbed}{...}
{p2col 5 25 29 2:Matrices}{p_end}
{synopt:{cmd:e(moment)}}1-by-q row-vector localized moment; columns in d() order{p_end}
{synopt:{cmd:e(moment_vcov)}}q-by-q full moment covariance; rows/columns in d() order{p_end}
{synopt:{cmd:e(theta_eval)}}1-by-q supplied evaluation vector{p_end}
{synopt:{cmd:e(b)}}Stata estimation matrix containing moments{p_end}
{synopt:{cmd:e(V)}}Stata estimation matrix containing moment covariance{p_end}

{pstd}
{cmd:e(moment)} has row name {cmd:moment}. All three primary matrices carry
the regressor names in their component dimension. Stata uses a row vector
for this single target, unlike R's component-by-target moment matrix.

{p2col 5 25 29 2:Scalars}{p_end}
{synopt:{cmd:e(N)}}complete-case estimation-sample size{p_end}
{synopt:{cmd:e(neff)}}DNN effective sample size{p_end}
{synopt:{cmd:e(k)}, {cmd:e(s)}}resolved localization parameters{p_end}
{synopt:{cmd:e(kx)}, {cmd:e(kd)}}reporting and moment dimensions{p_end}
{synopt:{cmd:e(tau)}}quantile index{p_end}

{p2col 5 25 29 2:Macros}{p_end}
{synopt:{cmd:e(theta_source)}}supplied{p_end}
{synopt:{cmd:e(inference_target)}}moment{p_end}
{synopt:{cmd:e(root_solver)}}none{p_end}
{synopt:{cmd:e(xvars)}, {cmd:e(x0)}}reporting covariates and supplied target coordinates{p_end}
{synopt:{cmd:e(yvar)}, {cmd:e(dvars)}}outcome and regressor names{p_end}
{synopt:{cmd:e(gvar)}}shape-variable name{p_end}
{synopt:{cmd:e(mdvars)}}residualizer variable names{p_end}
{synopt:{cmd:e(distance)}, {cmd:e(method)}}metric and rank-weight configuration{p_end}
{synopt:{cmd:e(cmd)}, {cmd:e(cmdline)}}command identity and original call{p_end}

{marker examples}{...}
{title:Examples}

{pstd}
A self-contained median-score example with known shape and regressor means.
These oracle inputs illustrate computation, not a fitted-nuisance guarantee.

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. set seed 20261017}{p_end}
{phang2}{cmd:. gen double x1 = rnormal()}{p_end}
{phang2}{cmd:. gen double x2 = rnormal()}{p_end}
{phang2}{cmd:. gen double md1 = .3*x1-.2*x2}{p_end}
{phang2}{cmd:. gen double md2 = -.1*x1+.4*x2}{p_end}
{phang2}{cmd:. gen double d1 = md1+.4*rnormal()}{p_end}
{phang2}{cmd:. gen double d2 = md2+.4*rnormal()}{p_end}
{phang2}{cmd:. gen double g_tau = .25*x1-.15*x2}{p_end}
{phang2}{cmd:. gen double Y = d1-.5*d2+g_tau+.5*rnormal()}{p_end}
{phang2}{cmd:. plqr_dnn Y, x(x1 x2) x0(0 0) d(d1 d2) g(g_tau) md(md1 md2) thetaeval(1 -.5) tau(.5) s(12)}{p_end}
{phang2}{cmd:. matrix m = e(moment)}{p_end}
{phang2}{cmd:. matrix C = e(moment_vcov)}{p_end}
{phang2}{cmd:. matrix theta_supplied = e(theta_eval)}{p_end}
{phang2}{cmd:. matrix list m}{p_end}
{phang2}{cmd:. matrix list C}{p_end}
{phang2}{cmd:. matrix list theta_supplied}{p_end}
{phang2}{cmd:. display "First moment's jackknife SD: " sqrt(C[1,1])}{p_end}

{pstd}
The source archive's {cmd:examples/plqr-moments.do} also reconstructs the
scores and checks the reported moments against explicit weighted sums.

{marker scope}{...}
{title:R companion and scope}

{pstd}
The R companion's experimental PLQR wrapper has additional density-backed
coefficient inference, multiple targets, tuning grids and PLQR-specific
reporting/plot methods. Those are not Stata options. Neither front end
currently solves a local PLQR root. Stata's full moment covariance is useful
for inspecting score variation, but it does not replace the missing
coefficient Jacobian. Any later coefficient estimator needs a separate,
explicit statistical and software contract.

{title:Also see}

{psee}{help rdistnn}, {help weights_dnn}, {help cate_dnn}, {help np_dnn}, {help clate_dnn}

{title:Saved moment reporting and reuse}
{pstd}
No local root solving was performed. Supplied coefficients remain separate
from evaluated score moments.
{pstd}
{help rdistnn_results} extracts standardized moment tables and supplied theta
separately. {help rdistnn_predict} evaluates the saved frozen scores at new
single/batch targets, optionally returning joint moment covariance. post and
store() create reusable derived fits. Current estimates save/use retains the
numeric model without original active data. Coefficient inference and root
solving remain outside this command's scope.
