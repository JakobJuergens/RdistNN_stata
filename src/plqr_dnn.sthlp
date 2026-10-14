{smcl}
{* *! version 0.1.0  2026-07-12}{...}
{viewerjumpto "Syntax"          "plqr_dnn##syntax"}{...}
{viewerjumpto "Description"     "plqr_dnn##description"}{...}
{viewerjumpto "Options"         "plqr_dnn##options"}{...}
{viewerjumpto "Stored results"  "plqr_dnn##results"}{...}
{viewerjumpto "Examples"        "plqr_dnn##examples"}{...}
{viewerjumpto "Future work"     "plqr_dnn##future"}{...}

{title:Title}

{phang}
{bf:plqr_dnn} {hline 2} Experimental localized PLQR moment evaluation at a
supplied coefficient vector

{pstd}
Stata 17 or later is required.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:plqr_dnn}
{it:outcome}
{ifin}{cmd:,}
{opth x(varlist)}
{opt x0(numlist)}
{opth d(varlist)}
{opth g(varname)}
{opth md(varlist)}
{opt thetaeval(numlist)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth x(varlist)}}covariates defining the local neighborhood{p_end}
{synopt:{opt x0(numlist)}}evaluation point; one value per covariate{p_end}
{synopt:{opth d(varlist)}}partially linear regressors D{p_end}
{synopt:{opth g(varname)}}supplied g_tau(X) nuisance variable{p_end}
{synopt:{opth md(varlist)}}supplied m_D(X) variables, in d() order{p_end}
{synopt:{opt thetaeval(numlist)}}supplied theta vector, one value per d() variable{p_end}

{syntab:Compatibility aliases}
{synopt:{opth q(varname)}}temporary alias for {cmd:g()}; do not combine them{p_end}
{synopt:{opt beta(numlist)}}temporary alias for {cmd:thetaeval()}; do not combine them{p_end}

{syntab:Model and localization}
{synopt:{opt tau(#)}}quantile index in (0,1); default {cmd:tau(0.5)}{p_end}
{synopt:{opt k(#)}}nearest observations selected per size-s subsample; default {cmd:k(1)}{p_end}
{synopt:{opt s(#)}}kernel order (< n); default uses floor(0.4*n^((kx+1)/(kx+2))){p_end}
{synopt:{opt distance(string)}}{cmd:l2} (default) or {cmd:mahalanobis}{p_end}
{synopt:{opt method(string)}}{cmd:exact} (default) or {cmd:asymp}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{bf:Experimental status.}  {cmd:plqr_dnn} evaluates the localized PLQR score
moment at the vector supplied in {opt thetaeval()}.  It is not a coefficient estimator.
It does not solve the local moment equation for a root,
estimate a conditional density, or report coefficient standard errors,
confidence intervals, p-values, or significance tests.

{pstd}
PLQR functionality is outside the CATE and nonparametric-regression focus of
the initial public release. Its API may change in later releases. The Stata
command's implemented contract is moment evaluation at supplied coefficients;
the broader R companion also has density-backed inference facilities that
are not provided by this command.

{pstd}
For observation i and score component j, the command forms

{phang2}
psi_ij(theta) = (tau - 1{Y_i <= D_i'theta + g_i}) * (D_ij - mD_ij),
{p_end}

{pstd}
then computes its DNN-weighted local mean at {opt x0()}.  The posted covariance
is the rank-shift jackknife covariance of this localized score moment.  It is
not a coefficient covariance matrix.  The displayed output explicitly states
that no local root solving was performed.

{pstd}
The canonical paper-aligned nuisance notation is {opt g()} for g_tau(X) and
{opt md()} for m_D(X).  Stata does not fit these nuisance functions.  The
canonical evaluation-vector option is {opt thetaeval()}.

{marker options}{...}
{title:Options}

{dlgtab:Moment inputs}

{phang}{opth x(varlist)} specifies the localization covariates X.

{phang}{opt x0(numlist)} supplies the target point, one value per variable in
{opt x()}.

{phang}{opth d(varlist)} supplies the partially linear regressors D.

{phang}{opth g(varname)} supplies g_tau(X), the additive nuisance entering the
indicator threshold.  {opt q()} is a temporary compatibility alias.  Supplying
both {opt g()} and {opt q()} is an error.

{phang}{opth md(varlist)} supplies m_D(X).  Its variables must match {opt d()}
in length and order.

{phang}{opt thetaeval(numlist)} supplies the coefficient vector at which the
score is evaluated.  {opt beta()} is a temporary compatibility alias.
Supplying both {opt thetaeval()} and {opt beta()} is an error.

{phang}{opt tau(#)} sets the quantile index in (0,1).

{dlgtab:Localization}

{phang}{opt k(#)} sets how many nearest observations are selected within each
size-s subsample, not the full-sample positive-support size.

{phang}{opt s(#)} sets the subset/kernel order, a localization analogue whose
direction is opposite to the usual bandwidth convention.  Exact top-k weights
have positive support through rank n - s + k (n - s + 1 for {cmd:k(1)}).
Larger s therefore tightens localization.  Because the moment covariance uses
a delete-1 jackknife, {cmd:s()} must satisfy k <= s < n.

{phang}{opt distance(string)} selects {cmd:l2} or {cmd:mahalanobis}.

{phang}{opt method(string)} selects exact DNN weights or the {cmd:asymp} Beta
order-statistic approximation.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:plqr_dnn} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2:Matrices}{p_end}
{synopt:{cmd:e(moment)}}localized score moment at supplied theta{p_end}
{synopt:{cmd:e(moment_vcov)}}rank-shift covariance of {cmd:e(moment)}{p_end}
{synopt:{cmd:e(theta_eval)}}supplied evaluation vector{p_end}
{synopt:{cmd:e(b)}}compatibility copy of {cmd:e(moment)} required by eclass{p_end}
{synopt:{cmd:e(V)}}compatibility copy of {cmd:e(moment_vcov)}{p_end}
{synopt:{cmd:e(beta)}}temporary compatibility copy of {cmd:e(theta_eval)}{p_end}

{p2col 5 22 26 2:Scalars}{p_end}
{synopt:{cmd:e(N)}}estimation-sample size{p_end}
{synopt:{cmd:e(neff)}}DNN effective sample size{p_end}
{synopt:{cmd:e(k)}, {cmd:e(s)}}DNN localization parameters{p_end}
{synopt:{cmd:e(kx)}, {cmd:e(kd)}}covariate and moment dimensions{p_end}
{synopt:{cmd:e(tau)}}quantile index{p_end}

{p2col 5 22 26 2:Macros}{p_end}
{synopt:{cmd:e(theta_source)}}{cmd:supplied}{p_end}
{synopt:{cmd:e(inference_target)}}{cmd:moment}{p_end}
{synopt:{cmd:e(root_solver)}}{cmd:none}{p_end}
{synopt:{cmd:e(gvar)}}resolved g_tau(X) variable{p_end}
{synopt:{cmd:e(g_option)}}{cmd:g} or compatibility alias {cmd:q}{p_end}
{synopt:{cmd:e(theta_option)}}{cmd:thetaeval} or compatibility alias {cmd:beta}{p_end}
{synopt:{cmd:e(mdvars)}}resolved m_D(X) variables{p_end}
{synopt:{cmd:e(distance)}, {cmd:e(method)}}weight configuration{p_end}

{marker examples}{...}
{title:Examples}

{pstd}Evaluate the localized score moment at a supplied theta.{p_end}
{phang2}{cmd:. plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) g(g_tau) md(md1 md2)}{...}
{phang2}{cmd:    thetaeval(1 -0.5) tau(0.5) k(5)}{p_end}

{pstd}Legacy aliases remain temporarily available.{p_end}
{phang2}{cmd:. plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(g_tau) md(md1 md2)}{...}
{phang2}{cmd:    beta(1 -0.5)}{p_end}

{marker future}{...}
{title:Future functionality}

{pstd}
Local root solving and density-backed coefficient inference are not implemented
by this experimental command.  They require a separate estimator contract and
are reserved for future functionality.
