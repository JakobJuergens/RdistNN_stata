{smcl}
{* *! version 0.0.1  2026-05-30}{...}
{viewerjumpto "Syntax"          "plqr_dnn##syntax"}{...}
{viewerjumpto "Description"     "plqr_dnn##description"}{...}
{viewerjumpto "Options"         "plqr_dnn##options"}{...}
{viewerjumpto "Stored results"  "plqr_dnn##results"}{...}
{viewerjumpto "Examples"        "plqr_dnn##examples"}{...}
{viewerjumpto "References"      "plqr_dnn##references"}{...}

{title:Title}

{phang}
{bf:plqr_dnn} {hline 2} Local PLQR orthogonal moments via
distributional nearest neighbors

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
{opth q(varname)}
{opth md(varlist)}
{opt beta(numlist)}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth x(varlist)}}covariates defining the local neighborhood{p_end}
{synopt:{opt x0(numlist)}}evaluation point; one value per covariate{p_end}
{synopt:{opth d(varlist)}}partially linear regressors D{p_end}
{synopt:{opth q(varname)}}supplied q_tau(X) nuisance variable{p_end}
{synopt:{opth md(varlist)}}supplied m_D(X) nuisance variables, in d() order{p_end}
{synopt:{opt beta(numlist)}}fixed beta vector, one value per d() variable{p_end}

{syntab:Model}
{synopt:{opt tau(#)}}quantile index in (0,1); default {cmd:tau(0.5)}{p_end}

{syntab:Tuning}
{synopt:{opt k(#)}}nearest neighbors (>= 1); default {cmd:k(5)}{p_end}
{synopt:{opt s(#)}}kernel order (< n); default uses floor(n^((kx+1)/(kx+2))){p_end}
{synopt:{opt distance(string)}}{cmd:l2} (default) or {cmd:mahalanobis}{p_end}
{synopt:{opt method(string)}}{cmd:exact} (default) or {cmd:asymp}{p_end}
{synopt:{opt level(#)}}confidence level; default {cmd:level(95)}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:plqr_dnn} computes localized partially linear quantile regression
(PLQR) orthogonal moments at a supplied coefficient vector {opt beta()}.
It does not estimate nuisance functions, solve for beta, or perform GMM.

{pstd}
For each observation i and component j, the command forms

{phang2}
g_ij = (tau - 1{Y_i <= D_i'beta + q_i}) * (D_ij - md_ij),
{p_end}

{pstd}
then reports the DNN-weighted local mean of each moment component at {opt x0}.
The variables in {opt q()} and {opt md()} are supplied nuisance estimates or
oracle values.  The variables in {opt md()} must have the same length and
order as {opt d()}.

{pstd}
Standard errors and covariances are computed via the rank-shift jackknife.
For {cmd:distance(l2)}, this is the exact delete-1 recomputation because
remaining distances do not change.  For {cmd:distance(mahalanobis)}, the
sample covariance metric is computed from the full estimation sample and
held fixed across jackknife deletions.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth x(varlist)} specifies the covariates X used to define the local
neighborhood around x0.

{phang}
{opt x0(numlist)} specifies the evaluation point.  Provide one number per
variable in {opt x()}.

{phang}
{opth d(varlist)} specifies the partially linear regressors D.

{phang}
{opth q(varname)} supplies the q_tau(X) nuisance variable.

{phang}
{opth md(varlist)} supplies the m_D(X) nuisance variables.  The list must
match {opt d()} in length and order.

{phang}
{opt beta(numlist)} supplies the fixed coefficient vector used in the
indicator threshold D_i'beta + q_i.  It must contain one number per variable
in {opt d()}.

{dlgtab:Model}

{phang}
{opt tau(#)} sets the quantile index.  It must lie strictly between 0 and 1.

{dlgtab:Tuning}

{phang}
{opt k(#)} sets the number of nearest neighbors used in the DNN weights.
Default is {cmd:k(5)}.

{phang}
{opt s(#)} sets the kernel order.  Default {cmd:s(0)} uses the undersmoothing
rule floor(n^((kx+1)/(kx+2))).  Because jackknife covariances are posted,
{cmd:s()} must satisfy k <= s < n.

{phang}
{opt distance(string)} selects {cmd:l2} or {cmd:mahalanobis}.  The
Mahalanobis metric uses the full estimation-sample covariance and requires
at least two covariates.

{phang}
{opt method(string)} selects {cmd:exact} DNN rank weights or the
{cmd:asymp} Beta order-statistic approximation.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:plqr_dnn} stores the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}1 x p_d row vector of localized PLQR moments{p_end}
{synopt:{cmd:e(V)}}p_d x p_d jackknife covariance matrix for the moments{p_end}
{synopt:{cmd:e(beta)}}1 x p_d row vector supplied in {opt beta()}{p_end}

{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(neff)}}effective sample size 1/sum(w_i^2){p_end}
{synopt:{cmd:e(k)}}neighbor count k{p_end}
{synopt:{cmd:e(s)}}kernel order s{p_end}
{synopt:{cmd:e(kx)}}covariate dimension{p_end}
{synopt:{cmd:e(kd)}}number of PLQR moment components{p_end}
{synopt:{cmd:e(tau)}}quantile index{p_end}

{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:plqr_dnn}{p_end}
{synopt:{cmd:e(cmdline)}}command line used for estimation{p_end}
{synopt:{cmd:e(x0)}}evaluation point{p_end}
{synopt:{cmd:e(xvars)}}covariate variable names{p_end}
{synopt:{cmd:e(dvars)}}regressor variable names{p_end}
{synopt:{cmd:e(mdvars)}}m_D(X) nuisance variable names{p_end}
{synopt:{cmd:e(qvar)}}q_tau(X) nuisance variable name{p_end}
{synopt:{cmd:e(yvar)}}outcome variable name{p_end}
{synopt:{cmd:e(distance)}}distance metric used{p_end}
{synopt:{cmd:e(method)}}weight method used{p_end}

{marker examples}{...}
{title:Examples}

{pstd}Simulate a small PLQR-style dataset with oracle nuisances.{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 1000}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. gen x1 = rnormal()}{p_end}
{phang2}{cmd:. gen x2 = rnormal()}{p_end}
{phang2}{cmd:. gen md1 = 0.3*x1 - 0.2*x2}{p_end}
{phang2}{cmd:. gen md2 = -0.1*x1 + 0.4*x2}{p_end}
{phang2}{cmd:. gen d1 = md1 + rnormal()*0.4}{p_end}
{phang2}{cmd:. gen d2 = md2 + rnormal()*0.4}{p_end}
{phang2}{cmd:. gen q = 0.25*x1 - 0.15*x2}{p_end}
{phang2}{cmd:. gen y = d1 - 0.5*d2 + q + rnormal()*0.5}{p_end}

{pstd}Compute localized PLQR moments at the origin for a supplied beta.{p_end}
{phang2}{cmd:. plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2)}{...}
{phang2}{cmd:    beta(1 -0.5) tau(0.5) k(5)}{p_end}

{pstd}With Mahalanobis distance.{p_end}
{phang2}{cmd:. plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2)}{...}
{phang2}{cmd:    beta(1 -0.5) distance(mahalanobis)}{p_end}

{marker references}{...}
{title:References}

{phang}
Juergens, J. (2025).
{it:Robust Inference for Conditional Z-Estimation via Distributional Nearest Neighbors}.
Working paper.
