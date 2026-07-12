{smcl}
{* *! version 0.0.3  2026-07-12}{...}
{viewerjumpto "Syntax"         "np_dnn##syntax"}{...}
{viewerjumpto "Description"    "np_dnn##description"}{...}
{viewerjumpto "Options"        "np_dnn##options"}{...}
{viewerjumpto "Stored results" "np_dnn##results"}{...}
{viewerjumpto "Examples"       "np_dnn##examples"}{...}

{title:Title}

{phang}
{bf:np_dnn} {hline 2} Nonparametric regression at a point via
distributional nearest neighbors

{pstd}
Stata 17 or later is required.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:np_dnn}
{it:outcome}
{ifin}{cmd:,}
{opth x(varlist)}
{opt x0(numlist)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth x(varlist)}}covariates (one or more){p_end}
{synopt:{opt x0(numlist)}}evaluation point; one value per covariate{p_end}

{syntab:Tuning}
{synopt:{opt k(#)}}nearest neighbors (>= 1); default {cmd:k(1)}{p_end}
{synopt:{opt s(#)}}kernel order (< n); default uses floor(n^((kx+1)/(kx+2))){p_end}
{synopt:{opt distance(string)}}{cmd:l2} (default) or {cmd:mahalanobis}{p_end}
{synopt:{opt method(string)}}{cmd:exact} (default) or {cmd:asymp}{p_end}
{synopt:{opt level(#)}}confidence level; default {cmd:level(95)}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:np_dnn} estimates E[Y | X = x0] at a fixed covariate point {opt x0}
using DNN rank weights.  This is the building block for {helpb cate_dnn}
when applied to pseudo-outcomes directly.

{pstd}
Standard errors are computed via the rank-shift jackknife.  For
{cmd:distance(l2)}, this is the exact delete-1 recomputation because the
remaining distances do not change.  For {cmd:distance(mahalanobis)}, the
sample covariance metric is computed from the full estimation sample and
held fixed across jackknife deletions.

{marker options}{...}
{title:Options}

{phang}{opth x(varlist)} covariates defining the neighborhood.
A single covariate is allowed with {cmd:distance(l2)}.

{phang}{opt x0(numlist)} evaluation point; one value per covariate.

{phang}{opt k(#)} number of nearest neighbors (integer >= 1).
Default {cmd:k(1)}, matching the paper theory and the R package.

{phang}{opt s(#)} kernel order; must satisfy k <= s < n.
Default uses undersmoothing rule floor(n^((kx+1)/(kx+2))).

{phang}{opt distance(string)} {cmd:l2} (default, Euclidean) or
{cmd:mahalanobis} (ellipsoidal, requires kx >= 2).

{phang}{opt method(string)} {cmd:exact} (default, hypergeometric formula)
or {cmd:asymp} (Beta order-statistic approximation).

{marker results}{...}
{title:Stored results}

{pstd}{cmd:np_dnn} stores the following in {cmd:e()}:

{synoptset 15}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}point estimate (1x1: mu_hat = E[Y|X=x0]){p_end}
{synopt:{cmd:e(V)}}delete-1 jackknife variance{p_end}

{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}sample size{p_end}
{synopt:{cmd:e(neff)}}effective sample size 1/sum(w_i^2){p_end}
{synopt:{cmd:e(k)}}{cmd:k}{p_end}
{synopt:{cmd:e(s)}}{cmd:s}{p_end}
{synopt:{cmd:e(kx)}}covariate dimension{p_end}

{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:np_dnn}{p_end}
{synopt:{cmd:e(xvars)}}covariate variable names{p_end}
{synopt:{cmd:e(distance)}}distance metric used{p_end}
{synopt:{cmd:e(method)}}weight method used{p_end}

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. np_dnn y, x(x1 x2) x0(0 0)}{p_end}
{phang2}{cmd:. np_dnn y, x(x1 x2) x0(0 0) k(10) s(200)}{p_end}
{phang2}{cmd:. np_dnn y, x(x1) x0(0) k(5)}{p_end}
{phang2}{cmd:. np_dnn y, x(x1 x2) x0(0 0) distance(mahalanobis)}{p_end}
