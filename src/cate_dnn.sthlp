{smcl}
{* *! version 0.0.2  2026-04-06}{...}
{viewerjumpto "Syntax"          "cate_dnn##syntax"}{...}
{viewerjumpto "Description"     "cate_dnn##description"}{...}
{viewerjumpto "Options"         "cate_dnn##options"}{...}
{viewerjumpto "Stored results"  "cate_dnn##results"}{...}
{viewerjumpto "Examples"        "cate_dnn##examples"}{...}
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
{it:outcome}
{ifin}{cmd:,}
{opth treat(varname)}
{opth x(varlist)}
{opt x0(numlist)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth treat(varname)}}binary treatment indicator (0/1){p_end}
{synopt:{opth x(varlist)}}covariates used to define the neighborhood{p_end}
{synopt:{opt x0(numlist)}}evaluation point; one value per covariate{p_end}

{syntab:Nuisances (propensity and both outcome models required)}
{synopt:{opth propensity(varname)}}estimated or oracle propensity score e(X){p_end}
{synopt:{opth m0(varname)}}estimated or oracle E[Y|X, D=0]{p_end}
{synopt:{opth m1(varname)}}estimated or oracle E[Y|X, D=1]{p_end}

{syntab:Tuning}
{synopt:{opt k(#)}}number of nearest neighbors; default {cmd:k(5)}{p_end}
{synopt:{opt s(#)}}kernel order; default {cmd:s(0)} uses
    floor(n^((kx+1)/(kx+2))){p_end}
{synopt:{opt distance(string)}}{cmd:l2} (default) or {cmd:mahalanobis}{p_end}

{syntab:Reporting}
{synopt:{opt level(#)}}set confidence level; default {cmd:level(95)}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:cate_dnn} estimates the conditional average treatment effect
tau(x0) = E[Y(1) - Y(0) | X = x0] at a fixed covariate point {opt x0}
using the distributional nearest-neighbor (DNN) estimator of Juergens (2025).

{pstd}
The estimator applies DNN rank weights w_i(x0) to per-observation
pseudo-outcomes, yielding a localized weighted mean.  The doubly-robust
AIPW pseudo-outcome is

{phang2}phi_i = D_i*(Y_i-m1_i)/eclip_i - (1-D_i)*(Y_i-m0_i)/(1-eclip_i) + (m1_i-m0_i){p_end}

{pstd}
where eclip_i = min(max(e_i, 1e-8), 1-1e-8).  The command requires
{opt propensity()}, {opt m0()}, and {opt m1()}.

{pstd}
Standard errors are computed via the rank-shift jackknife.  For
{cmd:distance(l2)}, this is the exact delete-1 recomputation because the
remaining distances do not change.  For {cmd:distance(mahalanobis)}, the
sample covariance metric is computed from the full estimation sample and
held fixed across jackknife deletions.

{pstd}
Nuisance functions (propensity score and outcome models) must be supplied
by the user as variables in the dataset.  They may be oracle values or
first-stage estimates computed on the same sample (Neyman orthogonality
controls first-stage bias) or from a cross-fitting procedure.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth treat(varname)} specifies the binary treatment indicator D (0 or 1).

{phang}
{opth x(varlist)} specifies the covariates X used to define the local
neighborhood around x0.

{phang}
{opt x0(numlist)} specifies the evaluation point.  Provide one number per
variable in {opt x()}, separated by spaces.  Example: {cmd:x0(0 0)} for
the two-dimensional origin.

{dlgtab:Nuisances}

{phang}
{opth propensity(varname)} provides the propensity score e(X) = P(D=1|X).
Values must lie in [0, 1].  Values at or very near the endpoints are clipped
only inside the residual denominators.

{phang}
{opth m0(varname)} provides the conditional outcome mean E[Y|X, D=0].
This option is required.

{phang}
{opth m1(varname)} provides the conditional outcome mean E[Y|X, D=1].
This option is required.

{dlgtab:Tuning}

{phang}
{opt k(#)} sets the number of nearest neighbors used in the DNN weights.
Default is {cmd:k(5)}.

{phang}
{opt s(#)} sets the kernel order (bandwidth analog).  Larger s includes
more observations with nonzero weight.  Default {cmd:s(0)} uses the
undersmoothing rule floor(n^((kx+1)/(kx+2))).

{phang}
{opt distance(string)} selects the distance metric.  {cmd:l2} (default)
uses Euclidean distance.  {cmd:mahalanobis} uses the sample covariance
matrix (via {cmd:invsym()}) to form an ellipsoidal neighborhood;
recommended when covariates are strongly correlated
({it:cf.} DGP 5 of Juergens 2025).  Requires at least 2 covariates.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:cate_dnn} stores the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(neff)}}effective sample size 1/sum(w_i^2){p_end}
{synopt:{cmd:e(k)}}neighbor count k{p_end}
{synopt:{cmd:e(s)}}kernel order s{p_end}
{synopt:{cmd:e(kx)}}covariate dimension{p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}point estimate (1x1: tau_hat){p_end}
{synopt:{cmd:e(V)}}jackknife variance (1x1){p_end}

{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:cate_dnn}{p_end}
{synopt:{cmd:e(x0)}}evaluation point{p_end}
{synopt:{cmd:e(xvars)}}covariate variable names{p_end}
{synopt:{cmd:e(distance)}}distance metric used{p_end}

{marker examples}{...}
{title:Examples}

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
{phang2}{cmd:. cate_dnn Y, treat(D) x(x1 x2) x0(0 0)}{...}
{phang2}{cmd:    propensity(e_oracle) m0(m0_true) m1(m1_true)}{p_end}

{pstd}With Mahalanobis distance (for correlated covariates).{p_end}
{phang2}{cmd:. cate_dnn Y, treat(D) x(x1 x2) x0(0.75 0)}{...}
{phang2}{cmd:    propensity(e_oracle) m0(m0_true) m1(m1_true)}{...}
{phang2}{cmd:    distance(mahalanobis)}{p_end}

{marker references}{...}
{title:References}

{phang}
Juergens, J. (2025).
{it:Robust Inference for Conditional Z-Estimation via Distributional Nearest Neighbors}.
Working paper.
