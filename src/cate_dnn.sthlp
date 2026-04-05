{smcl}
{* *! version 0.0.1  2026-04-05}{...}
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

{syntab:Nuisances (at least propensity + one outcome model required)}
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
The estimator applies DNN rank weights w_i(x0) to the Robinson (1988)
pseudo-outcomes rho_i = (D_i - e_i)(Y_i - m_i) / (e_i(1-e_i)),
yielding a localized weighted mean.  Standard errors are computed via
the fast leave-one-out jackknife.  The DNN weights exploit the
generalized-U-statistic / L-statistic structure to deliver valid
same-sample inference without cross-fitting.

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
Values should lie strictly in (0, 1).

{phang}
{opth m0(varname)} provides the conditional outcome mean E[Y|X, D=0].

{phang}
{opth m1(varname)} provides the conditional outcome mean E[Y|X, D=1].

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
matrix to form an ellipsoidal neighborhood; recommended when covariates
are strongly correlated ({it:cf.} DGP 5 of Juergens 2025).

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
{phang2}{cmd:. gen Y = m0_true + D*tau_true + rnormal()*0.5}{p_end}

{pstd}Estimate CATE at the origin using oracle nuisances.{p_end}
{phang2}{cmd:. cate_dnn Y, treat(D) x(x1 x2) x0(0 0)}{...}
{phang2}{cmd:    propensity(e_oracle) m0(m0_true) m1(m0_true)}{p_end}

{pstd}With Mahalanobis distance (for correlated covariates).{p_end}
{phang2}{cmd:. cate_dnn Y, treat(D) x(x1 x2) x0(0.75 0)}{...}
{phang2}{cmd:    propensity(e_oracle) m0(m0_true) m1(m0_true)}{...}
{phang2}{cmd:    distance(mahalanobis)}{p_end}

{marker references}{...}
{title:References}

{phang}
Juergens, J. (2025).
{it:Robust Inference for Conditional Z-Estimation via Distributional Nearest Neighbors}.
Working paper.

{phang}
Robinson, P. M. (1988).
Root-N-consistent semiparametric regression.
{it:Econometrica}, 56(4), 931-954.
