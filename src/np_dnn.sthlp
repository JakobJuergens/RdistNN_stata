{smcl}
{* *! version 0.0.1  2026-04-05}{...}
{viewerjumpto "Syntax"         "np_dnn##syntax"}{...}
{viewerjumpto "Description"    "np_dnn##description"}{...}
{viewerjumpto "Options"        "np_dnn##options"}{...}
{viewerjumpto "Stored results" "np_dnn##results"}{...}
{viewerjumpto "Examples"       "np_dnn##examples"}{...}

{title:Title}

{phang}
{bf:np_dnn} {hline 2} Nonparametric regression at a point via
distributional nearest neighbors

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
{synopt:{opth x(varlist)}}covariates{p_end}
{synopt:{opt x0(numlist)}}evaluation point{p_end}

{syntab:Tuning}
{synopt:{opt k(#)}}nearest neighbors; default {cmd:k(5)}{p_end}
{synopt:{opt s(#)}}kernel order; default uses floor(n^((kx+1)/(kx+2))){p_end}
{synopt:{opt distance(string)}}{cmd:l2} (default) or {cmd:mahalanobis}{p_end}
{synopt:{opt level(#)}}confidence level; default {cmd:level(95)}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:np_dnn} estimates E[Y | X = x0] at a fixed point {opt x0} using
DNN rank weights.  This is the building block for {helpb cate_dnn} when
applied to pseudo-outcomes directly.

{marker options}{...}
{title:Options}

{phang}{opth x(varlist)} covariates defining the neighborhood.

{phang}{opt x0(numlist)} evaluation point; one value per covariate.

{phang}{opt k(#)} number of nearest neighbors. Default {cmd:k(5)}.

{phang}{opt s(#)} kernel order. Default uses undersmoothing rule.

{phang}{opt distance(string)} {cmd:l2} or {cmd:mahalanobis}.

{marker results}{...}
{title:Stored results}

{pstd}{cmd:np_dnn} stores the following in {cmd:e()}:

{synoptset 15}{...}
{synopt:{cmd:e(b)}}point estimate (1x1: mu_hat = E[Y|X=x0]){p_end}
{synopt:{cmd:e(V)}}jackknife variance{p_end}
{synopt:{cmd:e(N)}}sample size{p_end}
{synopt:{cmd:e(neff)}}effective sample size{p_end}
{synopt:{cmd:e(k)}}{cmd:k}{p_end}
{synopt:{cmd:e(s)}}{cmd:s}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:np_dnn}{p_end}

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. np_dnn y, x(x1 x2) x0(0 0)}{p_end}
{phang2}{cmd:. np_dnn y, x(x1 x2) x0(0 0) k(10) s(200)}{p_end}
