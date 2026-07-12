{smcl}
{* *! version 0.0.5  2026-07-12}{...}
{viewerjumpto "Syntax"         "weights_dnn##syntax"}{...}
{viewerjumpto "Description"    "weights_dnn##description"}{...}
{viewerjumpto "Stored results" "weights_dnn##results"}{...}
{viewerjumpto "Examples"       "weights_dnn##examples"}{...}

{title:Title}

{phang}
{bf:weights_dnn} {hline 2} Compute DNN rank weights at a fixed point

{pstd}
Stata 17 or later is required.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:weights_dnn}
{it:x1} [{it:x2} {it:...}]
{ifin}{cmd:,}
{opt x0(numlist)}
[{opt k(#)}]
[{opt s(#)}
{opt generate(newvar)}
{opt distance(string)}
{opt method(string)}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:weights_dnn} computes the DNN rank weights w_i(x0) for each
observation.  These are the core localization weights used by
{helpb cate_dnn}, {helpb np_dnn}, and {helpb plqr_dnn}.

{pstd}
A single covariate is permitted when using {cmd:distance(l2)} (the default).
{cmd:distance(mahalanobis)} requires at least 2 covariates to form a
covariance matrix.

{pstd}
Use {opt generate()} to store the weights as a new variable for
inspection or downstream use.

{pstd}
For exact top-k DNN weights, positive support is global ranks 1 through
n - s + k, so the support size is n - s + k.  With the paper-default
{cmd:k(1)}, the edge and support size are n - s + 1.  The parameter k counts
the nearest observations selected within each size-s subsample; it is not the
full-sample positive-support size.

{pstd}
The order s is a localization or bandwidth analogue whose direction is
opposite to the usual bandwidth convention.  Holding n and k fixed, larger s
shrinks positive support and gives tighter localization.  Smaller s spreads
weight across more observations and increases the effective local sample size.

{marker options}{...}
{title:Options}

{phang}{opt x0(numlist)} evaluation point; one value per covariate.

{phang}{opt k(#)} number selected within each size-s subsample (integer >= 1);
default {cmd:k(1)}.

{phang}{opt s(#)} subset/kernel order; must satisfy k <= s <= n.  The endpoint
{cmd:s(n)} is valid because {cmd:weights_dnn} does not compute a jackknife.
Default {cmd:s(0)} uses the paper's finite-sample rule
floor(0.4*n^((kx+1)/(kx+2))).

{phang}{opt generate(newvar)} name of new variable to store the weights.

{phang}{opt distance(string)} {cmd:l2} (default, Euclidean) or
{cmd:mahalanobis} (ellipsoidal, requires >= 2 covariates).

{phang}{opt method(string)} {cmd:exact} (default) or {cmd:asymp}.

{marker results}{...}
{title:Stored results}

{pstd}{cmd:weights_dnn} stores the following in {cmd:r()}:

{synoptset 12}{...}
{synopt:{cmd:r(neff)}}effective sample size 1/sum(w_i^2){p_end}
{synopt:{cmd:r(n)}}sample size{p_end}
{synopt:{cmd:r(k)}}{cmd:k}{p_end}
{synopt:{cmd:r(s)}}{cmd:s}{p_end}
{synopt:{cmd:r(kx)}}covariate dimension{p_end}

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. weights_dnn x1 x2, x0(0 0) k(5) generate(w)}{p_end}
{phang2}{cmd:. summarize w}{p_end}
{phang2}{cmd:. weights_dnn x1, x0(0) k(5) generate(w)}{p_end}
{phang2}{cmd:. weights_dnn x1 x2, x0(0 0) k(5) distance(mahalanobis)}{p_end}
