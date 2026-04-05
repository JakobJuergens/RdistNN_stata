{smcl}
{* *! version 0.0.1  2026-04-05}{...}
{viewerjumpto "Syntax"         "weights_dnn##syntax"}{...}
{viewerjumpto "Description"    "weights_dnn##description"}{...}
{viewerjumpto "Stored results" "weights_dnn##results"}{...}
{viewerjumpto "Examples"       "weights_dnn##examples"}{...}

{title:Title}

{phang}
{bf:weights_dnn} {hline 2} Compute DNN rank weights at a fixed point

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:weights_dnn}
{it:x1 x2 ...}
{ifin}{cmd:,}
{opt x0(numlist)}
{opt k(#)}
[{opt s(#)}
{opt generate(newvar)}
{opt distance(string)}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:weights_dnn} computes the DNN rank weights w_i(x0) for each
observation.  These are the core localization weights used by
{helpb cate_dnn} and {helpb np_dnn}.

{pstd}
Use {opt generate()} to store the weights as a new variable for
inspection or downstream use.

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
