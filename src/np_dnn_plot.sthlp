{smcl}
{* *! version 0.1.0  2026-07-12}{...}
{title:Title}

{phang}{bf:np_dnn_plot} {hline 2} Plot NP estimates, diagnostics, tuning, weights, and outcomes

{title:Syntax}

{p 8 17 2}{cmd:np_dnn_plot}
[{cmd:,} {opt type(plot_type)} {opt view(view_type)} {opt style(style)}
{opt xindex(#)} {opt targets(numlist)} {opt candidates(numlist)}
{opt allcandidates} {opt level(#)} {opt logweights} {opt name(name)}
{opt replace} {opt nodraw} {opt saving(filename)}]

{title:Description}

{pstd}
{cmd:np_dnn_plot} provides Stata-native analogues of the R package's NP
estimate, diagnostics, tuning, weight, and dashboard plots.  Plot types are
{cmd:estimates} (default), {cmd:diagnostics}, {cmd:tuning}, {cmd:weights},
{cmd:outcome}, and {cmd:dashboard}.  Stored-result views do not refit.
The effective-kernel weight view validates and reconstructs weights from the
current estimation sample.

{title:Plot types and views}

{phang}{cmd:type(estimates)} plots target estimates and confidence intervals.
Styles are {cmd:auto}, {cmd:profile}, {cmd:pointrange}, and {cmd:forest}.{p_end}
{phang}{cmd:type(diagnostics)} combines ESS, maximum weight, and top-5 share.
Candidate selection switches to tuning diagnostics.{p_end}
{phang}{cmd:type(tuning)} plots available candidate estimates against s and
requires {cmd:tuningmode(grid)}. The built-in grid is lambda =
(.3, .4, .5, .75), with .4 fitted by default.{p_end}
{phang}{cmd:type(weights)} supports {cmd:view(profile)},
{cmd:view(cumulative)}, and {cmd:view(kernel)}.  Candidate comparison is
available for profiles.{p_end}
{phang}{cmd:type(outcome)} plots weighted outcome means with weighted 10th-90th
percentile intervals.{p_end}
{phang}{cmd:type(dashboard)} combines estimates, localization diagnostics,
outcome summaries, weights, and tuning when available.{p_end}

{title:Options}

{phang}{opt xindex(#)} selects a target coordinate for a profile axis.
Zero lets the command detect a single varying coordinate.{p_end}
{phang}{opt targets(numlist)} selects target row numbers.{p_end}
{phang}{opt candidates(numlist)} and {opt allcandidates} select grid rows.{p_end}
{phang}{opt logweights} uses a log scale for positive rank weights.{p_end}
{phang}{opt name(name)}, {opt replace}, and {opt nodraw} follow Stata graph
semantics.  {opt saving(filename)} exports or saves the graph.{p_end}

{title:Stored results}

{synoptset 18}{...}
{synopt:{cmd:r(type)}}plot type{p_end}
{synopt:{cmd:r(view)}}resolved view or style{p_end}
{synopt:{cmd:r(graph)}}graph name{p_end}
{synopt:{cmd:r(n_plotted)}}number of plotted rows, candidates, targets, or panels{p_end}
{synopt:{cmd:r(saved)}}saved filename, if any{p_end}

{title:Examples}

{phang2}{cmd:. np_dnn_plot, type(estimates)}{p_end}
{phang2}{cmd:. np_dnn_plot, type(tuning) allcandidates}{p_end}
{phang2}{cmd:. np_dnn_plot, type(weights) view(kernel) targets(1 3)}{p_end}
{phang2}{cmd:. np_dnn_plot, type(outcome)}{p_end}
{phang2}{cmd:. np_dnn_plot, type(dashboard) saving(np_dashboard.pdf) replace}{p_end}
