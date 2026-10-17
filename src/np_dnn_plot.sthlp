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

{pstd}
All single-view plots return {cmd:r(plot_data)} containing only rows used by
the graph, with canonical target/candidate identifiers. {cmd:r(n_rows)} counts
these data rows; {cmd:r(n_plotted)} retains its existing target, candidate, or
panel count. Empty estimate plots retain their legacy count and return
{cmd:r(n_rows)} = 0 without a plot-data matrix. Undefined CLATE ratios have
no plotted estimate; CLATE zero-variance ratios retain their stored zero SE and have no Wald interval.
Ordinary CATE/NP zero standard errors retain degenerate mean intervals.

{pstd}
{cmd:targets()} and candidate selection can be combined for estimates and
tuning. They are mutually exclusive for diagnostics. Rank-weight profiles and
cumulative plots are common across targets and reject {cmd:targets()}; use
{cmd:view(kernel)} for target-specific weights. Candidate options apply only
to estimates, tuning, diagnostics, and weight profiles. Dashboards reject
selection options. Options that do not affect the requested plot are rejected.

{pstd}
Observation-level kernel plots require unchanged original values,
estimation sample, and row order. Older saved fits without order signatures
must be refitted for these views; stored-result views still work after data
changes. {cmd:style()} applies to estimates and dashboard estimate panels;
{cmd:logweights} applies only to weight profiles.

{title:Stored results}

{synoptset 18}{...}
{synopt:{cmd:r(type)}}plot type{p_end}
{synopt:{cmd:r(view)}}resolved view or style{p_end}
{synopt:{cmd:r(graph)}}graph name{p_end}
{synopt:{cmd:r(n_plotted)}}number of plotted rows, candidates, targets, or panels{p_end}
{synopt:{cmd:r(n_omitted)}}selected estimate/tuning rows omitted because the estimate is undefined{p_end}
{synopt:{cmd:r(selected_targets)}}all selected estimate/tuning target IDs, including undefined targets{p_end}
{synopt:{cmd:r(n_rows)}}number of data rows used by a single-view plot; missing for dashboards{p_end}
{synopt:{cmd:r(targets)}}canonical target IDs present in the plot data, when applicable{p_end}
{synopt:{cmd:r(candidates)}}canonical candidate IDs present in the plot data; 0 denotes the fitted default{p_end}
{synopt:{cmd:r(n_panels)}}dashboard panel count{p_end}
{synopt:{cmd:r(panels)}}dashboard panel names in composition order{p_end}
{synopt:{cmd:r(panel_info)}}dashboard rows: panel_id, n_rows, n_plotted; no common plot-data matrix{p_end}
{synopt:{cmd:r(plot_data)}}selected data for estimates, tuning, diagnostics, weights, and other family-specific views{p_end}
{synopt:{cmd:r(saved)}}saved filename, if any{p_end}

{pstd}
Comparable estimate columns are target coordinates, candidate_id, s, is_default,
target_id, estimate, se, ci_lower, ci_upper. Tuning keeps its existing
candidate_id, s, lambda, is_default, target_id, estimate, se, ci_lower, ci_upper
schema. Diagnostic data retain plot_id and add target_id and candidate_id;
rank weights expose rank, weight, cum_share, candidate_id, s, is_default.
Kernel data expose target_id, rank, distance, weight, cum_share. Overlap names
retain the treatment/instrument distinction. NP outcome columns retain weighted
means, standard deviations, and weighted 10th/50th/90th percentiles.

{title:Examples}

{phang2}{cmd:. np_dnn_plot, type(estimates)}{p_end}
{phang2}{cmd:. np_dnn_plot, type(tuning) allcandidates}{p_end}
{phang2}{cmd:. np_dnn_plot, type(weights) view(kernel) targets(1 3)}{p_end}
{phang2}{cmd:. np_dnn_plot, type(outcome)}{p_end}
{phang2}{cmd:. np_dnn_plot, type(dashboard) saving(np_dashboard.pdf) replace}{p_end}
