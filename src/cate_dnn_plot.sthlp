{smcl}
{* *! version 0.1.1  2026-07-12}{...}
{viewerjumpto "Syntax" "cate_dnn_plot##syntax"}{...}
{viewerjumpto "Description" "cate_dnn_plot##description"}{...}
{viewerjumpto "Plot types" "cate_dnn_plot##types"}{...}
{viewerjumpto "Options" "cate_dnn_plot##options"}{...}
{viewerjumpto "Stored results" "cate_dnn_plot##results"}{...}
{viewerjumpto "Examples" "cate_dnn_plot##examples"}{...}

{title:Title}

{phang}
{bf:cate_dnn_plot} {hline 2} Plot CATE estimates, localization diagnostics,
tuning sensitivity, overlap, weights, and influence

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:cate_dnn_plot}
[{cmd:,} {opt type(plot_type)} {opt view(view_type)} {opt style(style)}
{opt targets(numlist)} {opt candidates(numlist)} {opt allcandidates}
{opt xindex(#)} {opt level(#)} {opt sharelimits(# #)} {opt nlabel(#)}
{opt logweights} {opt name(name)} {opt replace} {opt nodraw}
{opt saving(filename)}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:cate_dnn_plot} is the Stata-native plotting layer for active
{help cate_dnn} estimates.  It covers the numerical views exposed by the R
package's {cmd:plot_estimates()}, diagnostics plot, {cmd:plot_tuning()},
{cmd:plot_overlap()}, {cmd:plot_weights()}, {cmd:plot_influence()}, and
{cmd:plot_dashboard()} methods.

{pstd}
Estimate, diagnostics, tuning, and overlap plots use only stored {cmd:e()}
results.  Effective-kernel and influence plots reconstruct observation-level
weights from the current estimation data.  They verify that the original
covariates, score variables, and estimation sample are still available and
fail rather than silently use stale data.

{marker types}{...}
{title:Plot types}

{phang}
{cmd:type(estimates)} (default) plots CATE estimates with normal confidence
bounds and a zero reference.  {opt style(auto)} chooses a profile for at least
four targets varying in one coordinate, a forest layout for small non-profile
target sets, and pointranges otherwise.  Candidate comparisons consume the
stored tuning grid without refitting.

{phang}
{cmd:type(diagnostics)} combines effective sample size, maximum weight, and
top-five weight-share panels.  Add {opt allcandidates} or {opt candidates()}
to plot tuning-candidate diagnostics against kernel order {it:s}.  For target
plots, the command also annotates the 90-percent physical-radius and top-five
absolute score-contribution-share ranges.

{phang}
{cmd:type(tuning)} plots candidate estimates against {it:s}, marks the fitted
order, and reuses stored candidate standard errors.  It requires a grid fit;
the built-in sensitivity grid is lambda = (.3, .4, .5, .75), with .4 fitted
by default.

{phang}
{cmd:type(overlap)} supports {cmd:view(mass)} (default) for treated/control
weighted mass, {cmd:view(armess)} for normalized within-arm treated/control
ESS, and {cmd:view(propensity)} for weighted propensity means with q10-q90
ranges.  These views require a nuisance-based AIPW or Robinson fit.

{phang}
{cmd:type(weights)} supports {cmd:view(profile)} (default),
{cmd:view(cumulative)}, and {cmd:view(kernel)}.  Candidate comparison is
available for rank profiles.  Kernel views show target-specific weight against
distance and require the original estimation data.

{phang}
{cmd:type(influence)} plots each observation's localized contribution
{it:w_ij(score_i-tau_j)} against distance and labels the largest absolute
contributions.

{phang}
{cmd:type(dashboard)} combines estimates, available overlap panels,
cumulative weights, and tuning sensitivity into a two-column graph.

{marker options}{...}
{title:Options}

{phang}
{opt type(plot_type)} selects {cmd:estimates}, {cmd:diagnostics},
{cmd:tuning}, {cmd:overlap}, {cmd:weights}, {cmd:influence}, or
{cmd:dashboard}.

{phang}
{opt view(view_type)} selects the overlap or weight view described above.
For overlap, choose {cmd:mass}, {cmd:armess}, or {cmd:propensity}.

{phang}
{opt style(style)} selects {cmd:auto}, {cmd:profile}, {cmd:pointrange}, or
{cmd:forest} for estimate panels.

{phang}
{opt xindex(#)} selects a target-coordinate column for profile x axes.
The default 0 automatically detects exactly one varying coordinate.

{phang}
{opt targets(numlist)} restricts target-specific views.  Weight-kernel and
influence views default to all targets through eight, then first/middle/last.

{phang}
{opt candidates(numlist)} selects tuning-candidate row numbers.
{opt allcandidates} selects every available candidate.  These options require
{cmd:tuningmode(grid)} and cannot be combined.

{phang}
{opt level(#)} sets confidence-bound coverage; default 95.

{phang}
{opt sharelimits(# #)} sets treated-share overlap guides; default .1 .9.

{phang}
{opt nlabel(#)} selects how many largest absolute influence contributions to
label per target; default 3.  Zero suppresses labels.

{phang}
{opt logweights} uses a logarithmic y axis for positive rank weights.

{phang}
{opt name(name)}, {opt replace}, and {opt nodraw} follow standard Stata graph
semantics.  The default graph name is {cmd:cate_dnn_}{it:plot_type}.

{phang}
{opt saving(filename)} saves the final graph.  {cmd:.gph}, {cmd:.png},
{cmd:.pdf}, {cmd:.svg}, {cmd:.tif}, {cmd:.jpg}, and {cmd:.eps} are recognized.

{marker results}{...}
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
Observation-level kernel and influence plots require unchanged original values,
estimation sample, and row order. Older saved fits without order signatures
must be refitted for these views; stored-result views still work after data
changes. {cmd:style()} applies to estimates and dashboard estimate panels;
{cmd:logweights} applies only to weight profiles.

{title:Stored results}

{synoptset 18 tabbed}{...}
{synopt:{cmd:r(type)}}plot type{p_end}
{synopt:{cmd:r(view)}}resolved view or estimate style{p_end}
{synopt:{cmd:r(graph)}}named graph{p_end}
{synopt:{cmd:r(n_plotted)}}number of targets, candidates, or dashboard panels{p_end}
{synopt:{cmd:r(n_omitted)}}selected estimate/tuning rows omitted because the estimate is undefined{p_end}
{synopt:{cmd:r(selected_targets)}}all selected estimate/tuning target IDs, including undefined targets{p_end}
{synopt:{cmd:r(n_rows)}}number of data rows used by a single-view plot; missing for dashboards{p_end}
{synopt:{cmd:r(targets)}}canonical target IDs present in the plot data, when applicable{p_end}
{synopt:{cmd:r(candidates)}}canonical candidate IDs present in the plot data; 0 denotes the fitted default{p_end}
{synopt:{cmd:r(n_panels)}}dashboard panel count{p_end}
{synopt:{cmd:r(panels)}}dashboard panel names in composition order{p_end}
{synopt:{cmd:r(panel_info)}}dashboard rows: panel_id, n_rows, n_plotted; no common plot-data matrix{p_end}
{synopt:{cmd:r(plot_data)}}selected data for estimates, tuning, diagnostics, weights, and other family-specific views{p_end}
{synopt:{cmd:r(saved)}}saved filename, if supplied{p_end}

{marker examples}{...}
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

{phang2}{cmd:. cate_dnn_plot}{p_end}
{phang2}{cmd:. cate_dnn_plot, type(estimates) style(profile)}{p_end}
{phang2}{cmd:. cate_dnn_plot, type(diagnostics)}{p_end}
{phang2}{cmd:. cate_dnn_plot, type(tuning) allcandidates}{p_end}
{phang2}{cmd:. cate_dnn_plot, type(overlap) view(mass)}{p_end}
{phang2}{cmd:. cate_dnn_plot, type(overlap) view(armess)}{p_end}
{phang2}{cmd:. cate_dnn_plot, type(weights) view(kernel) targets(1 3)}{p_end}
{phang2}{cmd:. cate_dnn_plot, type(influence) nlabel(3)}{p_end}
{phang2}{cmd:. cate_dnn_plot, type(dashboard) saving(cate_dashboard.pdf) replace}{p_end}
