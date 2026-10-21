{smcl}
{* *! version 0.2.1  2026-10-17}{...}
{viewerjumpto "Syntax" "np_dnn_plot##syntax"}{...}
{viewerjumpto "Plot types" "np_dnn_plot##types"}{...}
{viewerjumpto "Applicable options" "np_dnn_plot##options"}{...}
{viewerjumpto "Data guards" "np_dnn_plot##guards"}{...}
{viewerjumpto "Returned data" "np_dnn_plot##results"}{...}
{viewerjumpto "Examples" "np_dnn_plot##examples"}{...}

{title:Title}

{phang}{bf:np_dnn_plot} {hline 2} Stata graphs for nonparametric regression fits

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}{cmd:np_dnn_plot} [{cmd:,} {opt type(plot_type)}
{opt view(view_type)} {opt style(style)} {opt xindex(#)}
{opt targets(numlist)} {opt candidates(numlist)} {opt allcandidates}
{opt level(#)}
{opt logweights} {opt name(name)} {opt replace} {opt nodraw}
{opt saving(filename)}]

{pstd}
Requires active {help np_dnn} estimates. Default is type(estimates).
The command preserves the data and active fit. Stored views do not refit;
kernel views reconstruct observation-level weights after
using the saved numeric model. See {help rdistnn} for the package overview
and {help np_dnn_diagnostics} for diagnostic interpretation. R plot names describe related
workflows; these commands produce named Stata graphs, without ggplot layers,
custom-theme arguments or arbitrary twoway-option passthrough.

{marker types}{...}
{title:Plot types and scientific interpretation}

{phang}{cmd:type(estimates)} draws conditional-mean estimates and normal pointwise
confidence limits with a zero reference. style(auto) chooses a profile for
at least four selected targets with a resolved coordinate axis, a forest for
2..12 targets without such an axis, and pointranges otherwise. Explicit
style(profile) requires exactly one varying target coordinate or xindex().
An explicit xindex() can trace a supplied path even when other coordinates vary;
it does not adjust or hold those other coordinates fixed. Forest rows use
original numeric target IDs; coordinate axes use fitted covariate names.{p_end}

{phang}{cmd:type(tuning)} compares stored candidate estimates against s in
separate target panels, marks the fitted order and uses retained SEs.
Requires tuningmode(grid); defaults to every status=1 candidate and all targets
through three, then first/middle/last. This sensitivity plot does not select s
or provide inference following data-dependent candidate selection.{p_end}

{phang}{cmd:type(diagnostics)} combines ESS, maximum weight and top-five
weight-share panels. Target mode uses target IDs; candidate mode uses s.
The .9*n ESS guide marks nearly uniform weights, not a validity test.
{p_end}

{phang}{cmd:type(weights)} supports view(profile), the default rank-weight
plot; view(cumulative), with .5/.9 cumulative-weight guides; and view(kernel),
weight against fitted-metric distance by target. Profiles/cumulative weights
are common to all targets. Under exact weights, the fitted profile marks
support edge n-s+k, not rank s; larger s tightens localization. Kernel data
are target-specific and use saved numeric observations in current fits. Kernel graphs have no
rank-s radius marker. Candidate comparison applies only to profiles.{p_end}

{phang}{cmd:type(outcome)} shows the local weighted outcome mean and weighted
10th-90th percentiles. These are outcome-distribution ranges, not confidence
intervals for the mean or fitted quantile regressions. They remain available
with nose; level() does not apply. Uses a ribbon with at least four targets
and a resolved coordinate, otherwise pointranges.{p_end}

{phang}{cmd:type(dashboard)} combines estimates, diagnostics, outcome, cumulative weights, and tuning for grid fits.
Dashboards use the full target batch and contain no kernel/influence panels.
style(), xindex() and level() affect applicable dashboard panels.{p_end}

{marker options}{...}
{title:Applicable options}

{pstd}
{opt name(name)}, {opt replace}, {opt nodraw}, and {opt saving(filename)}
apply to all types. The default graph name is {cmd:np_dnn_}{it:type}.
replace permits reusing the graph name and overwriting an output file.
Other options are accepted only for the following types/views; combinations
that cannot affect the requested plot error.

{phang}{bf:estimates}: targets(), candidates()/allcandidates, style(),
xindex(), level(). Target and candidate selection may be combined.{p_end}
{phang}{bf:tuning}: targets(), candidates()/allcandidates, level().{p_end}
{phang}{bf:diagnostics}: targets() OR candidates()/allcandidates,
never both target and candidate selection.{p_end}
{phang}{bf:weights/profile}: candidates()/allcandidates and logweights;
no targets(). logweights drops nonpositive weights and uses a log y axis.{p_end}
{phang}{bf:weights/cumulative}: no target/candidate or other view-specific options.{p_end}
{phang}{bf:weights/kernel}: targets(); no candidate selection.{p_end}
{phang}{bf:outcome}: targets(), xindex(). No style() or level().{p_end}
{phang}{bf:dashboard}: style(), xindex(), level(); no target/candidate selection.{p_end}

{pstd}
view() applies only to weights. style() is auto, profile,
pointrange or forest and applies only to estimate panels. xindex() is an
integer in 0..e(p): zero detects exactly one varying coordinate, otherwise
uses target ID; positive j chooses coordinate j in fitted xcols() order.
level() is a confidence percentage, default c(level) (normally 95), rather
than the R proportion .95. It recomputes bounds without changing e().


{pstd}
targets() and candidates() are positive integer numlists of original row IDs
in e(targets) and e(tuning_candidates); Stata ranges such as targets(1/3)
are allowed. Repeated plot IDs are removed. candidates() and allcandidates
are mutually exclusive and require tuningmode(grid). Explicit candidates
must have status=1; allcandidates plots only status=1 rows. In contrast,
{help np_dnn_diagnostics} can report infeasible candidate rows. Noncandidate estimate,
diagnostic and rank-weight data use candidate_id=0; grid comparisons retain
actual IDs and is_default marks the fitted candidate.

{pstd}
Default estimate/diagnostic/outcome selections use all targets.
Kernel defaults to all targets through eight, then
first/middle/last. Tuning has the separate three-target default above.


{pstd}
saving() recognizes .gph, .png, .pdf, .svg, .tif, .jpg, and .eps.
.gph uses graph save; the others use graph export after displaying the named
graph, even when its construction used nodraw. No manual GUI action is needed.
Missing/unsupported suffixes error. To customize a figure beyond these options,
copy r(plot_data) immediately, use svmat inside preserve/restore, and draw
with standard Stata graph commands.

{marker guards}{...}
{title:Original-data guards and uncertainty availability}

{pstd}
Estimate, diagnostic, tuning, outcome and dashboard views use stored
results. Rank profiles/cumulative weights use stored n/k/s/method. They
remain available after original data change, provided their stored tables
exist. Current kernel views use e(model_*) snapshots and work after original
inputs change or are cleared, including after estimates save/use. Fits without
the current required snapshot need refitting. See {help rdistnn_predict}.

{pstd}
Confidence limits are pointwise, even when the fitted covariance retains
cross-target dependence. nose keeps estimates with missing SEs and bounds.
A zero SE gives point bounds. Weight concentration or favorable overlap does not establish Gaussian approximation, negligible target bias or validity of nuisance reuse.

{marker results}{...}
{title:Returned data and metadata}

{pstd}
For single views, r(plot_data) contains selected numeric data, and r(n_rows)
counts its rows. Copy r() results before another rclass command. Let X0 in
the schemas below mean x0_1 ... x0_p in fitted covariate order.

{phang}{bf:estimates}: {cmd:X0 candidate_id s is_default target_id estimate
se ci_lower ci_upper}.{p_end}
{phang}{bf:tuning}: {cmd:candidate_id s lambda is_default target_id estimate
se ci_lower ci_upper}.{p_end}
{phang}{bf:target diagnostics}: {cmd:plot_id ess max_weight top5_weight_share is_default target_id candidate_id}.{p_end}
{phang}{bf:candidate diagnostics}: {cmd:plot_id ess max_weight top5_weight_share
is_default target_id candidate_id}.{p_end}
{phang}{bf:rank weights}: {cmd:rank weight cum_share candidate_id s is_default}.{p_end}
{phang}{bf:kernel}: {cmd:target_id rank distance weight cum_share}.{p_end}
{phang}{bf:outcome}: {cmd:X0 target_id weighted_mean_y weighted_sd_y
weighted_y_q10 weighted_y_q50 weighted_y_q90}.{p_end}
{pstd}
For target diagnostics, plot_id=target_id, candidate_id=0, is_default=0;
for candidate diagnostics, plot_id=s and target_id is missing. Candidate IDs
must therefore come from candidate_id, not plot_id. Rank/kernel data retain
zero weights except for logweights profiles. Kernel tables
are grouped by selected target and retain original estimation-observation
order within each target, carrying rank explicitly. Estimate target filtering
follows stored target order; diagnostic selection and tuning target loops
follow requested order. Join using explicit IDs rather than row number.

{synoptset 24}{...}
{synopt:{cmd:r(type)}}requested type{p_end}
{synopt:{cmd:r(view)}}resolved view/style; diagnostics gives targets/candidates{p_end}
{synopt:{cmd:r(graph)}}graph name{p_end}
{synopt:{cmd:r(saved)}}saved filename, if supplied{p_end}
{synopt:{cmd:r(plot_data)}}single-view data matrix; absent for empty estimates/tuning or dashboards{p_end}
{synopt:{cmd:r(n_rows)}}single-view data-row count; missing for dashboards{p_end}
{synopt:{cmd:r(targets)}}sorted unique IDs actually in data, when applicable{p_end}
{synopt:{cmd:r(candidates)}}sorted unique candidate IDs in data, when applicable{p_end}
{synopt:{cmd:r(n_omitted)}}estimate/tuning rows omitted for missing estimates{p_end}
{synopt:{cmd:r(selected_targets)}}estimate/tuning target IDs, including unavailable estimates{p_end}
{synopt:{cmd:r(n_plotted)}}historical count described below{p_end}
{synopt:{cmd:r(n_panels)}}dashboard panel count{p_end}
{synopt:{cmd:r(panels)}}dashboard panel names in composition order{p_end}
{synopt:{cmd:r(panel_info)}}dashboard columns panel_id n_rows n_plotted{p_end}

{pstd}
r(n_plotted) counts estimate rows before missing estimates are dropped;
candidates for tuning/rank weights; selected rows for diagnostics;
outcome targets; targets for kernel; or dashboard panels.
It can differ from r(n_rows): G candidates and T targets yield up to G*T
estimate/tuning rows; rank weights yield n*G and kernels yield n*T.
An all-undefined estimate selection draws a blank graph, returns n_rows=0,
and has no plot_data matrix. Partial selections retain original target IDs.

{pstd}
A dashboard returns panel metadata, no common plot_data, and n_rows=missing.
Panels are estimates diagnostics outcome weights_cumulative, then tuning for a grid fit.
Use individual types for selected target/candidate data or panel-specific options.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. np_dnn_plot, level(90) name(np_profile) replace nodraw}{p_end}
{phang2}{cmd:. matrix figure_data = r(plot_data)}{p_end}
{phang2}{cmd:. np_dnn_plot, type(estimates) targets(1 3) candidates(1 3) style(forest) replace nodraw}{p_end}
{phang2}{cmd:. np_dnn_plot, type(tuning) allcandidates replace nodraw}{p_end}
{phang2}{cmd:. np_dnn_plot, type(weights) view(kernel) targets(1 3) replace nodraw}{p_end}
{phang2}{cmd:. np_dnn_plot, type(outcome) targets(1 3) replace nodraw}{p_end}
{phang2}{cmd:. np_dnn_plot, type(dashboard) saving("np-dashboard.pdf") replace nodraw}{p_end}

{title:Also see}

{psee}{help rdistnn}, {help np_dnn}, {help np_dnn_diagnostics}{p_end}
