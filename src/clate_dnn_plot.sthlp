{smcl}
{* *! version 0.1.5  2026-10-25}{...}
{viewerjumpto "Syntax" "clate_dnn_plot##syntax"}{...}
{viewerjumpto "Plot types" "clate_dnn_plot##types"}{...}
{viewerjumpto "Applicable options" "clate_dnn_plot##options"}{...}
{viewerjumpto "Data guards" "clate_dnn_plot##guards"}{...}
{viewerjumpto "Returned data" "clate_dnn_plot##results"}{...}
{viewerjumpto "Examples" "clate_dnn_plot##examples"}{...}

{title:Title}

{phang}{bf:clate_dnn_plot} {hline 2} Stata graphs for CLATE fits

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}{cmd:clate_dnn_plot} [{cmd:,} {opt type(plot_type)}
{opt view(view_type)} {opt style(style)} {opt xindex(#)}
{opt targets(numlist)} {opt candidates(numlist)} {opt allcandidates}
{opt level(#)} {opt sharelimits(# #)} {opt nlabel(#)}
{opt logweights} {opt name(name)} {opt replace} {opt nodraw}
{opt saving(filename)}]

{pstd}
Requires active {help clate_dnn} estimates. Default is type(estimates).
The command preserves the data and active fit. Stored views do not refit;
kernel/influence views reconstruct observation-level weights after
using the saved numeric model. See {help rdistnn} for the package overview
and {help clate_dnn_diagnostics} for diagnostic interpretation. R plot names describe related
workflows; these commands produce named Stata graphs, without ggplot layers,
custom-theme arguments or arbitrary twoway-option passthrough.

{marker types}{...}
{title:Plot types and scientific interpretation}

{phang}{cmd:type(estimates)} draws CLATE estimates and normal pointwise
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
Target plots also print physical-radius-90 and top-five absolute score-share ranges; these are not additional graph facets.{p_end}

{phang}{cmd:type(weights)} supports view(profile), the default rank-weight
plot; view(cumulative), with .5/.9 cumulative-weight guides; and view(kernel),
weight against fitted-metric distance by target. Profiles and cumulative weights
use actual block-averaged weights and can differ by target. Under exact weights,
the fitted profile marks the unaveraged rank-table edge n-s+k, not rank s;
a crossing exact-distance group can extend nonzero weights beyond that marker. Kernel data
are target-specific and use saved numeric observations in current fits. Kernel graphs have no
rank-s radius marker. Candidate comparison applies only to profiles.{p_end}

{pstd}
Exact-distance averaging can extend support beyond the nominal rank cutoff.
Final observation weights determine realized support and concentration; see
{help weights_dnn##geometry:distance ties}.

{phang}{cmd:type(overlap)} supports view(mass), the default, for instrument A=1/A=0 masses;
view(armess) for normalized within-arm ESS; and view(propensity) for weighted
propensity means and q10-q90 ranges. Mass uses lines with a coordinate axis
and at least four targets, otherwise bars; armess always uses bars. Propensity
always uses stored clipped-propensity summaries, not a density plot.
sharelimits() adds arm-1 share guides and reports outlying shares.{p_end}

{phang}{cmd:type(influence)} plots w*(phi_y-theta*phi_d) against fitted-metric distance,
labelling the largest absolute contributions with dataset observation numbers.
These are fitted score contributions, not individual causal effects or
jackknife pseudo-values. Selected ratios must be defined; these contributions are not divided by the first stage.{p_end}

{phang}{cmd:type(dashboard)} combines estimates, available overlap mass/propensity, cumulative weights, and tuning for grid fits; diagnostics replace overlap when those inputs are absent.
Dashboards use the full target batch and contain no kernel/influence panels.
style(), xindex() and level() affect applicable dashboard panels.{p_end}

{marker options}{...}
{title:Applicable options}

{pstd}
{opt name(name)}, {opt replace}, {opt nodraw}, and {opt saving(filename)}
apply to all types. The default graph name is {cmd:clate_dnn_}{it:type}.
replace permits reusing the graph name and overwriting an output file.
Other options are accepted only for the following types/views; combinations
that cannot affect the requested plot error.

{phang}{bf:estimates}: targets(), candidates()/allcandidates, style(),
xindex(), level(). Target and candidate selection may be combined.{p_end}
{phang}{bf:tuning}: targets(), candidates()/allcandidates, level().{p_end}
{phang}{bf:diagnostics}: targets() OR candidates()/allcandidates,
never both target and candidate selection.{p_end}
{phang}{bf:weights/profile}: targets(), candidates()/allcandidates and logweights.
logweights drops nonpositive weights and uses a log y axis.{p_end}
{phang}{bf:weights/cumulative}: targets(); no candidate selection.{p_end}
{phang}{bf:weights/kernel}: targets(); no candidate selection.{p_end}
{phang}{bf:overlap/mass}: targets(), xindex(), sharelimits().{p_end}
{phang}{bf:overlap/armess}: targets(); no xindex().{p_end}
{phang}{bf:overlap/propensity}: targets(), xindex().{p_end}
{phang}{bf:influence}: targets(), nlabel(); no candidate selection.{p_end}
{phang}{bf:dashboard}: style(), xindex(), level(); no target/candidate selection.{p_end}

{pstd}
view() applies only to weights and overlap (mass, armess, propensity). style() is auto, profile,
pointrange or forest and applies only to estimate panels. xindex() is an
integer in 0..e(p): zero detects exactly one varying coordinate, otherwise
uses target ID; positive j chooses coordinate j in fitted xcols() order.
level() is a confidence percentage, default c(level) (normally 95), rather
than the R proportion .95. It recomputes bounds without changing e().
sharelimits() requires two increasing numbers in [0,1], default .1 .9. nlabel() is a nonnegative integer, default 3; zero suppresses labels.

{pstd}
targets() and candidates() are positive integer numlists of original row IDs
in e(targets) and e(tuning_candidates); Stata ranges such as targets(1/3)
are allowed. Repeated plot IDs are removed. candidates() and allcandidates
are mutually exclusive and require tuningmode(grid). Explicit candidates
must have status=1; allcandidates plots only status=1 rows. In contrast,
{help clate_dnn_diagnostics} can report infeasible candidate rows. Noncandidate estimate,
diagnostic and rank-weight data use candidate_id=0; grid comparisons retain
actual IDs and is_default marks the fitted candidate.

{pstd}
Default estimate/diagnostic/overlap selections use all targets.
Kernel/influence defaults to all targets through eight, then
first/middle/last. Tuning has the separate three-target default above.
Overlap requires supplied instrument nuisances; paired signals do not contain instrument/propensity information. Labels and column names refer to the instrument, not treatment uptake.

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
Estimate, diagnostic, tuning, overlap and dashboard views use stored
results. Every view validates the complete current fitted snapshot.
Rank profiles/cumulative weights reconstruct the actual weights from cached
ranks, distances and target-specific groups. Current kernel/influence views use
e(model_*) snapshots and work after
original inputs change or are cleared, including after estimates save/use.
Fits without the current required snapshot need refitting.
See {help rdistnn_predict}.

{pstd}
Confidence limits are pointwise, even when the fitted covariance retains
cross-target dependence. nose keeps estimates with missing SEs and bounds.
CATE/NP zero SEs give point bounds. CLATE zero SEs remain in the matrix but status wald_available=0 suppresses Wald tests/bounds; undefined CLATE ratios are omitted from graph data.

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
{phang}{bf:target diagnostics}: {cmd:plot_id ess max_weight top5_weight_share physical_radius_90 top5_score_share is_default target_id candidate_id}.{p_end}
{phang}{bf:candidate diagnostics}: {cmd:plot_id ess max_weight top5_weight_share
is_default target_id candidate_id}.{p_end}
{phang}{bf:rank weights}: {cmd:rank weight cum_share candidate_id s is_default target_id}.{p_end}
{phang}{bf:kernel}: {cmd:target_id rank distance weight cum_share}.{p_end}
{phang}{bf:overlap/mass}: {cmd:X0 target_id overlap_available instrument_1_mass
instrument_0_mass instrument_1_share}.{p_end}
{phang}{bf:overlap/armess}: {cmd:X0 target_id arm_ess_available
instrument_1_ess instrument_0_ess}.{p_end}
{phang}{bf:overlap/propensity}: {cmd:X0 target_id propensity_available
propensity_mean propensity_q10 propensity_q90}.{p_end}
{phang}{bf:influence}: {cmd:target_id obs_id rank distance weight score
contribution phi_y phi_d}. score is the target-specific phi_y-theta*phi_d residual.{p_end}
{pstd}
For target diagnostics, plot_id=target_id, candidate_id=0, is_default=0;
for candidate diagnostics, plot_id=s and target_id is missing. Candidate IDs
must therefore come from candidate_id, not plot_id. Rank/kernel data retain
zero weights except for logweights profiles. Kernel/influence tables
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
{synopt:{cmd:r(n_panels)}}dashboard panel count{p_end}
{synopt:{cmd:r(panels)}}dashboard panel names in composition order{p_end}
{synopt:{cmd:r(panel_info)}}dashboard columns panel_id n_rows{p_end}

{pstd}
Use r(n_rows) for the actual single-view data-row count, and r(n_omitted)
for estimate/tuning rows removed because their estimates are missing.
Use r(n_panels) for dashboard panels; target and candidate IDs identify the
selected content without a view-dependent count.
An all-undefined estimate selection draws a blank graph, returns n_rows=0,
and has no plot_data matrix. Partial selections retain original target IDs.

{pstd}
A dashboard returns panel metadata, no common plot_data, and n_rows=missing.
Panels are estimates, available overlap_mass and overlap_propensity (or diagnostics when overlap is unavailable), weights_cumulative, and tuning for a grid fit.
Use individual types for selected target/candidate data or panel-specific options.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. clate_dnn_plot, level(90) name(clate_profile) replace nodraw}{p_end}
{phang2}{cmd:. matrix figure_data = r(plot_data)}{p_end}
{phang2}{cmd:. clate_dnn_plot, type(estimates) targets(1 3) candidates(1 3) style(forest) replace nodraw}{p_end}
{phang2}{cmd:. clate_dnn_plot, type(tuning) allcandidates replace nodraw}{p_end}
{phang2}{cmd:. clate_dnn_plot, type(weights) view(kernel) targets(1 3) replace nodraw}{p_end}
{phang2}{cmd:. clate_dnn_plot, type(overlap) view(mass) replace nodraw}{p_end}
{phang2}{cmd:. clate_dnn_plot, type(influence) targets(1 3) nlabel(0) replace nodraw}{p_end}
{phang2}{cmd:. clate_dnn_plot, type(dashboard) saving("clate-dashboard.pdf") replace nodraw}{p_end}

{title:Also see}

{psee}{help rdistnn}, {help clate_dnn}, {help clate_dnn_diagnostics}{p_end}
