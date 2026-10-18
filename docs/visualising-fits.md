# Visualising fits

`cate_dnn_plot`, experimental `clate_dnn_plot`, and `np_dnn_plot` turn the
active fit into Stata graphs. Choose a panel based on the question: estimates
for pointwise uncertainty, diagnostics for concentration, weights for local
geometry, overlap for arm support, influence for score contributions, and
NP outcome panels for the local outcome distribution. The commands return
numeric graph data as well as named graphs, so a report can use the precise
rows behind a panel.

The [complete example](../examples/diagnostics-and-plots.do) generates seeded
data and runs every supported plot type without GUI interaction or downloaded
inputs. See [diagnostics and reporting](diagnostics-and-reporting.md) for the
underlying quantities, [getting started](getting-started.md) for installation,
and `help rdistnn` for the installed overview.

## Gallery from the runnable example

These plots use the seeded tutorial data in
[diagnostics-and-plots.do](../examples/diagnostics-and-plots.do). They illustrate
reporting choices, not empirical evidence or nuisance-learner performance.

![CATE profile with normal pointwise confidence ribbon](figures/cate-profile.png)

The CATE profile connects five target estimates along x1 with x2 held at zero.
Its ribbon is pointwise uncertainty, not a simultaneous confidence band.

![CATE dashboard combining estimates, localization, overlap, and rank weights](figures/cate-dashboard.png)

The dashboard places the estimate beside concentration, arm-overlap and rank
weight summaries. The panels answer different questions; none verifies the
conditions for valid inference by itself. The nested tuning facets are compact
in a dashboard; use a standalone tuning panel when the target-specific paths
need detailed reading.

![NP localized outcome means and weighted tenth to ninetieth percentiles](figures/np-outcome.png)

The NP outcome band is a weighted 10th–90th percentile range. It describes
local outcome dispersion and is different from uncertainty about the mean.

![CLATE forest plot for supplied paired tutorial signals](figures/clate-forest.png)

The paired-signal CLATE forest plot emphasizes target identities and ratio
uncertainty. These synthetic paired signals illustrate computation; they do
not by themselves identify a causal complier effect.

## Draw an estimate panel

After a CATE batch fit:

```stata
cate_dnn_plot, type(estimates) level(90) name(cate_profile) replace nodraw
matrix figure_data = r(plot_data)
```

`type(estimates)` is the default. `level()` is a percentage, defaulting to
Stata's current `c(level)`; it recomputes normal bounds using stored SEs.
NP labels estimates as conditional means and CLATE labels them as CLATE.
Their covariance is already on the estimator scale. Intervals are pointwise,
even though the fitted covariance retains cross-target dependence.

`style(auto)` uses the selected targets to choose:

| Selected targets | Style |
| --- | --- |
| At least four, with exactly one varying coordinate (or an explicit `xindex()`) | `profile`: line and confidence ribbon |
| More than one and at most twelve, without a resolved varying coordinate | `forest`: horizontal ranges with target IDs |
| Otherwise | `pointrange`: estimates with vertical ranges |

Choose `style(profile)`, `style(forest)`, or `style(pointrange)` explicitly
when that layout serves the report. `style(profile)` requires an identified
x coordinate. `xindex(j)` selects coordinate j in the `xcols()` order;
`xindex(0)` automatically detects exactly one varying coordinate and otherwise
uses target ID. Explicit `xindex()` lets the user choose an axis even when
several coordinates vary; the resulting line is a path through the supplied
profiles, not an automatically adjusted conditional response curve. Axis
labels use the covariate name; forest rows use original numeric target IDs.
There are no custom-label or arbitrary `twoway` passthrough options.

`nose` plots estimates without intervals. CATE/NP zero SEs give point bounds.
CLATE zero SEs remain in the returned data but have no Wald bounds when its
status mask rules out inference. Undefined CLATE estimates are omitted from
the graph data while their original IDs remain in full fit/diagnostic tables.
An entirely undefined selection produces a blank estimate graph,
`r(n_rows)=0`, and no `r(plot_data)` matrix. Check `r(n_omitted)` and
`r(selected_targets)` when interpreting an estimate or tuning panel.

## Compare tuning choices

A grid fit retains candidate estimates and marginal SEs. Plotting does not
choose a candidate or update the fit:

```stata
cate_dnn_plot, type(estimates) targets(1 3) candidates(1 3) ///
    style(forest) name(cate_comparison) replace nodraw
cate_dnn_plot, type(tuning) targets(1 3) allcandidates ///
    name(cate_sensitivity) replace nodraw
```

In `type(estimates)`, candidates share the same target axes; the fitted
candidate is emphasized. `type(tuning)` uses s as its x axis, marks the
fitted order and separates targets into panels. It defaults to all available
candidates and, when there are more than three targets, to first/middle/last
targets. With at most three targets it uses all of them. The built-in
multiplier grid is `.3 .4 .5 .75` with `.4` fitted by default; custom grids
are described in [weights and tuning](weights-and-tuning.md).

`targets()` and `candidates()` are integer numlists of original target and
candidate row IDs, including Stata range notation such as `targets(1/3)`.
Plot commands remove repeated IDs. `candidates()` and `allcandidates` are
mutually exclusive. Candidate plotting requires `tuningmode(grid)` and only
`status=1` candidates: explicit infeasible IDs error, whereas
`allcandidates` skips them. Diagnostics tables can retain infeasible rows.
The fitted noncandidate view uses `candidate_id=0`; grid comparisons use
actual candidate IDs, even for the fitted candidate. Inspect `is_default`
rather than assuming candidate 1 is fitted.

## Weight geometry and concentration

All three estimators share:

```stata
cate_dnn_plot, type(weights) view(profile) logweights replace nodraw
cate_dnn_plot, type(weights) view(cumulative) replace nodraw
cate_dnn_plot, type(weights) view(kernel) targets(1 3) replace nodraw
cate_dnn_plot, type(diagnostics) targets(1 3) replace nodraw
cate_dnn_plot, type(diagnostics) allcandidates replace nodraw
```

`view(profile)` plots DNN weight against rank. Under exact weights, the
fitted profile marks the support edge `n-s+k`; it is not rank s.
`logweights` removes zero weights and changes the y scale. Candidate
comparisons are available only for profiles. `view(cumulative)` shows the
cumulative rank weight with .5/.9 guides. Rank profiles and cumulative
weights are common across targets, so `targets()` is rejected there.

`view(kernel)` plots target-specific weight against distance under the fitted
L2 or Mahalanobis metric. It reconstructs weights from the original data;
the returned matrix stays in original estimation-observation order and
carries `rank` explicitly. It defaults to all targets up to eight, then
first/middle/last. Stata's kernel plot does not add the R rank-s radius marker.

`type(diagnostics)` combines ESS, maximum weight, and top-five weight share.
It uses target ID, or s for candidate scope, as its x axis. The ESS panel's
`.9*n` guide marks near-uniform weights. CATE/CLATE target views additionally
print the range of 90% physical radii and score-contribution shares; these
are not extra graph facets. See [diagnostic interpretation](diagnostics-and-reporting.md#inspect-the-fitted-localization)
before treating a large or small diagnostic as evidence of inferential validity.

## Treatment and instrument overlap

After a CATE AIPW/Robinson fit or CLATE supplied-IV-nuisance fit:

```stata
cate_dnn_plot, type(overlap) view(mass) sharelimits(.1 .9) replace nodraw
cate_dnn_plot, type(overlap) view(armess) targets(1 3) replace nodraw
cate_dnn_plot, type(overlap) view(propensity) replace nodraw
```

CATE labels arms treated/control. CLATE uses instrument=1/instrument=0 and
instrument propensity; treatment uptake is a different quantity. Replace
the command prefix with `clate_dnn_plot` for those fits.

`mass` displays local arm masses. It uses lines when a coordinate axis is
resolved and at least four targets are selected, and bars otherwise.
`sharelimits()` supplies two increasing values in [0,1], with .1/.9 as the
default guides; it also reports targets whose arm-1 share lies outside them.
`armess` always uses bars of normalized within-arm ESS, so it rejects
`xindex()`. `propensity` displays weighted means and weighted q10–q90 ranges
of clipped propensities, using stored summaries. It does not draw an
observation-level propensity density or an overlap test. Precomputed CATE
scores and paired CLATE signals have no arm/propensity data; overlap views
error with an availability message.

## Trace score contributions

```stata
cate_dnn_plot, type(influence) targets(1 3) nlabel(3) replace nodraw
clate_dnn_plot, type(influence) targets(1 3) nlabel(0) replace nodraw
```

Run each command while its matching fit is active. CATE plots
`w_i*(rho_i-tau_hat)` against distance. CLATE plots
`w_i*(phi_y_i-theta_hat*phi_d_i)` and also returns both component signals.
These are score contributions, not individual causal effects, deletion
pseudo-values, or influence divided by the first stage. CLATE requires
defined selected ratios. `nlabel()` labels that many largest absolute
contributions per target with original dataset observation numbers; its
default is 3 and zero disables labels. Default target selection matches the
kernel view. NP has no influence view.

## Distinguish an outcome distribution from a confidence interval

```stata
np_dnn_plot, type(estimates) level(90) replace nodraw
np_dnn_plot, type(outcome) targets(1/3) replace nodraw
```

`type(outcome)` shows the weighted outcome mean with weighted 10th–90th
percentiles. These ranges describe outcome heterogeneity in the local
weighted sample. They are neither confidence intervals for the conditional
mean nor conditional quantile-regression estimates. They can remain wide
when the mean is precisely estimated, and are still available with `nose`.
The panel uses a ribbon when a coordinate axis is resolved and at least
four targets are selected, otherwise pointranges. There is no `level()`
option for outcome bands: the displayed probabilities are fixed.

## Applicable options

`name()`, `replace`, `nodraw`, and `saving()` apply to every type. Other
options are accepted only where they can affect the requested graph; invalid
combinations return an error.

| Type/view | `targets()` | Candidate selection | `style()` | `xindex()` | `level()` | Additional option |
| --- | --- | --- | --- | --- | --- | --- |
| `estimates` | Yes | Yes, jointly with targets | Yes | Yes | Yes | — |
| `tuning` | Yes | Yes, jointly with targets | — | — | Yes | Requires grid fit |
| `diagnostics` | Yes | Yes, exclusive of targets | — | — | — | — |
| `weights`, `profile` | — | Yes | — | — | — | `logweights` |
| `weights`, `cumulative` | — | — | — | — | — | — |
| `weights`, `kernel` | Yes | — | — | — | — | Original data required |
| `overlap`, `mass` (CATE/CLATE) | Yes | — | — | Yes | — | `sharelimits()` |
| `overlap`, `armess` (CATE/CLATE) | Yes | — | — | — | — | — |
| `overlap`, `propensity` (CATE/CLATE) | Yes | — | — | Yes | — | — |
| `influence` (CATE/CLATE) | Yes | — | — | — | — | `nlabel()`, original data |
| `outcome` (NP) | Yes | — | — | Yes | — | — |
| `dashboard` | — | — | Yes | Yes | Yes | Full target batch |

`view()` applies only to weights and overlap; it defaults to `profile`
and `mass`, respectively. `style()` changes estimate panels, including the
estimate panel inside a dashboard. `xindex()` must be an integer in 0..p.
`nlabel()` must be a nonnegative integer. Plot `level()` uses Stata's usual
confidence-level validation.

## Preserve and export a figure

```stata
cate_dnn_plot, type(estimates) name(cate_main) replace nodraw
cate_dnn_plot, type(dashboard) name(cate_overview) replace nodraw ///
    saving("cate-overview.pdf")
```

Default graph names are `cate_dnn_TYPE`, `clate_dnn_TYPE`, or `np_dnn_TYPE`.
A repeated graph name requires `replace`; the same option permits overwriting
an output file. Supported `saving()` suffixes are `.gph`, `.png`, `.pdf`,
`.svg`, `.tif`, `.jpg`, and `.eps`. `.gph` saves an editable Stata graph;
other suffixes call `graph export`. Raster/vector export displays the named
graph before export even if creation used `nodraw`; it does not require
manual GUI interaction. Unsupported or missing suffixes error. The runnable
example keeps its optional exports under ignored `build/documentation/`.

Stata panels are named graphs, without R's ggplot object, theme, layer, or
patchwork list interfaces. To create a custom Stata figure, immediately
copy `r(plot_data)`, transfer it into a temporary dataset with `svmat`, and
use standard graph commands. Preserve/restore the caller's data around that
conversion. The graph data schemas below are Stata's interface, rather than
R's `plot_data_*()` frame schema.

## Returned graph data and counts

Let `X0` below denote `x0_1 ... x0_p`. Columns are in the order shown.

| Single-view result | `r(plot_data)` columns |
| --- | --- |
| Estimates | `X0 candidate_id s is_default target_id estimate se ci_lower ci_upper` |
| Tuning | `candidate_id s lambda is_default target_id estimate se ci_lower ci_upper` |
| Target diagnostics, CATE/CLATE | `plot_id ess max_weight top5_weight_share physical_radius_90 top5_score_share is_default target_id candidate_id` |
| Target diagnostics, NP | `plot_id ess max_weight top5_weight_share is_default target_id candidate_id` |
| Candidate diagnostics, all families | `plot_id ess max_weight top5_weight_share is_default target_id candidate_id` |
| Rank weights, profile/cumulative | `rank weight cum_share candidate_id s is_default` |
| Kernel | `target_id rank distance weight cum_share` |
| CATE overlap, mass | `X0 target_id overlap_available treated_mass control_mass treated_share` |
| CATE overlap, arm ESS | `X0 target_id arm_ess_available treated_ess control_ess` |
| Propensity overlap, CATE/CLATE | `X0 target_id propensity_available propensity_mean propensity_q10 propensity_q90` |
| CATE influence | `target_id obs_id rank distance weight score contribution` |
| CLATE influence | `target_id obs_id rank distance weight score contribution phi_y phi_d` |
| NP outcome | `X0 target_id weighted_mean_y weighted_sd_y weighted_y_q10 weighted_y_q50 weighted_y_q90` |

CLATE mass/arm-ESS columns replace `treated_mass control_mass treated_share`
and `treated_ess control_ess` with `instrument_1_mass instrument_0_mass
instrument_1_share` and `instrument_1_ess instrument_0_ess`. CLATE influence
`score` is `phi_y-theta_hat*phi_d` for that target; CATE `score` is rho.
For target diagnostics `plot_id=target_id`, `candidate_id=0`, and
`is_default=0`; for candidate diagnostics `plot_id=s` and `target_id` is
missing. Thus `plot_id` is an axis value, not a substitute for a candidate ID.

`r(n_rows)` is the number of returned data rows. With G candidates and T
targets, estimate/tuning tables can contain G*T rows; rank-weight tables
contain n*G, and kernel/influence tables contain n*T. Zero weights remain
in those matrices except that `logweights` removes them from profiles.
Estimate/tuning tables omit missing estimates; overlap tables omit
unavailable rows. Distance-based matrices are grouped by requested target
and retain original estimation-observation order within that target.

`r(n_plotted)` has a different historical meaning: selected estimate rows
before undefined estimates are dropped; number of candidates for tuning
or rank-weight views; number of target/candidate diagnostic rows; number of
overlap/outcome targets; number of kernel/influence targets; or number of
dashboard panels. Use `r(n_rows)` for tabulation rather than assuming these
counts agree. Estimate/tuning also return `r(n_omitted)` and
`r(selected_targets)`, including targets with undefined estimates.

`r(type)`, `r(view)`, `r(graph)` and `r(saved)` describe the graph.
`r(targets)` and `r(candidates)`, when applicable, identify IDs actually
present in the graph data in sorted unique order. Copy the matrix and any
metadata before another rclass command. Use explicit IDs to join tables;
estimate/overlap/outcome target filtering follows stored target order,
while diagnostic selection and tuning target loops follow requested order.

## Dashboards

```stata
np_dnn_plot, type(dashboard) name(np_overview) replace nodraw
local panels "`r(panels)'"
matrix panel_counts = r(panel_info)
```

CATE/CLATE dashboards combine estimates, available overlap mass/propensity,
cumulative weights, and tuning when the fit used a grid. When overlap inputs
are absent, diagnostics replace overlap panels. NP dashboards combine
estimates, diagnostics, outcome summaries, cumulative weights, and tuning
for grid fits. Their graphs retain CATE, CLATE, instrument or conditional-mean
labels as appropriate. Dashboards do not reconstruct kernels or influence.

`r(panels)` lists composed panels in order, `r(n_panels)` counts them, and
`r(panel_info)` has columns `panel_id n_rows n_plotted`. Its counts retain
the definitions for each panel. A dashboard returns no common
`r(plot_data)` matrix and has `r(n_rows)=.`. Select individual views when
you need their data, target restrictions, or panel-specific options.

## Original-data guards and restored fits

Stored-result estimate, diagnostics, tuning, overlap, outcome and dashboard
views remain usable after data changes. Rank profiles/cumulative weights
are recreated from stored n/k/s/method, without observation inputs. Kernel
and influence views instead validate the original estimation sample,
referenced covariates and score variables, their values, and row order.
Changed values, changed sample, missing variables, or reordered observations
require a refit; an unchanged sample count alone is insufficient.

`estimates store`/`estimates restore` can restore a fit within a session;
observation-level views additionally need its original verified data and
sample. Older saved CATE/NP fits without canonical summary fields retain
coefficient/covariance fallbacks for replay and estimates. Stored diagnostic
views require the tables needed for that view. Older fits without ordered
signatures must be refitted for kernel/influence views. A legacy object
cannot gain missing state by merely restoring it.

The [R visualization article](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/vignettes/visualising-fits.Rmd)
provides the related R workflow. R-only plot conditioning/calibration,
two-scale regression and generic GMM plotting are not Stata commands;
see the [documentation crosswalk](r-documentation-crosswalk.md) and
[notation and conventions](notation-and-conventions.md).
