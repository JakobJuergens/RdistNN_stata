# Diagnostics and reporting

After fitting `cate_dnn`, `clate_dnn`, or `np_dnn`, use the
matching `_diagnostics` command to inspect what each local estimate averages.
Use estimator replay for coefficient tables, and use diagnostics to assess
weight concentration, locality, and the inputs relevant to that estimator.
The [runnable example](../examples/diagnostics-and-plots.do) creates all of
its data and checks selection, saved tables, graph data, and data guards.

Start with [getting started](getting-started.md), [notation and
conventions](notation-and-conventions.md), or [CATE and
regression](cate-and-regression.md). Installed reference help is available
through `help rdistnn` and each command's help.

Use [reusing fits](reusing-fits.md) for graph-free extraction with
`rdistnn_results`, saved candidate fits and predictions after `estimates use`.
PLQR now shares the reporting tables while retaining moment-only inference.

## Report estimates and uncertainty

Each estimator displays a metadata header and pointwise normal confidence
intervals. Replay at a different confidence level without refitting:

```stata
cate_dnn, level(90)
matrix coefficient_table = r(table)
matrix summary_table = r(summary)
```

Use `clate_dnn, level(90)` or `np_dnn, level(90)` after the corresponding fit.
`level()` uses percentages (90 rather than R's .90) and defaults to Stata's
current `c(level)`, normally 95.
Replay changes reported bounds without changing `e()`.

| Result | Shape and interpretation |
| --- | --- |
| `r(table)` | 9 by m, columns in original target order; rows `b se z pvalue ll ul df crit eform` |
| `r(summary)` | m by 6; columns `estimate std_error z pvalue ci_lower ci_upper` |
| `r(level)` | Confidence level in percent |
| `e(estimates)`, `e(se_targets)` | 1 by m; all requested targets, including unavailable values |
| `e(estimate_vcov)` | m by m estimator-scale covariance, retaining cross-target covariance; missing entries under `nose` |
| `e(b)`, `e(V)` | Stata coefficient/covariance posting; CLATE may post only the defined target subset |

The `df` row is missing because reporting uses the normal distribution;
`crit` contains the normal critical value and `eform` is zero. Copy needed
`r()` results immediately, since later rclass commands replace them.
The returned covariance is already on the estimate scale: do not divide it
by effective sample size again. Joint covariance across targets does not make
the displayed intervals simultaneous confidence bands.

`nose` leaves standard errors, tests, and confidence limits missing. For
CATE/NP, a zero standard error gives point bounds; a nonzero estimate has
p=0 while its infinite mathematical z statistic is represented as Stata
missing. If both estimate and SE are zero, z and p are missing. CLATE retains
a zero SE but suppresses Wald tests and bounds when `wald_available=0`.
Undefined ratios stay visible in full-target tables. Read the components
and statuses before interpreting a CLATE interval.

## Inspect the fitted localization

These commands use stored tables, preserve the data and active estimates,
and perform no refit:

```stata
cate_dnn_diagnostics
cate_dnn_diagnostics, targets(3 1) detail
matrix local_table = r(localization)
matrix scientific_table = r(diagnostics)
```

`targets(3 1)` means target rows 3 and 1 in `e(targets)`, in that order.
Use `np_dnn_diagnostics` or `clate_dnn_diagnostics` for their active fits.
The default compact display is a selection of columns; `r(diagnostics)`
always contains the full selected table. `detail` changes only the display.

| Diagnostic | What it answers |
| --- | --- |
| `ess` | With normalized weights, `1/sum(w_i^2)` measures effective information from weight concentration. It is neither a raw observation count nor a complier count. |
| `max_weight`, `top1_weight_share`, `top5_weight_share`, `top10_weight_share` | How much mass is carried by the largest weights? Top shares use at most n observations. |
| `n_nonzero_weight` | Number of weights larger in absolute value than `sqrt(machine epsilon)`; tiny positive weights can be excluded numerically. |
| `radius_rank_k`, `radius_rank_s` | Distance to the k-th and s-th nearest observations under the fitted metric; these are rank radii, not support edges. |
| `weighted_mean_distance`, `weighted_median_distance` | Typical distance of observations actually receiving weight. |
| `physical_radius_90` (CATE/CLATE) | Ordinary Euclidean radius enclosing 90% of weight in the original covariate coordinates, even when ranks use Mahalanobis distance. Its units depend on covariate scaling. |
| Arm mass/share (CATE/CLATE) | Whether both relevant binary arms carry local weight. CATE uses treatment arms; CLATE uses instrument arms. |
| Within-arm ESS (CATE/CLATE) | `(sum_arm w)^2/sum_arm(w^2)`, after normalizing within each arm. An arm with zero weight mass has missing ESS. |
| `top5_score_share` (CATE/CLATE) | Top-five share of absolute weighted centered score contributions, which can concentrate much more than the raw weights. |
| Weighted outcome summaries (NP) | Mean, dispersion and quantiles of the local outcome distribution, rather than sampling uncertainty of its mean. |

For exact top-k DNN weights the mathematical positive support is
`n-s+k`, and for k=1 it is `n-s+1`. Increasing s tightens localization;
it is opposite in direction to an ordinary bandwidth. A low ESS can reveal
concentration; a high ESS can reveal averaging over a wide region. Neither
establishes Gaussian approximation, negligible target bias, credible nuisance
predictions, or causal identification. The candidate `degenerate` flag
marks `ess >= .9*n`, a nearly uniform-weight guide rather than a hypothesis
test. See [weights and tuning](weights-and-tuning.md).

For CATE the score contributions are `w_i*(rho_i-tau_hat)`. For CLATE they
are `w_i*(phi_y_i-theta_hat*phi_d_i)`. These are fitted score residuals,
not individual treatment effects or deletion pseudo-values. If every
contribution is zero, `score_contribution_available` can equal one while
`top5_score_share` is missing because its total is zero. CLATE marks score
concentration unavailable when its ratio is undefined or contributions are
nonfinite. These fields remain available with `nose` when the relevant
scores are defined.

## Availability follows the estimation route

| Active route | Weight/distance | Arm mass, arm ESS, propensity | Score concentration | Outcome distribution |
| --- | --- | --- | --- | --- |
| CATE AIPW or Robinson | Yes | Treatment and clipped treatment propensity | Yes | Not an NP diagnostic |
| CATE `pseudooutcome` | Yes | Unavailable | Yes | Not an NP diagnostic |
| CLATE supplied IV nuisances | Yes | Instrument and clipped instrument propensity | For defined finite paired scores | Not an NP diagnostic |
| CLATE paired signals | Yes | Unavailable | For defined finite paired scores | Not an NP diagnostic |
| NP | Yes | Not applicable | Not supplied | Yes |

Current Stata fits compute and store distance diagnostics automatically;
there is no R-style `return=c("W","dx")` retention option. Availability
flags distinguish absent inputs from small measured values. Raw
propensities must lie in [0,1]; propensity summaries use the values after
clipping to `[clip,1-clip]`. Clipping does not demonstrate overlap. A
precomputed score does not allow the package to recover treatment or
instrument inputs it never received. PLQR has no Stata diagnostics/plot
command; inspect its moment results as explained in [PLQR
moments](plqr-moments.md).

## Full diagnostic and saved-table schemas

For target scope, `r(localization)` has these columns, in order:

```text
target_id x0_1 ... x0_p
ess max_weight top1_weight_share top5_weight_share top10_weight_share
n_nonzero_weight distance_available radius_rank_k radius_rank_s
weighted_mean_distance weighted_median_distance
```

`r(diagnostics)` starts with that common block, followed by the estimator's
scientific columns in the following order.

CATE:

```text
overlap_available treated_mass control_mass treated_share
propensity_available weighted_mean_propensity weighted_propensity_q10
weighted_propensity_q90 physical_radius_90 arm_ess_available
treated_ess control_ess score_contribution_available top5_score_share
```

CLATE has the same order, replacing `treated_mass`, `control_mass`,
`treated_share`, `treated_ess`, and `control_ess` by `instrument_1_mass`,
`instrument_0_mass`, `instrument_1_share`, `instrument_1_ess`, and
`instrument_0_ess`. Its underlying `e(diagnostics)` retains CATE's internal
arm column names; the public diagnostics and plot commands apply the
scientific instrument names.

NP:

```text
outcome_available weighted_mean_y weighted_sd_y
weighted_y_q10 weighted_y_q50 weighted_y_q90
```

NP's weighted standard deviation is the square root of the weighted second
central moment with divisor `sum(w)`, without a sample correction. Quantiles
are the first observed value whose cumulative normalized weight reaches the
requested probability; they are not interpolated. With normalized weights,
`weighted_mean_y` equals the NP estimate.

The compact CATE/CLATE table contains `target_id ess top5_weight_share
physical_radius_90`, the two arm ESS fields, and `top5_score_share`. NP's
compact table contains `target_id ess max_weight top5_weight_share
weighted_mean_y weighted_sd_y`.

To save full selected diagnostics:

```stata
cate_dnn_diagnostics, targets(3 1) saving("cate-report.dta") replace
```

The dataset contains numeric variables named after the full table columns;
target IDs and coordinates are explicit variables. Matrix row labels do not
become variables. CATE/NP save `r(diagnostics)`. For CLATE target scope,
the file appends the eight component columns after `target_id` and all
five status columns described next; `r(diagnostics)` itself remains the
scientific diagnostic matrix. All saved variables are double precision.
`replace` permits overwriting; the caller's data are restored even if saving
fails. A saved `.dta` is a report, not a complete fit for `estimates restore`.

## Report CLATE components and statuses together

`clate_dnn_diagnostics, targets(3 1)` also returns aligned
`r(components)` and `r(status)`, including rows with undefined ratios.
`r(components)` has nine columns:

```text
target_id reduced_form first_stage numerator_se denominator_se
component_covariance component_correlation denominator_stat estimate
```

`denominator_stat` is the signed first-stage estimate divided by its SE.
Correlation requires both component SEs to be positive; the statistic
requires a positive denominator SE. Zero or unavailable SEs give missing
values. The paired covariance is estimator-scale, and the ratio uses common
weights for both signals. `e(component_vcov)` retains the full covariance
in interleaved target order `Y1,D1,Y2,D2,...`.

| `r(status)` column | Values |
| --- | --- |
| `ratio_status` | 0 computed, 1 zero denominator, 2 numerical failure |
| `negative_denominator` | 1 when the finite first stage is negative |
| `near_zero_denominator` | 1 when the nonzero first stage is small relative to the largest absolute input denominator signal at the configured tolerance |
| `covariance_status` | 0 computed, 1 not computed, 2 zero score variance, 3 numerical failure |
| `wald_available` | 1 when pointwise Wald inference is available |

A negative first stage remains a computed signed ratio. The near-zero flag
does not truncate, regularize, select or suppress a ratio. A zero denominator
preserves the components while making the ratio missing. These diagnostics
have no automatic weak-instrument cutoff and provide no weak-IV-uniform
guarantee. The Stata package does not implement R's ratio confidence sets,
jackknife pseudo-value concentration fields, or confidence-set geometry.
See [conditional LATE](conditional-late.md) for the causal and inference
conditions.

## Inspect tuning candidates

After `tuningmode(grid)`, candidate diagnostics summarize the weight profile
for each order, independently of targets:

```stata
np_dnn_diagnostics, candidates(3 1) detail
np_dnn_diagnostics, allcandidates saving("np-candidates.dta") replace
```

`candidates()` selects original rows of `e(tuning_candidates)` in the
supplied order. `allcandidates` includes every stored row, including
infeasible candidates. Both require a grid fit, are mutually exclusive,
and cannot accompany `targets()`. In contrast, plotting with
`allcandidates` uses available (`status=1`) candidates only.

Candidate `r(diagnostics)` and its saved dataset have these 13 columns:

```text
candidate_id computed_id is_default feasible reason status lambda s
ess max_weight top_n top5_weight_share degenerate
```

`computed_id` identifies duplicate resolved orders sharing computation.
`is_default` marks the fitted candidate. `status` is 1 for available and
0 for infeasible rows; `reason` is 0 for feasible, 1 for s<k, 2 for s>=n
when inference was requested, or 3 for s>n under `nose`. `lambda` is missing
for direct `sgrid()` choices. `top_n=min(5,n)` records the top-share count;
unavailable numerical fields stay missing for infeasible rows. The compact
candidate table has `candidate_id is_default feasible status lambda s ess
max_weight top5_weight_share degenerate`.

All diagnostics commands return `r(n_rows)`, `r(scope)` (`targets` or
`candidates`), `r(selected)` (original IDs in returned order), and `r(saved)`.
Only target scope returns `r(localization)` and CLATE component/status
matrices. Target IDs and row identities survive selection; use them rather
than treating returned row number as a new target number.

## Transfer to a report or figure

A practical reporting set is the coefficient table, selected full diagnostics,
and, for CLATE, components and statuses. Plot each source of concern with
[visualising fits](visualising-fits.md): weight profiles for concentration,
kernels/radii for locality, overlap for arm support, and influence for score
concentration. The [R documentation crosswalk](r-documentation-crosswalk.md)
records R features without a Stata implementation; the corresponding
[R diagnostics article](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/vignettes/diagnostics-reporting.Rmd)
is a conceptual reference, while Stata syntax and availability above describe
this implementation.
