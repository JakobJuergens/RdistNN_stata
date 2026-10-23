# Conditional LATE from supplied predictions

`clate_dnn` estimates a conditional local average treatment effect using a
binary instrument `A`, binary treatment `D`, and externally supplied nuisance
predictions. Supplied-prediction CLATE and ordinary paired-signal ratios are
supported workflows in package version **0.1.3** under the documented pointwise
inference conditions.
Start with the [documentation index](index.md) for installation
and the [runnable oracle example](../examples/clate.do) for a complete script.
Installed users can run `help clate_dnn` without this article.

The reporting covariates in `xcols()` can be a subset of the adjustment
variables used by the nuisance models. The causal target and identifying
assumptions are given in the
[paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
Conditional LATE generally differs from population CATE and residualized IV.

At each target, observations with exactly equal represented selected-metric
distance share the mean of their assigned rank weights, using the same groups
for numerator and denominator. Each deletion removes the actual observation
before averaging the remaining group's `n-1` weights. Treatment, instrument and
nuisance values do not define the groups. The paper assumes a continuous
distance distribution, under which old and new rules agree almost surely;
inference with persistent distance ties is not established.

## Inputs and the augmented signals

All prediction variables must be aligned with their observations. The package
fits no nuisance learners and cannot detect a prediction vector merged onto
the wrong records. Prepare and merge predictions using stable identifiers
before estimation.

| Stata option | Quantity evaluated at each observation's full adjustment state |
| --- | --- |
| `ecol(e)` | `Pr(A=1 \| V)`, instrument propensity |
| `my0col(my0)` | `E[Y \| A=0,V]` |
| `my1col(my1)` | `E[Y \| A=1,V]` |
| `md0col(md0)` | `E[D \| A=0,V]` |
| `md1col(md1)` | `E[D \| A=1,V]` |

The arm indices refer to the **instrument**, not treatment. Both adjusted
signals use the same propensity and DNN weights; the fit divides their local
averages. Averaging individual signal ratios would give a different statistic.
See `help clate_dnn` for the signal definition.

`D` and `A` must be coded 0/1 on the complete estimation sample. Raw `e`,
`md0`, and `md1` predictions must be finite and within `[0,1]` throughout the
requested `if`/`in` sample, including rows that might otherwise be excluded
for a missing outcome or covariate. Missing outcome, covariate, instrument,
treatment, or outcome-prediction values are excluded. Use `e(sample)` and
`e(N)` to audit the retained sample.

Only the instrument propensity is clipped. `clip()` defaults to `1e-8` and
must lie strictly between zero and `.5`. Treatment predictions are validated
but remain unchanged. Clipping does not establish instrument overlap or
bound an unbounded observed outcome. Nonfinite constructed signals stop the
fit; rescaling outcomes can be necessary for extreme inputs.

## One or several reporting profiles

Supply exactly one of `x0()` and `x0matrix()`. In a target matrix, each row is
a reporting profile and columns follow `xcols()` order.

```stata
matrix X0 = (-.5,0 \ 0,0 \ .5,0)
clate_dnn Y, xcols(x1 x2) x0matrix(X0) dcol(D) acol(A) ///
    ecol(e) my0col(my0) my1col(my1) md0col(md0) md1col(md1) ///
    training(supplied) tuningmode(grid)
matrix components = e(components)
matrix statuses = e(status)
matrix ratio_covariance = e(theta_vcov)
```

The variables in this fragment must already exist; the complete seeded data
generation is in [examples/clate.do](../examples/clate.do). Nuisances may use
additional controls absent from `xcols()` because only their evaluated
predictions enter the command.

`training(supplied)` is the default declaration. Alternatives are `oracle`,
`full_sample`, `oob`, and `cross_fit`. These labels record provenance; they
neither train a learner nor verify out-of-sample construction or stability.

With an explicit one-sided design, add `onesided` and omit `md0col()`, which
is then interpreted as zero. A supplied `md0col()` must be identically zero,
and observed treatment must satisfy `D<=A`. This observed check cannot prove
the structural restriction that treatment is impossible under instrument 0.
Augmented first-stage signals can still be negative for eligible
nonparticipants. They are neither probabilities nor compliance indicators.

## Supply the paired signals directly

An external score-construction workflow can supply two frozen, aligned
variables instead:

```stata
gen double e_used = min(max(e,1e-8),1-1e-8)
gen double phi_y = my1-my0 + A*(Y-my1)/e_used ///
    - (1-A)*(Y-my0)/(1-e_used)
gen double phi_d = md1-md0 + A*(D-md1)/e_used ///
    - (1-A)*(D-md0)/(1-e_used)
clate_dnn phi_y, pseudooutcome denominator(phi_d) ///
    xcols(x1 x2) x0matrix(X0) tuningmode(grid)
```

`pseudooutcome denominator()` cannot be combined with the nuisance-input
options or `onesided`. Complete outcome-signal, denominator-signal and
reporting-covariate rows form the estimation sample. No automatic propensity
clipping or causal validation occurs on this route. Instrument overlap
information is unavailable because instrument and propensity variables are
not supplied. This is the Stata route for generic paired affine ratios;
there is no separate Stata `affine_ratio_dnn` command.

## Localization and tuning

CLATE shares CATE/NP controls for exact/asymptotic weights, distances and
orders. Larger `s()` tightens localization. Use `tuningmode(grid)` for
sensitivity comparisons or specify `s()` directly. The grid retains a default
fit rather than selecting a preferred candidate. See
[weights and tuning](weights-and-tuning.md) for the default rule, custom grids
and feasibility constraints.

Uncertainty requires `1 <= k <= s < n`; `nose` permits `s=n` and leaves
uncertainty missing. All comparisons reuse the supplied signals and fitted
geometry. Supporting a weight option does not extend the paper's theory to it.

## Paired covariance, including covariance across targets

Covariance retains the dependence between numerator and denominator and
across reporting profiles, with supplied signals and geometry held fixed.
It does not jackknife deleted ratios or refit nuisance learners. Dividing
the numerator SE by the first stage discards denominator uncertainty.

| Stored matrix | Dimensions for `m` targets | Meaning/order |
| --- | --- | --- |
| `e(targets)`, `e(x0_matrix)` | `m` by `p` | Targets in requested row order |
| `e(components)` | `m` by 2 | `reduced_form first_stage` |
| `e(estimates)` | 1 by `m` | Full ratio table, including missing ratios |
| `e(se_targets)` | 1 by `m` | Ratio SEs; missing when unavailable |
| `e(theta_vcov)`, `e(estimate_vcov)` | `m` by `m` | Joint ratio covariance |
| `e(component_vcov)` | `2m` by `2m` | Interleaved `y1,d1,y2,d2,...` |
| `e(status)` | `m` by 5 | Numerical statuses listed below |
| `e(neff_targets)`, `e(ess_targets)` | 1 by `m` | `1/sum_i w_i(x)^2` |
| `e(diagnostics)`, `e(target_diagnostics)` | `m` by 25 | Shared localization/overlap schema |

Both covariance matrices are already on the estimator scale. Do not divide
them by ESS again. The Stata matrices include covariance across targets;
the pinned R ratio reference instead stores marginal covariance lists.
Cross-target covariance alone does not create simultaneous confidence bands.
ESS measures reporting-weight concentration and is not a complier count.

Stata cannot post missing coefficients in `e(b)`. Consequently, `e(b)`
contains only defined ratios, named `theta` for one target or `theta1`,
`theta2`, etc. for a batch; gaps preserve target numbers. `e(V)` is posted
only if covariance is available for every posted ratio. If all ratios are
undefined, full targets, components, statuses and sample remain available
without coefficients. Use the canonical full matrices for batch export.

## Undefined, signed and numerically small first stages

`e(status)` has one row for every requested target, in this column order:

| Column | Values |
| --- | --- |
| `ratio_status` | `0=ok`, `1=zero_denominator`, `2=numerical_failure` |
| `negative_denominator` | 1 for a finite negative aggregate denominator |
| `near_zero_denominator` | 1 for the numerical cancellation condition below |
| `covariance_status` | `0=ok`, `1=not_computed`, `2=zero_score_variance`, `3=numerical_failure` |
| `wald_available` | 1 only for an available positive-variance ratio calculation |

An exactly zero first stage gives a missing ratio while retaining the
components and their available covariance. A finite negative first stage
gives its signed ratio and is flagged, without replacement or trimming.
Numerical failures remain explicit. Ratio and covariance status must be read
together: a component covariance can be computed even when no ratio exists.
Finite component covariance can also remain available when ratio inference
is numerically unavailable.

`near_zero_denominator` flags numerical cancellation. `denominatortol()`
controls it; the flag does not change the estimate or Wald mask and is not a
statistical weak-instrument threshold.

Exact zero paired-score variance is a computed zero SE with
`covariance_status=2` and no Wald statistic, p-value, or confidence bounds.
`nose` instead gives `covariance_status=1` and missing uncertainty. Neither
case is evidence that the causal effect is known without uncertainty.

## Report, diagnose and plot

Fit and replay return a 9-by-`m` `r(table)` with rows
`b se z pvalue ll ul df crit eform` and an `m`-by-6 `r(summary)` with columns
`estimate std_error z pvalue ci_lower ci_upper`. Copy these temporary
r-class results before running diagnostics or another r-class command.
`level()` uses percentage units and defaults to Stata's `c(level)`, normally
95. `clate_dnn, level(90)` recomputes pointwise normal bounds without changing
the fitted arrays. Every target retains its place, including undefined ratios.

```stata
clate_dnn, level(90)
matrix report = r(summary)
clate_dnn_diagnostics, targets(1 3) detail
matrix component_report = r(components)
matrix status_report = r(status)
clate_dnn_plot, type(estimates) style(forest) replace
clate_dnn_plot, type(overlap) view(mass) replace
clate_dnn_plot, type(influence) targets(2) replace
clate_dnn_plot, type(tuning) allcandidates replace
```

Target diagnostics return component estimates, their SEs/covariance and
correlation, and the signed first-stage statistic when its SE is positive.
Instrument mass, within-arm ESS and propensity summaries describe `A`,
rather than treatment uptake `D`. The underlying 25-column `e(diagnostics)`
retains shared CATE slot names; `clate_dnn_diagnostics` labels these as
instrument arms. Influence shows `w*(phi_Y-theta*phi_D)`. These are weighted
paired-score contributions, not individual effects or jackknife pseudo-values.

Grid fits additionally retain `G`-by-`m` `e(tuning_estimates)`,
`e(tuning_se)`, `e(tuning_reduced_form)` and `e(tuning_first_stage)`.
`e(tuning_target_status)` has `G*m` rows ordered by candidate then target,
with `candidate_id target_id` followed by the five status columns.
Skipped candidates retain missing component/status entries. Candidate
diagnostic selection uses the shared grid tables; first-stage component and
status reports come from target selection for the default fit.

Stored estimate, overlap, diagnostics and tuning views use retained tables.
Current physical-kernel and influence views use the saved numeric sample and
signals after active data changes. `rdistnn_predict` relocalizes saved frozen signals with optional paired
ratio covariance; `rdistnn_candidate` creates a standalone grid fit. Both
preserve active estimates by default. See [reusing fits](reusing-fits.md),
[diagnostics and reporting](diagnostics-and-reporting.md) and
[visualising fits](visualising-fits.md) for selection, saving and exports.

## Inference conditions and limits

For the causal target, IV assumptions and conditions on nuisance learning and
localization, see the [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
`clate_dnn` accepts supplied predictions and holds their signals fixed for
uncertainty calculations; it does not fit or refit nuisance learners.

Reporting is pointwise Wald inference. Inspect the local first stage alongside
the ratio: a finite interval or numerical status is not a weak-instrument
coverage guarantee. R additionally offers experimental
[ratio confidence sets](https://JakobJuergens.github.io/RdistNN/articles/conditional-late-inference.html);
Stata does not implement that inversion. The
[crosswalk](r-documentation-crosswalk.md) records the remaining R-only features.
