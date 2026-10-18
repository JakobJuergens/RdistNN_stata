# Conditional LATE from supplied predictions

`clate_dnn` estimates a conditional local average treatment effect using a
binary instrument `A`, binary treatment `D`, and externally supplied nuisance
predictions. It is experimental functionality on `feature/clate`, outside
release 0.1.0. Start with the [documentation index](index.md) for installation
and the [runnable oracle example](../examples/clate.do) for a complete script.
Installed users can run `help clate_dnn` without this article.

The reporting covariates `X` can be a subset of the pre-instrument adjustment
state `V=(X,Z)`. Nuisance predictions must adjust for all variables needed for
instrument exogeneity; DNN weights use only `xcols()`. At a reporting profile
`X=x`, the numerator and denominator average their respective adjusted
instrument-arm contrasts over the same conditional distribution of `Z`.
Their ratio is a compliance-weighted conditional effect. It generally differs
from population CATE, an unweighted average of stratum-specific LATEs, or a
residualized-IV estimand with instrument-variance weights.

Causal interpretation requires consistency, exclusion, conditional instrument
independence given `V`, monotonicity, instrument overlap, and positive local
compliance. Without these assumptions, the calculation is a standardized
observed-data Wald functional wherever its denominator is nonzero. Supplying
predictions does not verify identification.

## Inputs and the augmented signals

All prediction variables must be aligned with their observations. The package
fits no nuisance learners and cannot detect a prediction vector merged onto
the wrong records. Prepare and merge predictions using stable identifiers
before estimation.

| Stata option | Quantity evaluated at each observation's full adjustment state |
| --- | --- |
| `ecol(e)` | `Pr(A=1 | V)`, instrument propensity |
| `my0col(my0)` | `E[Y | A=0,V]` |
| `my1col(my1)` | `E[Y | A=1,V]` |
| `md0col(md0)` | `E[D | A=0,V]` |
| `md1col(md1)` | `E[D | A=1,V]` |

The arm indices refer to the **instrument**. The outcome predictions do not
add realized treatment as a conditioning variable. For each response `R=Y,D`,
the command constructs

```text
e_used = min(max(e,clip),1-clip)
phi_R  = m_R1 - m_R0 + A*(R-m_R1)/e_used
         - (1-A)*(R-m_R0)/(1-e_used)
Delta_Y(x) = sum_i w_i(x)*phi_Yi
Delta_D(x) = sum_i w_i(x)*phi_Di
theta(x)   = Delta_Y(x)/Delta_D(x)
```

Both signals use the same propensity and localization weights. This solves
the affine score equation `sum_i w_i(x)*(theta*phi_Di-phi_Yi)=0`; averaging
individual `phi_Y/phi_D` ratios gives a different statistic.

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
a reporting profile and columns follow `xcols()` order. The compatible
`x()` alias can replace `xcols()`, but they cannot be combined.

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

Defaults are `k(1)`, `method(exact)`, `distance(l2)`, and
`s=floor(.4*n^((p+1)/(p+2)))`, where `n=e(N)` and `p` is the number of
reporting covariates. `s()` is the subsample order; `k()` counts observations
averaged within each subsample. Exact weights can be positive through rank
`n-s+k`, rather than only through rank `s`. Larger `s` tightens localization.
`k=s` gives uniform exact weights, providing no target-specific localization.

Use `s()` for a deterministic order or `lambda()` to change the multiplier
in the paper rule. `method(asymp)` uses the normalized Beta-rank
approximation. `distance(mahalanobis)` uses the full estimation sample's
covariance geometry and keeps it fixed during deletions. Distances and ties
follow the [shared conventions](notation-and-conventions.md).

`tuningmode(grid)` evaluates multipliers `.3 .4 .5 .75`, retaining the
ordinary `.4` fit by default. `lambdagrid()` replaces the multiplier list;
`lambda()` still identifies the default order, inserted if absent. Direct
`sgrid()` instead specifies integer orders, with its first order as the
default unless `sdefault()` supplies another. `sdefault()` is an order, not a
candidate row number; an absent default order is inserted. Grid options
require grid mode, `lambdagrid()` and `sgrid()` cannot be combined, and
explicit `s()` cannot be combined with grid mode. Infeasible nondefault
candidates remain in the metadata; an infeasible default stops the fit.
Rounding can give several candidate rows the same resolved order.

These comparisons reuse the supplied signals and frozen geometry. They
describe sensitivity and do not automatically select an order or justify
inference after data-dependent selection. See
[weights and tuning](weights-and-tuning.md) for candidate metadata and
feasibility codes. Covariance requires `1<=k<=s<n`; `nose` permits `s=n` and
leaves uncertainty missing. Supported weighting options do not automatically
extend the exact `k=1` inference theory.

## Paired covariance, including covariance across targets

For every deletion, CLATE removes the same observation from both signal
arrays, shifts the remaining ranks, and recomputes weights with the same
`s`, `k`, method, targets and metric. Supplied predictions remain frozen.
The joint delete-one component averages give an estimator-scale covariance
for all reduced forms and first stages. Ratio covariance linearizes each
target at its full-sample ratio, retaining both component variances and their
cross-covariance. It does not jackknife deleted ratios or refit nuisances.

For one target with components `(Delta_Y,Delta_D)` and paired covariance `C`,

```text
V_theta = (C_YY + theta^2*C_DD - 2*theta*C_YD)/Delta_D^2
```

Omitting `C_YD` or dividing the numerator SE by the first stage discards
essential paired uncertainty. A four-observation calculation in the runnable
example has weights `(1/2,1/3,1/6,0)`, components `(2,1/2)`, ratio `4`, and
ratio variance `34`; dropping cross-covariance gives `50`. One deleted first
stage equals zero, while the paired component calculation remains finite.

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

For nonzero first stages, the near-zero flag checks
`abs(Delta_D)/max_i abs(phi_Di) <= denominatortol()`. The default is
`sqrt(machine epsilon)`, approximately `1.49e-8`; zero tolerance disables
this flag. It does not change estimates, covariance, or the Wald mask. This
is a numerical diagnostic, not a statistical weak-instrument threshold.

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
Physical-kernel and influence views require the original sample, values and
observation order; changing them requires refitting. Stata has no CLATE
`predict` method: call `clate_dnn` again with new reporting targets and the
same predictions, metric convention and explicit resolved `s()` to
relocalize. A new estimation command replaces the active results. See
[diagnostics and reporting](diagnostics-and-reporting.md) and
[visualising fits](visualising-fits.md) for selection, saving and exports.

## Inference conditions and limits

Orthogonality controls first-order population nuisance drift. Valid full-sample
inference also requires the affine theorem's design, oracle, nuisance accuracy,
quantitative deletion stability, envelope and exceptional-event controls,
fixed relevance, positive score variance, and point-target undersmoothing.
A `training(full_sample)` declaration does not establish these conditions.
The proof's deleted nuisance states and the software's frozen-array variance
deletions serve different purposes. Oracle examples illustrate the computation,
rather than repeated-sampling coverage for arbitrary learned predictions.

Pointwise normal Wald intervals do not have weak-instrument-uniform validity,
and the jackknife does not remove localization bias. The R companion has a
separate `ratio_confidence_set()` helper that can produce bounded intervals,
rays, all-real or empty sets; Stata does not implement this inversion. Its
algebra alone also does not establish uniform weak-IV coverage. Do not
interpret a finite Wald interval or a numerical status as such a guarantee.

The source-pinned R counterparts are [Conditional LATE](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/vignettes/conditional-late.Rmd),
[paired inference](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/vignettes/conditional-late-inference.Rmd),
and [ratio reporting](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/man/affine_ratio_result.Rd).
The [crosswalk](r-documentation-crosswalk.md) identifies other R-only APIs.
