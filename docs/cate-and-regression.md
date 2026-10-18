# CATE and nonparametric regression

`cate_dnn` and `np_dnn` share localization, tuning, batch-target and reporting
conventions. Their inputs determine what the local average means: CATE
averages an observation-level treatment-effect signal, while NP regression
averages an outcome. This article documents package **0.1.0.9000** on `master`.
See [getting started](getting-started.md) for installation and a compact tour.

## Targets and samples

The positional argument is a numeric Stata variable. `xcols()` lists the
numeric covariates used to rank observations. Supply exactly one target
specification:

```stata
np_dnn Y, xcols(x1 x2) x0(0 .25)
matrix X0 = (0,.25 \ .5,.25 \ -.5,.25)
np_dnn Y, xcols(x1 x2) x0matrix(X0)
```

Each matrix row is a target; its columns must follow `xcols()` order. The
matrix must already exist, have at least one row, exactly one column per
localization variable, and no missing entries. Column names do not match or
reorder covariates. Resolved targets are stored in `e(targets)` with rows
`target1`, `target2`, and so forth and covariate column names. Targets need
not be observed training rows. Distance and support diagnostics help assess
how much extrapolation a proposed profile entails.

`if` and `in` restrict the candidate training sample; they do not select
target rows. The estimation sample then follows the selected input route:

| Route | Variables required to be nonmissing |
| --- | --- |
| NP regression | Outcome and `xcols()` |
| Precomputed CATE score | Positional score and `xcols()` |
| AIPW CATE | Outcome, treatment, `xcols()`, propensity, both arm means |
| Robinson CATE | Outcome, treatment, `xcols()`, propensity, pooled mean |

There is an additional propensity validation step for nuisance-based CATE:
**every raw `ecol()` value selected by `if`/`in` must be nonmissing and in
`[0,1]` before complete-case exclusion for other columns**. An invalid
propensity raises an error even on a row whose outcome or covariates are
missing. Exclude such rows explicitly when that is the intended analysis
sample. Remaining treatment values must be coded `0` or `1`. Other required
missing columns exclude observations. Empty samples are rejected.

`e(N)` reports the resulting complete-case count; `e(sample)` identifies
those training rows while the matching data are active. The default order,
distance metric and diagnostics all use this sample. Different score routes
can use different rows if their required columns have different missingness;
restrict them to a common sample before comparing estimates. Keep supplied
nuisance predictions aligned by observation, including after merges or sorts.

## CATE: three signal routes

Let `D` be binary treatment, `Y` the outcome, and `U` the full state used to
adjust for treatment assignment. Localization can use a subset `X` of `U`.
Under the required causal assumptions, the target is
`tau(x) = E[Y(1)-Y(0)|X=x]`. Nuisances are evaluated externally; the package
does not fit propensity or outcome regressions.

### Full AIPW nuisance inputs

Supply treatment propensity `e(U)`, control mean `m0(U)=E[Y|D=0,U]` and
treated mean `m1(U)=E[Y|D=1,U]`:

```stata
cate_dnn Y, xcols(x1 x2) x0matrix(X0) tcol(D) ///
    ecol(ehat) m0col(m0hat) m1col(m1hat)
```

The constructed pseudo-outcome is

```text
rho_i = m1_i - m0_i
        + D_i*(Y_i-m1_i)/eclip_i
        - (1-D_i)*(Y_i-m0_i)/(1-eclip_i).
```

The returned estimate at each target is `sum_i w_i(x)*rho_i`. This route
retains treatment-mass, propensity and within-arm effective-sample-size
diagnostics in addition to the common weight and distance summaries.

### Robinson-form nuisance inputs

Supply the pooled mean `m(U)=E[Y|U]` with propensity and treatment:

```stata
cate_dnn Y, xcols(x1 x2) x0matrix(X0) tcol(D) ///
    ecol(ehat) mcol(mhat)
```

The implemented signal is

```text
rho_i = (D_i-eclip_i)*(Y_i-m_i)/(eclip_i*(1-eclip_i)).
```

The command averages this signal with DNN weights. It does not estimate a
residual least-squares coefficient or divide by the local mean of
`(D-eclip)^2`. This distinction matters when comparing an R-learner from
another package. `mcol()` is the pooled outcome mean, not an arm-specific
mean; it cannot be combined with `m0col()` or `m1col()`.

### Precomputed pseudo-outcome

When a separate pipeline constructs a CATE signal, store it as a variable and
use `pseudooutcome`:

```stata
cate_dnn rho, xcols(x1 x2) x0matrix(X0) pseudooutcome
```

No treatment or nuisance options may accompany this route. The supplied
signal is localized unchanged. Its interpretation and first-stage validity
are the caller's responsibility. Treatment-overlap and arm-ESS fields are
unavailable because the signal alone does not reveal treatment or propensity.
Weight, distance, physical-radius and score-contribution diagnostics remain
available. The runnable
[CATE and regression example](../examples/cate-and-regression.do) shows all
three routes on the same observations.

### Clipping and aliases

For both nuisance routes,
`eclip = min(max(e,clip),1-clip)`, with default `clip(1e-8)` and
`0 < clip < .5`. Raw values outside `[0,1]` are rejected. Clipped propensities
are used in the signal, stored propensity summaries and reconstructed
postestimation scores. `clip()` does not modify a precomputed signal.
Clipping controls numerical denominators; it does not establish treatment
overlap or repair insufficient local support.

Preferred names mirror the R arguments with underscores removed:

| Preferred Stata option | Accepted legacy alias | Quantity |
| --- | --- | --- |
| `xcols()` | `x()` | Localization variables |
| `tcol()` | `treat()` | Binary treatment |
| `ecol()` | `propensity()` | Treatment propensity |
| `m0col()` | `m0()` | Control mean |
| `m1col()` | `m1()` | Treated mean |
| `mcol()` | `mhat()` | Pooled mean |

Do not supply both spellings of one quantity. `pseudooutcome`, `x0matrix()`
and `tuningmode()` correspond to R's `pseudo_outcome`, `X0` and `tuning_mode`.
The accepted CATE compatibility flag `noconstant` has no effect; there is no
intercept to remove from a weighted mean.

## Nonparametric regression

```stata
np_dnn Y, xcols(x1 x2) x0matrix(X0)
```

The estimate is `mu(x)=sum_i w_i(x)*Y_i`, a localized conditional mean.
No nuisance columns are needed. `Y` may also be a caller-defined scalar
signal, provided the interpretation is appropriate. This command uses the
same target, `if`/`in`, tuning, method, distance, `level()` and `nose` controls
as CATE. Its only legacy input alias is `x()` for `xcols()`.

NP diagnostics include the weighted outcome mean, standard deviation and
10th, 50th and 90th percentiles. The standard deviation describes outcome
spread; it is not the standard error of the mean. Weighted outcome quantiles
describe that distribution rather than uncertainty for the conditional mean.
`np_dnn_plot, type(outcome)` displays these distributional summaries;
`type(estimates)` displays mean estimates and their normal intervals.

The R package's two-scale regression `np_tdnn()` is not implemented in the
Stata companion. See the [R documentation crosswalk](r-documentation-crosswalk.md)
for this and other interface differences.

## Localization, distance and tuning

Both commands default to `k(1)`, `method(exact)`, `distance(l2)` and
`tuningmode(default)`. The order is
`floor(lambda*n^((kx+1)/(kx+2)))`, with `lambda(.4)`, complete-case count `n`
and covariate dimension `kx`. Explicit positive integer `s()` overrides the
rule; `s(0)` is the default-rule sentinel. Orders must satisfy `1 <= k <= s`.
Jackknife fits require `s<n`; estimate-only fits permit `s=n`. An infeasible
default order raises an error, including in very small samples.

`k` counts neighbors inside each size-`s` subsample, rather than the number
of full-sample observations receiving weight. For exact top-`k` DNN weights,
positive support extends through rank `n-s+k`. Increasing `s` tightens
localization. `s=1` or `k=s` gives uniform weights and no localization.
The accompanying paper's stated theory uses `k=1`; larger fixed `k` is
supported by the implementation. `method(asymp)` uses beta-order-statistic
increments, rather than the finite-population exact weights. The exact
support formula does not apply to that approximation.

L2 distances depend on covariate scales. Mahalanobis distance uses the
estimation sample's covariance inverse, with an SVD Moore-Penrose fallback
using R's relative cutoff of `1e-10`. It requires at least two complete
observations and permits one localization
covariate. A rank-deficient covariance is accepted; the generalized inverse
can leave directions without distance contribution. The metric is held
fixed through jackknife deletions.

For sensitivity analysis:

```stata
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix list e(tuning_candidates)
matrix list e(tuning_estimates)
np_dnn_diagnostics, allcandidates
```

Default multipliers are `.3 .4 .5 .75`; main results stay at the `.4`
candidate. `lambdagrid()` changes the multipliers. Alternatively,
`sgrid(8 12 20) sdefault(12)` supplies direct orders. Grid controls require
`tuningmode(grid)`; `lambdagrid()` and `sgrid()` are mutually exclusive, and
explicit `s()` cannot be combined with grid mode. The package preserves
infeasible candidate rows and duplicate orders in the candidate table; it
does not select the best candidate. See [weights and tuning](weights-and-tuning.md)
for default insertion, feasibility codes and standalone weight extraction.

## What the uncertainty measures

By default, the commands calculate a delete-1 rank-shift jackknife while
holding supplied outcomes, scores and nuisance predictions fixed. L2 ranks
under deletion coincide with recomputing distances on the remaining sample.
For Mahalanobis, uncertainty is conditional on the full-sample covariance
metric; it does not re-estimate that metric on each deletion. No nuisance
learner is refit.

For CATE, causal identification, overlap, nuisance accuracy and localization
conditions are required. Orthogonality alone does not validate arbitrary
same-sample nuisance predictions. Full-sample reuse additionally requires
the paper's nuisance-rate and full-state stability conditions. Oracle,
cross-fitted/out-of-fold and full-sample predictions are accepted, but input
acceptance is not a diagnostic of their theoretical validity. NP inference
similarly needs regression and localization-bias conditions. The jackknife
does not remove target bias for either estimator.

`level(95)` controls normal pointwise intervals; replay accepts another
level. Batch `e(V)` contains the joint cross-target covariance, so it can
support a prespecified linear contrast. A collection of pointwise intervals
does not provide simultaneous coverage over all targets.

```stata
matrix B = e(b)
matrix V = e(V)
matrix L = (1,-1,0)
matrix target_difference = B*L'
matrix target_difference_variance = L*V*L'
```

This example assumes a three-target fit and compares its first two local
estimates. A substantive contrast still needs valid inference for its
targets. `nose` skips the jackknife: `e(b)` is retained, `e(V)` is omitted,
and standard errors, tests and bounds in common summary results are missing.
The all-target `e(estimate_vcov)` remains an all-missing matrix in that case.

## Retain, replay and extract results

Copy results you need before another estimation command replaces `e()`:

```stata
cate_dnn Y, xcols(x1 x2) x0matrix(X0) tcol(D) ///
    ecol(ehat) m0col(m0hat) m1col(m1hat) tuningmode(grid)
matrix summary95 = r(summary)
matrix b_cate = e(b)
matrix V_cate = e(V)
estimates store cate_fit
cate_dnn, level(90)
matrix summary90 = r(summary)
cate_dnn_diagnostics, targets(1 3) detail
matrix selected_diagnostics = r(diagnostics)
```

For either command, the main result names and shapes are:

| Result | Shape and meaning |
| --- | --- |
| `e(b)`, `e(estimates)` | 1-by-m estimates; single-target coefficient `tau` or `mu`, otherwise `tau1 ... taum` or `mu1 ... mum` |
| `e(V)`, `e(estimate_vcov)` | m-by-m joint jackknife covariance when requested |
| `e(se_targets)` | 1-by-m marginal standard errors |
| `e(targets)`, `e(x0_matrix)` | m-by-kx resolved target profiles |
| `e(ess_targets)`, `e(neff_targets)` | 1-by-m inverse weight concentration |
| `e(diagnostics)`, `e(target_diagnostics)` | One row per target; 25 CATE fields or 17 NP fields |
| `e(tuning_candidates)` | Candidate metadata, available in grid mode |
| `e(tuning_estimates)`, `e(tuning_se)` | Candidate-by-target estimates and marginal standard errors |
| `e(tuning_diagnostics)` | One row per candidate with ESS and concentration summaries |
| `r(table)` | 9-by-m report with rows `b se z pvalue ll ul df crit eform` |
| `r(summary)` | m-by-6 report with estimate, SE, z, p-value and bounds |

`ereturn list` shows sample size, order, method, score route and tuning
metadata. `e(ess)`/`e(neff)` exist for one target; the corresponding `_min`
and `_max` scalars summarize all targets. Effective sample size measures
concentration, not physical distance, nonzero support or a theorem's sample
size. Candidate tables retain marginal SEs rather than a covariance matrix
for every candidate.

`estimates restore cate_fit` makes its CATE postestimation commands available
again; use the corresponding NP commands after restoring an NP fit. Current
fits retain numeric observations/signals and geometry, so summaries, kernels,
influence views and `rdistnn_predict` work independently of the active data.
Legacy kernel/influence views retain original-input and observation-order
guards. `rdistnn_candidate` extracts a standalone candidate using saved ranks.
`rdistnn_results` exports tables without graphs. See [reusing fits](reusing-fits.md).

For saved diagnostic tables, use `saving()` on the diagnostics command;
see [diagnostics and reporting](diagnostics-and-reporting.md). For plot
selection, graph export and which views require observations, continue with
[visualising fits](visualising-fits.md). [Notation and conventions](notation-and-conventions.md)
connects these quantities to the paper. The R reference descriptions are
[CATE](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/man/cate_dnn.Rd)
and [NP regression](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/man/np_dnn.Rd);
the installed Stata helps specify the Stata syntax.
