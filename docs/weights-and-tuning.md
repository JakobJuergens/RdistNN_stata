# DNN weights and kernel-order tuning

The same localization geometry underlies Stata's CATE, CLATE, NP and PLQR
commands. This guide combines the applicable material from RdistNN's DNN
weight and kernel-order articles. [Notation](notation-and-conventions.md)
defines sample and target shapes; [weights-and-tuning.do](../examples/weights-and-tuning.do)
is a standalone runnable example.

## What k and s control

Imagine averaging a nearest-neighbor estimator over every size-s subset of
n observations. Within a subset, average its k nearest observations to x0.
The resulting complete average is a weighted average of the full sample.
`k(1)` is the default and the setting written out in the companion theory;
fixed larger k is implemented, but that implementation support alone is not
a new theorem guarantee for every configuration.

The order s has the opposite direction to a usual bandwidth. For fixed n
and k, increasing s concentrates weight nearer x0; decreasing s spreads it
out. Exact positive support reaches global rank n-s+k. In particular,
`k(1)` gives n-s+1 positive ranks, not s positive ranks. At k=s the exact
rule is uniform over n observations. The constraint k<=s<=n applies to
weight construction; delete-one uncertainty additionally needs s<n.

## Exact and asymptotic rules

For global rank r (nearest is 1), the exact finite-population coefficient is

\[
w_r=\frac{1}{k\binom ns}\sum_{\ell=1}^{k}
\binom{r-1}{\ell-1}\binom{n-r}{s-\ell},
\]

where impossible binomial choices contribute zero. The code evaluates the
combination terms on a log scale and normalizes for floating-point safety.
This is the default `method(exact)`; it does not enumerate all subsets.

`method(asymp)` replaces the without-replacement order probabilities by
Beta order-statistic increments. With F_l the Beta(l,s-l+1) CDF, it uses
`sum_l [F_l(r/n)-F_l((r-1)/n)]/k`, then normalizes. This is generally a
full-rank approximation; floating-point arithmetic may turn very small tail
weights into zero. Do not apply the exact finite-population support cutoff
to the asymptotic rule. Changing methods changes weights and estimates;
compare the methods explicitly when the approximation matters.

Euclidean `distance(l2)` is the default. Its ranking is sensitive to the
units of the coordinates; standardization can be a substantive choice.
`distance(mahalanobis)` uses the full estimation sample covariance and a
generalized inverse for rank-deficient covariance. The metric is fixed
through deletion covariance calculations. Both use deterministic distance,
then observation-order tie breaking. Changing a metric changes the
neighborhood; it is not a cosmetic plot option.

## Generate weights at one target

```stata
weights_dnn x1 x2, x0(0 0) k(1) s(8) generate(w)
scalar fitted_ess = r(neff)
summarize w
```

Run this on numeric covariates in the active dataset. `if`/`in` and complete
covariates determine the weight sample. `generate()` creates a new double
variable in original observation order, with missing values outside that
sample. Use a new variable name. Without `generate()`, the command reports
ESS but does not return a weight vector/matrix. It is r-class and posts no
estimation results; save `r(neff)`, `r(n)`, `r(k)`, `r(s)` and `r(kx)` before
another r-class command overwrites them.

`weights_dnn` supports a single `x0()` and no grid options. It requires at
least two covariates for Mahalanobis, PLQR also requires at least two covariates for Mahalanobis. CATE, CLATE and
NP permit a one-coordinate metric. All require at least two usable
observations for that metric. Use a loop with separate generated variables
for multiple standalone weight targets. Unlike R's `weights_matrix()` and
weight-result extras, Stata does not expose a public batch W/rank/order/
distance return object.

To reuse weights for a fitted signal, match its complete-case sample as well
as X, x0, method, distance, k and s. Computing weights over observations that
the estimator dropped can change n, ranks and the metric. Copy `e(sample)`
and fitted controls immediately after estimation; the runnable example
checks that a matched weighted mean equals the NP estimate.

## Default orders

CATE, CLATE and NP use

\[
s=\left\lfloor\lambda n^{(p+1)/(p+2)}\right\rfloor,
\qquad \lambda=0.4,
\]

unless an explicit order or grid is supplied. Here n is the usable sample
size and p is the number of localization covariates. For these estimators,
`lambda()` changes the multiplier. A positive `s()` fixes the order instead;
`e(tuning_source)` distinguishes `s` from a paper-resolved order. Read the
resolved `e(s)` rather than treating `e(lambda)` as a multiplier used in an
explicit-order fit. Small samples or large k can make the default infeasible;
supply a feasible order or change k. Orders are checked rather than silently
clamped into the allowed interval.

Standalone weights and PLQR also resolve `s(0)` with multiplier .4, but have
no `lambda()` or integrated grid options. For them supply `s()` explicitly
to use a different multiplier or a chosen order. PLQR always computes moment
covariance and requires s<n. `weights_dnn` can use s=n because it performs
no jackknife. `nose` allows s=n in CATE, CLATE and NP estimate-only fits.

## Integrated candidate grids

For CATE, CLATE and NP, `tuningmode(grid)` evaluates sensitivity candidates
while retaining one default fit:

```stata
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix list e(tuning_candidates)
np_dnn_diagnostics, allcandidates
np_dnn_plot, type(tuning) allcandidates name(order_sensitivity) replace
```

The default multiplier grid is .3, .4, .5 and .75 with default multiplier .4.
Use `lambdagrid()` for alternative positive multipliers, or `sgrid()` for
positive integer orders. They cannot be combined. Direct-order `sdefault()`
chooses a value from sgrid, not a candidate row number; without it the first
listed order is default. If the requested default order is absent it is
inserted. A multiplier default may likewise be inserted when the supplied
grid has no matching resolved default. This means returned row numbers should
be read from metadata instead of assumed from input positions.

```stata
np_dnn Y, xcols(x1 x2) x0matrix(X0) ///
    tuningmode(grid) sgrid(4 8 12) sdefault(8)
```

`lambda()` controls the paper-default fit in a multiplier grid. `s()` cannot
be combined with grid mode. `lambdagrid()`, `sgrid()` and `sdefault()` require
grid mode, and `sdefault()` additionally requires `sgrid()`. There is no R-style
tuning-spec registry, arbitrary tuning function, cross-validation optimizer,
or automatic winner selection. A sensitivity graph compares the specified
orders; it does not select the statistically best one or account for a
subsequent selection decision in its confidence intervals.

### Candidate metadata and availability

Each row of `e(tuning_candidates)` has eight columns:

| Column | Meaning |
|---|---|
| `candidate_id` | Original returned candidate row ID |
| `computed_id` | Shared computation group for identical resolved orders |
| `is_default` | 1 for the fit posted as the default result |
| `feasible` | Whether the order meets the sample/k/jackknife constraints |
| `reason` | 0=ok, 1=k>s, 2=jackknife requires s<n, 3=s>n on the no-SE path |
| `status` | 0=infeasible, 1=available |
| `lambda` | Multiplier for a multiplier candidate; missing for direct orders |
| `s` | Resolved integer order |

Duplicate resolved orders stay visible as distinct candidate rows, sharing
computation. Infeasible nondefault rows stay visible with missing estimates
and SEs. An infeasible default prevents fitting; the command does not silently
choose another default. Available CLATE candidate computations may still have
undefined target-specific ratios or unavailable Wald bounds: inspect
`e(tuning_target_status)` separately from candidate feasibility.

`e(tuning_estimates)` and `e(tuning_se)` are candidate-by-target matrices;
`e(tuning_diagnostics)` holds candidate weight summaries. Candidate diagnostics
can display every row, including infeasible rows. Candidate plots select
available candidates (`status=1`); explicitly selecting an infeasible row
returns an error. `allcandidates` in a plot means all available rows, whereas
`allcandidates` in the diagnostic table retains the complete candidate list.
The [reporting](diagnostics-and-reporting.md) and [plotting](visualising-fits.md)
guides show how to preserve original IDs and combine target/candidate selections.

## Jackknife and the localized target

For a scalar local average, the rank-shift implementation computes the
n delete-one estimates efficiently and uses the delete-one covariance scale
`(n-1)/n * sum_i (L_i-mean(L))*(L_i-mean(L))'`. Batch CATE/NP covariance uses
all targets together. CLATE uses paired numerator/denominator deletions and
the ratio linearization, rather than independent component variances.
The full-sample Mahalanobis metric and supplied scores/nuisances remain fixed.

The delete-one calculation cannot retain an order s=n after removing one
observation; this explains the s<n inference constraint. It does not imply
that s=n is invalid for standalone weight construction or estimate-only fits.
The posted SE measures variance around a localized target. It does not remove
localization bias, refit supplied nuisances, or turn a pointwise interval into
a simultaneous band. Orthogonality also does not establish inference with
arbitrary same-sample nuisance fits. See the estimator-specific conditions in
[CATE and regression](cate-and-regression.md) and [CLATE](conditional-late.md).

ESS is `1/sum(w_i^2)`, not positive support and not a coverage diagnostic.
A rank-s radius, weighted distance, 90-percent-mass physical radius and
within-arm ESS likewise describe different features. Inspect them alongside
candidate sensitivity rather than treating one threshold as a theorem check.

## Native acceleration and R boundaries

The optional platform-specific plugin accelerates standalone L2
`weights_dnn` calls. The estimator wrappers use their Mata routines; plugin
availability is not evidence that estimator fitting was accelerated.
By default a usable L2 plugin is attempted, with Mata fallback on load or
execution failure. `global rdistnn_use_plugin 0` forces Mata; `1` requires a
working plugin on eligible L2 calls. The global does not make Mahalanobis a
plugin-supported metric. Restore default dispatch with
`macro drop rdistnn_use_plugin`. See the [native guide](../src/plugin/README.md)
for platform assets, source builds and parity checks.

R additionally has custom distance/rank-weight/tuning registries, two-scale
TDNN weights, generic weight matrices and exposed jackknife primitives.
These have no Stata public command equivalents. The
[R documentation crosswalk](r-documentation-crosswalk.md) routes advanced
readers to the R sources without implying those APIs exist in Stata.
