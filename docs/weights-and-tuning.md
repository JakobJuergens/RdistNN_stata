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
weights into zero. Interval masses select lower or direct upper tails to avoid
subtracting CDF values rounded to one; no additional tail cutoff is imposed.
Do not apply the exact finite-population support cutoff
to the asymptotic rule. Changing methods changes weights and estimates;
compare the methods explicitly when the approximation matters.

Euclidean `distance(l2)` is the default. Its ranking is sensitive to the
units of the coordinates; standardization can be a substantive choice.
`distance(mahalanobis)` uses the full estimation sample covariance and a
centered-design SVD whitening factor. Retain singular values greater than
`max(n,p) * 2^-52 * largest_singular_value`; null directions contribute zero.
This avoids explicitly forming or inverting covariance. The metric is fixed
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
sample. Use a new variable name. Numeric weights, ranks and distances are always
returned in `r(weights)`, `r(ranks)` and `r(distances)`, with complete-sample
rows and target columns. `r(row_ids)` maps rows to original observation numbers.
The command is r-class and posts no estimation results; copy its returns before
another r-class command overwrites them.

Use `x0matrix()` instead of `x0()` for multiple targets; `generate()` requires
one target. Standalone weights have no grid options. Mahalanobis requires at
least two covariates in `weights_dnn` and PLQR; CATE, CLATE and NP permit one.
All require at least two usable observations for that metric.

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

### Numerical centering

All jackknife engines, including explicit deletions, use a numerical reference
chosen separately for each target and signal from observations with nonzero
weight in at least one deletion. This numerical shift precedes the usual
centering of deletion estimates around their mean. Prefer the first active
observation in localization-rank order when subtraction is in range and exact.
A range safeguard retains the same-sign endpoint/zero policy; use zero if a
subtraction would lose a separately represented signal component. The same
reference is used in every deletion; raw restored doubles can still round
away sub-ULP variation, so covariance uses the shifted deletions.

A nonzero entry at rank r in the n-1 weight table can use original ranks r
and r+1: deleting an earlier observation shifts the latter into position r.
For exact DNN, this union ends at n-s+k (n-s+1 for k=1). Unused observations
are excluded before subtraction, so a zero-weight outlier cannot affect the
reference or cause overflow. All n deletion identities and the original
`(n-1)/n` covariance scale are retained.

Covariance uses shifted deletions; raw deletion calculations restore the common
weighted offset. Point estimates use the original signals. This preserves the
estimator and its jackknife covariance in exact arithmetic without refitting
or statistically recentering the scores. Asymptotic weights retain every
represented nonzero entry, with no additional cutoff or renormalization.

As a shared R/Stata check, take four observations in increasing distance order,
signal `(1e16+4, 1e16+8, 1e16+12, -1)`, k=1 and s=3. The last observation
always has zero weight. Subtracting the active reference gives deletion means
`(4,0,0,0)`, hence variance `(3/4)*(9+1+1+1)=9` and SE 3. The R and Stata
unit suites check this analytical result across their fast and naive paths.

## Native acceleration and R boundaries

The optional platform-specific plugin computes scalar/batched L1, L2, scaled
L1 and scaled L2 distances for weight commands, all five estimators and
saved-fit prediction. Scale fitting, stable ranks and rank-weight tables use
Mata. The plugin also accelerates frozen-score jackknife deletions. Deletion calculations
reuse one n-1 weight table with prefix/suffix sums in both C++ and Mata.

Automatic distance dispatch keeps small jobs in Mata and attempts a usable
plugin for substantial batches, with Mata fallback if the API is absent or
execution fails. The capability probe is cached after its first use. Set
`global rdistnn_distance_engine plugin` to require this API, or `mata` to
force the portable distance path. Dropping that global restores `auto`, which
obeys `rdistnn_use_plugin`: `0` skips native calls and `1` requires them.
The distance selector takes precedence over that shared control;
explicit requirements bypass the size policy described in the native guide.
`rdistnn_last_distance_engine` records the last actual distance route.
Mahalanobis factorization remains in Mata; its projected norms use the
selected L2 engine, and its frozen-rank jackknife can use the compiled engine. Restore default dispatch with
`macro drop rdistnn_use_plugin`. For jackknife parity checks, set
`global rdistnn_jackknife_engine` to `plugin`, `mata` or `naive`; drop it to
restore automatic compiled -> fast Mata -> naive dispatch. The explicit fast
selectors report failures. Default auto mode uses fast Mata when a plugin is old,
missing, incompatible or fails during execution, and retains naive deletions
as the final fallback. Native bridge matrices created through Mata can exceed
`c(matsize)`; dispatch preserves that setting. See the [native guide](../src/plugin/README.md)
for platform assets, source builds and parity checks.

The [numerical stability guide](numerical-stability.md) describes stable TDNN
coefficients, Beta interval masses and exceptional-range L2 recovery, together
with compensated weighted/deletion sums and stable diagnostic calculations.

R additionally has custom distance/rank-weight/tuning registries, generic
weight specification objects and exposed jackknife primitives.
These have no Stata public command equivalents. The
[R documentation crosswalk](r-documentation-crosswalk.md) routes advanced
readers to the R sources without implying those APIs exist in Stata.

## Distance geometry

All estimation and weight commands accept `distance(l1)`, `distance(l2)`,
`distance(l1_scaled)`, `distance(l2_scaled)` and `distance(mahalanobis)`.
Scaled metrics fit one scale per covariate on the complete estimation sample.
`scalerule(sd)` uses the sample standard deviation; `mad` uses median absolute
deviation multiplied by 1.4826; `iqr` uses R's type-7 quartiles; `range` uses
the observed range. `scalerule(custom) scalevec(...)` supplies the scales
directly in covariate order. Absolute custom scales and the `scaleeps()`
floor follow R's contract. The epsilon defaults to 1e-12.

Scales and Mahalanobis geometry stay fixed across all deletion calculations
and saved-model predictions. A new target does not refit its metric.

`weights_dnn` accepts either `x0()` or `x0matrix()`. It returns complete-sample
rows by target columns in `r(weights)`, `r(ranks)` and `r(distances)`, plus
`r(row_ids)` mapping each row to the original dataset. `generate()` is available
for one target; it leaves excluded observations missing.

For signed two-scale tables use `weights_tdnn`; see
[two-scale regression](two-scale-regression.md). Probability-weight summaries
from single-scale DNN must not be applied to signed weights.
