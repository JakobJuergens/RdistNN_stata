# DNN weights and kernel-order tuning

The same localization geometry underlies Stata's CATE, CLATE, NP and PLQR
commands. This guide combines the applicable material from RdistNN's DNN
weight and kernel-order articles. [Notation](notation-and-conventions.md)
defines sample and target shapes; [weights-and-tuning.do](../examples/weights-and-tuning.do)
is a standalone runnable example.

## What k and s control

`s` is the subsample order and `k` is the neighbor count within a subsample.
The default is `k(1)`. Increasing `s` tightens localization,
opposite to the usual bandwidth convention.
The unaveraged exact rank table reaches global rank `n - s + k`, or `n - s + 1` with
`k(1)`; neither `k` nor `s` counts all full-sample observations with weight.
At `k=s`, exact weights are uniform. Orders must satisfy `k <= s <= n`;
jackknife uncertainty additionally needs `s < n`.

The support formula refers to the unaveraged rank table. Exact-distance
averaging can give observations beyond that edge nonzero weight when their
group crosses it. Calculate realized support and concentration from the final weights.

## Exact and asymptotic rules

`method(exact)` is the default finite-population DNN rule. It computes
weights without enumerating subsets. `method(asymp)` uses a different,
generally full-support approximation. Changing the method changes weights
and can change estimates; the exact support cutoff does not apply to the
approximation. The [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
gives the construction and theoretical results. Floating-point asymptotic
tails can have zero entries; full support is not a numerical guarantee.

### Numerical limits

At subnormal scales, rounded host Beta endpoints can lose a tail mass or
produce different zero versus nonzero support in R and Stata. Order-statistic
averaging and signed TDNN blending can also lose tiny masses. This remains an
unresolved limitation and can affect estimates and covariance for extreme
signals. Normalization or a native deletion engine cannot restore probability
already lost while constructing the table. Exact weights use a different
finite-population rule, so changing methods changes the calculation.
See the [numerical reference](numerical-stability.md#asymptotic-tail-limits).

## Distance geometry

Euclidean `distance(l2)` is the default. Its ranking is sensitive to the
units of the coordinates; standardization can be a substantive choice.
`distance(mahalanobis)` uses the full estimation sample covariance and a
centered-design SVD whitening factor. Retain singular values greater than
`max(n,p) * 2^-52 * largest_singular_value`; null directions contribute zero.
This avoids explicitly forming or inverting covariance. The metric is fixed
through deletion covariance calculations. High-level ranks use represented
distance then original row index; observations sharing a distance receive the
same averaged rank weight. Changing a metric changes the neighborhood and its
target-specific groups; it is not a cosmetic plot option.

### Exact distance ties

For a group at ranks `a:b`, each observation receives
`sum(w[a:b])/(b-a+1)`. Equality is checked on the actual represented output of
the selected metric at each target, with no tolerance or rounding. Distinct
covariate profiles at equal distance belong to the same group. Average assigned
rank weights over observations, including zero entries beyond nominal support;
do not evaluate a rank rule at a mean rank or average unique profiles. Signed
and unnormalized group masses are preserved. Individual weights, support and
ESS can change; a group crossing the nominal support edge can give later
members nonzero averaged weight.

Each jackknife deletion removes the actual observation before applying the
`n-1` table to the survivors. Deleting a member of ranks `a:b` leaves `b-a`
survivors at ranks `a:(b-1)`, each receiving
`sum(w_(n-1)[a:(b-1)])/(b-a)`; a deleted singleton contributes zero. Both the
table and group size change. Existing block kernels run per target and retain
common deleted observations for paired signals and cross-target covariance.
Stata freezes fitted geometry on every deletion route. R's NP/TDNN/CATE/PLQR
`loo` route can refit sample-dependent geometry; its fast route and all paired
ratio routes freeze it.

For `x=(-1,1,2,3,4)`, target zero and `k(1) s(3)`, the exact table
`(.6,.3,.1,0,0)` becomes `(.45,.45,.1,0,0)` because the first two observations
have equal distance. With outcome `(0,10,0,0,0)`, the estimate is `4.5` and the
five deleted estimates are `(7.5,0,5,5,5)`, whose mean is `4.5`.

For exact DNN, averaging corresponds to uniform orderings within distance ties
and uniform subsets of size `s`. Within a given subset, a boundary tie with `b`
included observations and `r` remaining top-`k` slots assigns each member
inclusion probability `r/b`, retaining the usual `1/k` averaging. For `k=2`, a subset containing one closer observation
and two tied boundary observations gives weights `(.5,.25,.25)`. Including all
boundary observations with equal weight defines a different rule.

If every represented distance to a target is equal and weights are normalized,
each observation receives `1/n` and each survivor `1/(n-1)`. Increasing `s`
cannot distinguish group members. Groups can differ across targets and metrics.
Moving covariate rows requires moving their aligned outcomes and nuisances too;
estimates and covariance are then invariant up to numerical precision.

The paper assumes a continuous distance distribution, under which the earlier
and current rules agree almost surely. Commands warn because inference with
persistent distance ties is not established. Averaging does not supply a
uniform mean-squared-error improvement.

A rank table alone cannot reconstruct the weights without the target-specific
groups. Bare low-level ranks/orders retain strict-order semantics unless
supplied distances accompany them. Saved fits cache ranks, distances, groups
and the selected policy; refit objects created under an earlier policy,
including version 0.1.2.

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
returned in `r(weights)`, `r(ranks)`, `r(distances)` and `r(groups)`, with
complete-sample rows and target columns. `r(groups)` contains integer IDs for
exact represented distance equality separately in each target column. `r(row_ids)` maps rows to original observation numbers.
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

Jackknife calculations hold supplied outcomes, scores and fitted geometry
fixed. Batch CATE/NP covariance uses common deletions across targets; CLATE
also retains dependence between numerator and denominator. An order `s=n`
is available for weights or estimate-only fits but cannot be retained after
a deletion, so uncertainty requires `s<n`.

Reported uncertainty does not remove localization bias or refit nuisance
learners. Intervals remain pointwise. The paper gives the inference conditions;
[CATE and regression](cate-and-regression.md) and [CLATE](conditional-late.md)
explain the returned covariance and status fields.

ESS, support size and physical radii describe different aspects of localization.
Use [diagnostics](diagnostics-and-reporting.md) to inspect them alongside tuning
sensitivity, rather than treating any one threshold as a validity test.

### Numerical centering

The engines center signals numerically before accumulation to protect small
variation around large levels. This does not redefine the score or refit
nuisances. The [numerical implementation reference](numerical-stability.md)
describes the shared arithmetic and its limits.

## Native acceleration and R boundaries

The portable Mata path is available without a plugin. Optional native
acceleration and engine controls are documented in the
[plugin guide](../src/plugin/README.md). Ordinary users can retain automatic
dispatch; performance and parity investigations can select engines explicitly.

R's custom distance, rank-weight and tuning registries and generic jackknife
primitives have no Stata public-command equivalents. The
[R crosswalk](r-documentation-crosswalk.md) records those boundaries.

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
rows by target columns in `r(weights)`, `r(ranks)`, `r(distances)` and
`r(groups)`, plus
`r(row_ids)` mapping each row to the original dataset. `generate()` is available
for one target; it leaves excluded observations missing.

For signed two-scale tables use `weights_tdnn`; see
[two-scale regression](two-scale-regression.md). Probability-weight summaries
from single-scale DNN must not be applied to signed weights.
