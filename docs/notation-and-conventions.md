# Notation and conventions

This is the Stata counterpart of RdistNN's notation guide. It describes the
interface in package **0.1.0.9002** in the development snapshot. Use [getting started](getting-started.md)
for a first fit and [the R crosswalk](r-documentation-crosswalk.md) when
translating an R script.

## Observations, localization and targets

| Symbol | Stata input/result | Meaning |
|---|---|---|
| n | `e(N)`; standalone weights use `r(n)` | Complete observations used by the fitted route |
| p, kx | `e(p)`, `e(kx)`; weights use `r(kx)` | Number of localization covariates |
| X | `xcols(x1 x2)` for CATE/CLATE/NP | Variables defining the neighborhood |
| x0 | `x0(0 .5)` | One target in the covariate order |
| X0 | `matrix X0=(0,.5 \ 1,.5)` and `x0matrix(X0)` | m-by-p matrix; one target per row |
| m | `e(m)` | Number of fitted targets for CATE/CLATE/NP |
| k | `k()`; default 1 | Nearest observations selected within each size-s subsample |
| s | `s()` or resolved tuning rule | DNN kernel order/subsample size |
| lambda | `lambda()` in CATE/CLATE/NP | Multiplier for the paper order rule |
| w_i(x) | `weights_dnn ..., generate(w)` | One observation's normalized localization weight |
| ESS/neff | `e(ess_targets)`, `e(neff_targets)` | Inverse weight concentration, 1/sum(w_i^2), for each target |

The active Stata dataset replaces R's explicit `X`, signal and nuisance
vectors/data frames. Each variable supplies one column. CATE, CLATE and NP
accept exactly one of `x0()` and `x0matrix()`. The target matrix has p columns,
no missing entries, and one row per evaluation point. Its row numbers are
original target IDs used by diagnostic and plot selectors. Target values
must follow the order in `xcols()`; names on an arbitrary Stata matrix do
not reorder its columns.

`weights_dnn` takes its localization variables positionally and supports one
`x0()` at a time. PLQR takes `x()` and one `x0()`; its coefficient vector has
one entry per `d()` regressor, independently of p. It has no batch target or
integrated grid interface. Current command differences are documented,
rather than inferred from the common estimator names.

Localization X can be smaller than the state used to fit nuisance functions.
For example, a nuisance prediction may use U=(X,Z), while an effect is
reported conditional on X alone. Supplying a predicted column does not tell
Stata which variables or training rule produced it. Choose the conditioning
state, adjustment state and target interpretation before fitting.

## Sample and variable conventions

`if`/`in` restrict observations. Missing required data are excluded according
to the estimator's selected route. Copy `e(sample)` after fitting if a later
analysis needs the exact sample. Standalone weights are missing outside the
usable sample in a generated weight variable; they do not post `e(sample)`.

Raw propensity validation has a stronger rule for nuisance-based CATE and
CLATE: a selected raw propensity must be finite and in [0,1] before clipping.
Do not expect an invalid propensity to be silently removed because another
required variable is missing. The estimator articles specify each route's
inputs and sample checks. For valid raw propensities,
`e_used=min(max(e,clip),1-clip)` is used consistently in score formation and
applicable diagnostics. `clip()` defaults to 1e-8; this numerical operation
does not establish population overlap. Supplied-score routes use their
signals directly and do not clip them.

Exact equality in distances is resolved by estimation-observation order.
Rankings are deterministic for the same data/order. Sorting tied observations
can change which row receives a particular rank, so use the same observation
order when comparing with R or reconstructing weights.

## Names that mean different things

| Estimator | Localized object | Interpretation of posted covariance |
|---|---|---|
| NP | sum(w_i Y_i) | Covariance of localized conditional-mean estimates |
| CATE | sum(w_i rho_i) | Covariance of localized score averages |
| CLATE | sum(w_i phi_yi)/sum(w_i phi_di) | Paired frozen-score ratio linearization covariance |
| PLQR | sum(w_i psi_i(theta_eval)) | Covariance of the score moment at the supplied coefficient |

In CATE, `tcol()` is the binary treatment. In CLATE, `dcol()` is binary
treatment and `acol()` is the binary instrument: arm-specific nuisance means
condition on the instrument, and overlap plots describe the instrument.
PLQR `d()` supplies partially linear regressors, which need not be binary.
Its `g()` is the additive threshold nuisance, and `md()` contains conditional
means of the `d()` columns in that same order.

A supplied CATE `pseudooutcome` is one score signal. A CLATE paired-signal fit
also requires `denominator()` and takes a ratio of two local averages.
Individual signal ratios are a different calculation. PLQR `thetaeval()` is
a coefficient at which a score is evaluated; it is not an estimated local
root. Compatibility aliases do not change these meanings.

## Shapes and result lifetime

For CATE, CLATE and NP, the canonical shapes are:

| Result | Shape | Order |
|---|---|---|
| `e(targets)`, `e(x0_matrix)` | m by p | Original target rows, covariate columns |
| `e(estimates)`, `e(se_targets)` | 1 by m | Every original target |
| `e(estimate_vcov)` | m by m | Same target order on both axes |
| `e(ess_targets)`, `e(neff_targets)` | 1 by m | Same target order |
| `r(table)` after fit/replay | 9 by m | Statistics by original target |
| `r(summary)` after fit/replay | m by 6 | Original targets by summary field |
| `e(tuning_estimates)`, `e(tuning_se)` | g by m | Candidate rows by original targets |

The [reporting guide](diagnostics-and-reporting.md) names the summary fields
and diagnostic schemas. Confidence levels are percentages in Stata
(`level(90)`), rather than proportions such as R's `level=.90`.
`r()` results are temporary: copy them before another command replaces them.
A new fit replaces `e()`; use `estimates store`/`estimates restore` to keep
multiple fitted results.

For ordinary CATE/NP, `e(b)` contains all estimates and an available `e(V)`
is joint target covariance. CLATE keeps undefined targets in the canonical
matrices but omits them from `e(b)`, because Stata rejects missing posted
coefficients. CLATE `e(V)` is posted only when all posted coefficients have
available covariance; use the full canonical matrices for batch alignment.
`nose` produces missing canonical uncertainty, not precision zero.

The `e(m)` count for CATE/CLATE/NP is unrelated to PLQR's nuisance name m_D.
PLQR has its own moment/coefficient dimensions and explicit
`e(moment)`, `e(moment_vcov)` and `e(theta_eval)`; see [PLQR moments](plqr-moments.md).
Do not apply coefficient-interpretation defaults to its compatibility `e(b)`.

## Localization and inference

Larger s tightens localization for fixed n and k; k does not count all
observations with positive averaged DNN weight. Exact support reaches rank
n-s+k, whereas the asymptotic rule uses Beta-CDF increments. ESS, nonzero
support, rank radii and physical radii answer different questions. The
[weight/tuning guide](weights-and-tuning.md) explains those distinctions and
the rule `floor(lambda*n^((p+1)/(p+2)))`.

The fast jackknife measures sampling variability of the localized estimator
with frozen supplied nuisances and a fixed full-sample Mahalanobis metric.
It does not measure localization bias or uncertainty from unrestricted
nuisance fitting. Joint covariance is not a simultaneous confidence band.
Keep the estimator-specific [CATE/NP](cate-and-regression.md),
[CLATE](conditional-late.md) and [PLQR](plqr-moments.md) inference scopes distinct.
