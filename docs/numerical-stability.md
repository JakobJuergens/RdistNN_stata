# Numerical implementation reference

RdistNN and its Stata companion share numerical contracts for distance-localized
estimation. C++, R and Mata implementations of a selected method must agree on
rank/tie order, support, target identities and computational statuses, with
justified tolerances for continuous results. Independent references establish
accuracy in addition to agreement between engines. See the
[contributor contract](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md#numerical-parity-contract)
and [native guide](../src/plugin/README.md) for verification and engine controls.

## Stable elementary calculations

TDNN coefficients use equivalent logarithm/exponential formulas that avoid
subtracting a rounded value near one when orders are close. Asymptotic DNN and
TDNN probability increments select lower or direct upper Beta tails, avoiding
subtraction of two CDF values that have both rounded to one. These choices
preserve the selected mathematical weighting rule and its normalization.
Very small weights may still underflow the supported numeric representation;
the code does not impose an additional tail cutoff.

The legacy C++ `weights_l2` API requires its incomplete-Beta continued
fraction to converge to a finite positive value within 200 iterations.
Failure raises Stata numerical error 498 before any output-variable writes,
with guidance to use `weights_dnn`. Standard estimation and weight commands
already evaluate Beta weights in Mata. This safeguard detects failed
iteration; it does not certify an error bound for every converged probability.
Exact endpoints and legitimate tail underflow retain their existing behavior.

Fitted IQR/MAD scales use R's guarded type-7 interpolation, including its
weighted-endpoint arithmetic. SD fitting uses common conservative range gates
before a quadratic sum can overflow or underflow. Exceptional columns use an
exact binary scale, shift by an observed value, and compensate the normalized
mean and squared-residual sums. This preserves variation around large common
offsets and aligns R/Mata arithmetic; ordinary columns retain their existing
SD calculations. Nonfinite fitted scales are rejected before the epsilon
floor is applied; only finite small or constant scales are clamped. Fitted
scales must themselves fit Stata's finite nonmissing range.

Scaled L1 and L2 share binary-half subtraction recovery when a raw difference
reaches `2^1023` and the inverse scale is below one. Ordinary differences retain
subtract-then-scale rounding and ordered accumulation in every engine.

L2 and fitted-scale L2 retain ordered coordinate accumulation with rounded
squares for ordinary rows. All engines recover exceptional rows with scaled
sum of squares when the accumulated squared sum reaches `2^1023`, a nonzero
coordinate's square is subnormal or zero, or a scaled raw difference reaches
`2^1023`. A huge raw difference with inverse scale below one uses binary-halved
operands before multiplying by twice that inverse scale. The shared boundary
accounts for Stata's narrower finite nonmissing range; using each host's own
overflow trigger could change neighbors. Recovery rescans exceptional rows
with constant accumulator state. Finite representable norms remain usable;
genuine out-of-range outputs retain failure/unavailable-value behavior.

Shared support-aware centering protects deletion covariance from large common
signal levels. Paired CLATE covariance also uses scaled linearized scores;
these operations retain numerator/denominator dependence. Numerical stability
does not change inference assumptions or the interpretation of a near-zero
first stage.

## Fitted Mahalanobis factorization

`_rdistnn_geometry.ado` centers each column safely on its own binary scale,
then chooses a common scale from centered variation before SVD. A huge
constant column cannot erase a small varying column. It retains singular values
greater than `max(n,p) * 2^-52 * largest_singular_value`, matching R. Projecting
rows through the resulting whitening factor and applying the existing native
or Mata L2 norm gives the sample-covariance Moore-Penrose metric without
forming a covariance inverse or squaring the design's condition number.
Numerical null directions contribute zero; a rank-zero design gives zero
distances. The factorization uses Mata; transformed norms honor the distance
engine selector and report its actual route. Projection skips exactly null
coordinates and recovers overflowing raw differences on the fitted scale.
Exceptional projection rows normalize before multiplication and restore the
binary scale ratio afterward, so a factor below one can recover a finite
projection after division alone would overflow.
Whitening projections use the shared weighted-signal reduction on separately
rounded double products. Reducing each row independently prevents scalar and
batch BLAS kernels from choosing different cancellation results; transformed
norms continue to use the selected distance engine. R applies the analogous
ordered compensated projection in its native and interpreted routes.

The explicit design-scale cutoff replaces the old inverse-attempt and
covariance-truncation policy. New fits may choose different neighbors when the
old arithmetic lost a resolvable direction. SVD libraries can still differ at
the rank boundary; universal bitwise factorization equivalence is not promised.
Saved geometry stores a header with scale and rank (rank when `p>=2`), a
raw observed origin and the zero-padded square factor: `(p+2) x p`.
Legacy `p x p` inverse snapshots and old fits without stored geometry retain
their original metric. The outer model version remains compatible. Target
reordering, subsetting and serialized prediction always reuse fitted geometry.

## Shared covariance arithmetic

Ordinary covariance entries retain vectorized cross-products. Off-diagonals
whose absolute value is at most `sqrt(eps)` times the two column norms are
recomputed with compensated double products, including the two-difference
residual of mean subtraction. This recovers small represented cross terms
lost by centering or BLAS cancellation. The residual workspace is one matrix
the size of the signals; diagonal range and normalization rules stay unchanged.

`_rdistnn_covariance.ado` owns centering, exceptional binary normalization,
cross-products, safe restoration and symmetry. Estimator adapters retain their
statistical factors, classes and failure-status policy. All native, Mata and
naive deletion engines return the same aligned deletion layout to this module;
there is no separate native covariance policy. TDNN, multi-component moments,
paired CLATE and saved-fit prediction call the same helper. Its R counterpart
is `R/utils-covariance.R`.

Ordinary columns use vectorized cross-products. Normalize before centering if
amplitude reaches `sqrt(2^1020/n/max(1,factor))` or is positive and below
`sqrt(DBL_MIN)`; small centered spreads also normalize. An observed-value shift
protects differences at large offsets. Binary exponent restoration occurs in
two steps, avoiding an overflowing or underflowing scale product. A shared
diagonal reduction keeps marginal and joint variances aligned. Symmetric
assembly preserves diagonal values and never doubles a large covariance.

Exact constant columns have zero variance. Nonzero variation with a zero,
negative or nonfinite final variance is unavailable. Ordinary estimators raise
a rescaling error instead of reporting false-zero or infinite uncertainty;
`nose` keeps point-only estimation available. CLATE requests failure masks and
preserves its existing numerical-failure status while keeping valid targets.
The final value must fit the host representation: some variances are finite
in R but exceed Stata's smaller finite nonmissing range. Numerical safeguards
do not establish the assumptions for inference.

Shared reporting keeps estimates, standard errors and critical values in
numeric scalars when constructing returned matrices. Avoid a round trip
through text macros: their formatting can round extreme numeric values even
when the original covariance and estimates are accurate. Display formatting
does not reduce the precision of returned tables or saved-fit predictions.

## Compensated weighted and deletion sums

`_rdistnn_accumulation.ado` owns weighted signal products, rank-weight mass and
retained correction arithmetic. The standalone native header
`compensated_arithmetic.hpp` is mirrored in RdistNN. Existing estimators supply
their signals, ranks, tables and statistical factors, while covariance policy
stays in `_rdistnn_covariance.ado`. Candidate fits and saved prediction reuse
these same adapters.

Sensitive weighted reductions follow Neumaier compensated addition. Each
product is rounded to double before addition; native code prevents fused
multiply-add from changing that contract. Mata uses native `quadsum()` with
explicit missing-value propagation on separately rounded product vectors for
weighted signal means; each output uses one vector, preserving valid components
when another component is missing. Shared R/native/Mata fixtures check these primitives against the
canonical arithmetic. General covariance matrix algebra remains vectorized.

The prefix engine stores the forward leading sum and its correction, then
streams the reverse state. It combines prefix leading value, suffix leading
value, prefix correction, suffix correction and restored offset before
rounding the result. R/Mata reverse sweeps use blocks of at most 2,048 rows,
so merge temporaries are bounded rather than full-size vectors per addition.
Support-aware centering, signed tables and all deletion rows are preserved.

Use `rdistnn_jackknife_engine` to force `plugin`, `mata` or `naive` paths in
verification; this selector changes computation, not the statistical method.
The native plugin must be updated with its ado/Mata source. These are ordinary
double-precision safeguards, without arbitrary-precision runtime dependencies.
A single correction cannot guarantee correctly rounded results for every
severely ill-conditioned sum, recover multiplication error, or replace the
existing range and host-representability checks.

Tests cover literal cancellation identities, explicit deletions, scalar/paired
signals, small samples, zero support, unused extremes, buffer boundaries,
reordered targets, estimator covariance and saved prediction. Native acceptance
checks require the actual plugin; discrete decisions remain exact, while
library-dependent continuous values use justified tolerances. The
[test-suite map](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/README.md)
and [crosswalk](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/R_TEST_CROSSWALK.md)
identify the relevant suites.

## Weighted quantiles and diagnostic ranges

Weighted quantiles stably sort by value then input position. R, native C++ and
Mata use the same serial Neumaier total and prefix calculation in sorted order,
selecting the first observation whose rounded prefix divided by the rounded
total reaches the requested probability. There is no boundary tolerance. At
probability zero this retains the first sorted observation, even with zero
weight. Weights must be finite and nonnegative with finite positive mass.
Exact represented DNN tables have shared selected-index fixtures, checked
against independently rounded exact totals and prefixes during generation.

NP outcome SD is a scaled weighted norm. Binary scaling and an observed origin
retain small differences around large outcome levels. Zero-weight outcomes
cannot determine the scale, and the SD is restored directly: its variance need
not itself be representable. Physical CATE/CLATE radii reuse stable ordinary
Euclidean distances in original coordinates, regardless of the fitted metric.
Fresh fits and saved-fit prediction use the same helpers. Diagnostic distance
work preserves the recorded localization engine. All calculations use ordinary
double precision and linear workspace, without a new runtime dependency.

## Portable exact rank tables

`_rdistnn_rank_weights.ado` and the mirrored native `rank_weight_tables.hpp`
share R's exact order-statistic PMF recurrence. Anchor each PMF at its mode,
recurse in both directions with adjacent combinatorial ratios, normalize in
rank order with serial Neumaier addition, and average with retained correction
state. A separate mantissa and binary exponent prevents a subnormal recurrence
from sticking at the smallest positive double. Work is O(n + k*(n-s+1)), workspace O(n),
and k=s uses uniform weights directly. Exact means the finite-subset law,
computed in practical double precision. Native, R and Mata represented-table
fixtures and independent subset/combinatorial oracles cover both accuracy and
boundary decisions.

Asymptotic Beta tables retain their established host library routines. Their
continuous weights can differ in rounding, so an exceptionally close quantile
boundary can still select different observations across hosts. Identical
represented tables use the shared exact index rule; no epsilon hides a discrete
difference. Ordinary SD and covariance now shift by an observed origin before
mean centering, rather than waiting for an exceptional numeric amplitude.
Deletion origins additionally check exact subtraction to retain cancellation
contributions such as 2^54 + 1 - 2^54.

Numerical module ownership is separate from numeric plot-data assembly in
`_rdistnn_plot_mata.ado`. Private Mata names use the package prefix, including
`rdistnn_log_binomial()` for the existing combinatorial reference formula.
Reloading core modules therefore preserves unrelated Mata functions.

The shared `_rdistnn_mata.ado` entry point loads arithmetic, covariance and
whitening helpers before distance/rank evaluation, weight maps and deletion
adapters. Those responsibilities live in `_rdistnn_accumulation.ado`,
`_rdistnn_covariance.ado`, `_rdistnn_geometry.ado`, `_rdistnn_distance_mata.ado`,
`_rdistnn_rank_weights.ado` and `_rdistnn_rank_loo.ado`. The native distance bridge
and plot-data adapters remain separate. Moving an implementation does not
change its numerical law or engine-selection contract.
