# Numerical implementation reference

RdistNN and its Stata companion share numerical contracts for distance-localized
estimation. C++, R and Mata implementations of a selected method must agree on
rank/tie order, support, target identities and computational statuses, with
justified tolerances for continuous results. Independent references establish
accuracy in addition to agreement between engines. See the
[contributor contract](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md#numerical-parity-contract)
and [native guide](../src/plugin/README.md) for verification and engine controls.

High-level ranks use exact represented selected-metric distance then original
row index. Target-specific groups contain every observation with the same
represented distance, including distinct raw covariate profiles and zero-weight
ranks beyond nominal support. Shared helpers average assigned rank weights over
observations in each group, preserving signed and unnormalized mass. No tolerance,
additional rounding, outcome, treatment or nuisance defines groups. R, C++ and
Mata share this rule; the Stata plugin consumes Mata's ranks and per-target
group IDs.
Bare supplied low-level ranks/orders retain strict-order semantics unless
distances accompany them. Current snapshots validate and reuse cached groups,
distances, ranks and `exact_distance_block_averaging_v1`; older policies require
refitting.

`_rdistnn_blocks.ado` owns block identification, mass averaging and original-signal
block deletions. Existing shared weight/deletion paths call it, so estimator adapters
do not maintain independent implementations. For a group at ranks `a:b`, average
the corresponding full-sample table entries. After deleting a member, average
the `n-1` entries at `a:(b-1)` across its remaining members; earlier groups retain
their ranks and later groups shift by one. Signed TDNN tables follow the same
rule without taking absolute values or renormalizing the block.

The optimized route accumulates original signals with unshifted, shifted and
survivor-averaged block weights. Prefix and suffix states outside the deleted
block merge with its actual survivor states, retaining all leading values and
corrections until the final sum. Rounded group-mean signals are not used.
Canonical contiguous groups retain original member order and require O(n*m*q)
work. Supplied ranks that interleave groups or reorder members within a group
use direct regrouping within the selected engine, requiring O(n^2*m*q) work.
Centering expands the union of deletion support to
entire intersecting target-specific distance groups before choosing a numerical origin.
Direct regrouping reduces block masses in original observation order and
signal products in localization-rank order. The native transport retains the
original identities even when signals are centered in rank order.
Mata and native-plugin routes are checked against explicit deletion references.
Public commands
warn that inference with persistent distance ties is not established. The paper
assumes a continuous distance distribution, under which the old and new rules
agree almost surely; no uniform MSE improvement is established.

## Stable elementary calculations

TDNN coefficients use equivalent logarithm/exponential formulas that avoid
subtracting a rounded value near one when orders are close. Asymptotic DNN and
TDNN probability increments select lower or direct upper Beta tails, avoiding
subtraction of two CDF values that have both rounded to one. These choices
preserve the selected mathematical weighting rule and its normalization.
The code does not impose an additional tail cutoff. Host Beta evaluation and
subsequent mixtures nevertheless have the unresolved numerical limitation
described under [asymptotic tails](#asymptotic-tail-limits).

Fitted IQR/MAD scales use R's guarded type-7 interpolation, including its
weighted-endpoint arithmetic. SD fitting uses common conservative range gates
before a quadratic sum can overflow or underflow. Exceptional columns use an
exact binary scale, shift by an observed value, and compensate the normalized
mean and squared-residual sums. This preserves variation around large common
offsets and aligns R/Mata arithmetic; ordinary columns retain their existing
SD calculations. Nonfinite fitted scales are rejected before the epsilon
floor is applied; only finite small or constant scales are clamped. Fitted
scales must themselves fit Stata's finite nonmissing range.

Ordinary L1/L2 differences retain subtract-then-scale rounding and ordered
coordinate accumulation; L2 rounds each square before adding it. Unscaled L2
recovers exceptional rows with scaled sum of squares when the accumulated
squared sum reaches `2^1023` or a nonzero coordinate's square is subnormal or
zero. Scaled L1/L2 also recover rows when a raw difference reaches `2^1023` in
magnitude or a nonzero scaled coordinate product becomes subnormal, including
zero. This common boundary covers both signs and Stata's narrower finite
nonmissing range; using each host's own overflow trigger could change neighbors.

Exceptional scaled rows retain coordinate products as binary mantissas and
exponents until the complete norm is restored. Thus individually unrepresentable
scaled coordinates can contribute to a representable L1 or L2 norm. For huge
raw differences, this recovery uses separately rounded binary quarters before
subtraction. Final binary powers are restored in two balanced steps to avoid
premature subnormal rounding. A coordinate is active when its raw difference and inverse
scale are both nonzero. With exactly one active coordinate, retain its directly
rounded scaled magnitude; the squared R L2 API then rounds its square separately.
With multiple active coordinates, an unsquared recovered norm is at least the
largest directly rounded scaled magnitude. No additional lower bound is applied
to multi-coordinate squared results. Recovery does not guarantee correct rounding
for every norm. It rescans only exceptional rows with constant accumulator state.
Genuine out-of-range outputs retain the documented failure or unavailable-value
behavior.

Shared support-aware centering protects deletion covariance from large common
signal levels. Paired CLATE covariance also uses scaled linearized scores;
these operations retain numerator/denominator dependence. Numerical stability
does not change inference assumptions or the interpretation of a near-zero
first stage.

## PLQR threshold comparisons

`_rdistnn_plqr_mata.ado` owns the supplied-coefficient indicator
`Y <= D*theta + g`; R and its native score kernels implement the same policy.
Even finite ordinary thresholds can lose a decisive residual through
cancellation. With `L = cols(D)+2`, `eps = 2^-52`, `eta = 2^-1074`, accumulate
`A = abs(g)+abs(Y)+sum(abs(rounded D_j*theta_j))` in column order. Keep the
ordinary comparison only when its threshold and margin are finite,
`A < 2^1022`, and the absolute rounded margin strictly exceeds
`8*(L+1)*eps*A + 32*L*eta`. This conservative envelope covers ordinary
product/sum rounding and gradual underflow, with headroom below Mata's
positive missing-value boundary. Ordinary dot-product error analysis requires
that overflow be excluded; see [Ogita, Rump and Oishi (2005)](https://ogilab.w.waseda.jp/ogita/math/doc/2005_OgRuOi.pdf).

All remaining comparisons use exact accumulation of the finite input products,
including products outside the floating-point range. Positive and negative
fixed-point vectors each contain 170 base-`2^26` limbs with bit origin `-2252`.
Significand splitting and immediate carries keep every integer operation exact
in binary64; capacity covers all terms allowed by matrix dimensions. The
highest differing limb determines the sign, and exact equality satisfies the
inclusive indicator. This applies the
[adaptive-predicate approach](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf)
with an [exact sum-of-products accumulator](https://publikationen.bibliothek.kit.edu/1000018822/1375602),
without a new runtime dependency. It assumes round-to-nearest binary64,
gradual underflow and no unsafe reassociation. BLAS evaluation may choose a
different fast route, but every resulting discrete comparison must agree.

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
old arithmetic lost a resolvable direction. SVD fitting and projection rounding
can change represented distance equality or ordering across packages or platforms,
even for well-conditioned full-rank designs far from the cutoff. Ranks,
exact-distance groups, weights and estimates can then differ materially.
The same projected norm kernel does not make independently fitted geometry
bitwise equivalent. Group averaging uses each fit's computed distances, with
no tolerance or additional rounding; it cannot reconcile different computed
distances. For example, centered rows `(2,1)`, `(-2,-1)`, `(1,2)`, `(-1,-2)`
are mathematically equidistant from zero, with squared Mahalanobis distance
`3/2`. A one-ULP computed distance difference can split this group: with
outcomes `(0,0,10,10)`, `k=1` and `s=3`, weights `(.5,.5,0,0)` give estimate
zero, while uniform weights give five.

Saved geometry stores a header with scale and rank (rank when `p>=2`), a
raw observed origin and the zero-padded square factor: `(p+2) x p`.
Target reordering, subsetting and serialized prediction reuse this fitted
geometry. Earlier inverse-metric snapshots require refitting.

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

Binary column normalization can nevertheless erase a small original entry or
cross-product. Sensitive off-diagonals involving normalized columns therefore
retry original-scale centering and compensated high/residual products, caching
each column once. Accept that value only when its products, sum and final
scaled result are finite and every nonzero-operand product remains normal;
otherwise keep the normalized answer. Separately rounded nonzero subnormal
products can cancel to zero even when their normalized cross-entry remains
representable, so a zero-only product check is insufficient. The common normal
boundary is `2^-1022`, matching
[R's `double.xmin`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/zMachine.html).
Ordinary covariance has no second product pass.
Exceptional full covariance retains O(n*q) scratch; marginal variances keep
their existing arithmetic. This repair is shared with R and does not promise
correctly rounded covariance for every ill-conditioned input. See the
[shared geometry/covariance design](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/18_geometry_and_covariance_arithmetic.md)
for its source basis and independent regression contract.

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
retain small differences around large outcome levels. Square-root weights are
formed before normalization to retain representable root contributions even
when a normalized input weight underflows. Zero-weight outcomes
cannot determine the scale, and the SD is restored directly: its variance need
not itself be representable. Physical CATE/CLATE radii reuse stable ordinary
Euclidean distances in original coordinates, regardless of the fitted metric.
Fresh fits and saved-fit prediction use the same helpers. Diagnostic distance
work preserves the recorded localization engine. All calculations use ordinary
double precision and linear workspace, without a new runtime dependency.

Fresh and restored CATE/CLATE concentration and influence views share
`rdistnn_score_share` in `_rdistnn_diagnostic_mata.ado`. Zero-weight rows are
inactive before risky arithmetic. Ordinary finite weighted contributions keep
their represented values. If finite residual operands overflow during
subtraction or weighting, binary mantissa/exponent descriptors preserve the
weighted result and permit a dimensionless concentration even when its
unscaled total is out of range. The same descriptors recover dimensional
influence contributions when representable. A genuinely nonfinite already
rounded `theta*phi_D` remains unavailable; this helper does not recover arbitrary
products of three inputs. R's paired jackknife-concentration diagnostics retain
their separate squared-pseudo-value definition.

## Portable exact rank tables

`_rdistnn_rank_weights.ado` shares R's exact order-statistic PMF recurrence.
The Stata plugin consumes the resulting tables; it has no separate weight
constructor. Anchor each PMF at its mode,
recurse in both directions with adjacent combinatorial ratios, normalize in
rank order with serial Neumaier addition, and average with retained correction
state. A separate mantissa and binary exponent prevents a subnormal recurrence
from sticking at the smallest positive double. Work is O(n + k*(n-s+1)), workspace O(n),
and k=s uses uniform weights directly. Exact means the finite-subset law,
computed in practical double precision. Native, R and Mata represented-table
fixtures and independent subset/combinatorial oracles cover both accuracy and
boundary decisions.

### Asymptotic tail limits

Asymptotic tables also use uniform weights when k=s, including the s=1
single-scale component of TDNN. Other Beta tables retain their established
host library routines. Rounded
interval endpoints can erase a tail mass or give R and Mata different zero
versus nonzero support at subnormal scales. Averaging order-statistic tables,
signed TDNN blending and normalization can lose further tiny masses. This
remains unresolved. Extreme signals can make such differences visible in
estimates and jackknife covariance. Fixing component probabilities alone would
not establish correctness of the final mixed table.

The plugin consumes Mata's constructed tables. Selecting a native deletion
engine or renormalizing a table cannot recover probability already lost in
construction. Exact tables use a different finite-population rule; switching
methods changes the calculation. Continuous weights can also differ at an
exceptionally close quantile boundary. Identical represented tables use the
shared exact index rule; no epsilon hides a discrete difference. Bounded
fixtures do not certify all tail scales, and this limitation does not relax
the numerical parity requirement.

### Origins for weighted deletions

Ordinary SD and covariance shift by an observed origin before
mean centering, rather than waiting for an exceptional numeric amplitude.
Deletion origins additionally check exact subtraction to retain cancellation
contributions such as 2^54 + 1 - 2^54. Supported columns whose entire amplitude
is subnormal use origin zero, preventing loss of an unrepresentable weighted
residual before the common origin is restored. When deletion weights are
available, zero is also used if the maximum absolute coefficient times the
supported amplitude is below the normal range, if the table mass is nonzero
and subnormal, or if nonzero mass times that amplitude is below the normal
range. Preserve separately rounded raw products and divide by mass afterward;
normalizing coefficients first changes the rounding. Exact zero signed mass
alone does not disable centering. Ordinary tiny tails with normal amplitude
and mass retain their observed origin. Shared weight summaries avoid repeated
table scans. Fast and direct routes use the same policy.

These safeguards do not guarantee exact product distributivity after arbitrary
block averaging or signed cancellation. Custom signed subnormal masses can
still overflow during normalization even when raw sums are finite; finite
double arithmetic does not guarantee partial recovery of every aggregate.

Numerical module ownership is separate from numeric plot-data assembly in
`_rdistnn_plot_mata.ado`. Private Mata names use the package prefix, including
`rdistnn_log_binomial()` for the existing combinatorial reference formula.
Reloading core modules therefore preserves unrelated Mata functions.

The shared `_rdistnn_mata.ado` entry point loads arithmetic, covariance and
whitening helpers before distance/rank evaluation, weight maps and deletion
adapters. Those responsibilities live in `_rdistnn_accumulation.ado`,
`_rdistnn_covariance.ado`, `_rdistnn_geometry.ado`, `_rdistnn_distance_mata.ado`,
`_rdistnn_rank_weights.ado`, `_rdistnn_blocks.ado` and `_rdistnn_rank_loo.ado`.
`_rdistnn_native_transport.ado` owns balanced packed deletion transfer and
current capability validation; numerical kernels retain the same logical inputs.
The native distance bridge and plot-data adapters remain separate. Moving an implementation does not
change its numerical law or engine-selection contract.

## Localization, diagnostics and saved-state ownership

The private `_rdistnn_tuning_plan.ado` resolves scalar/grid orders, default
insertion, duplicate computed IDs and feasibility. NP, CATE and CLATE use the
same plan. `_rdistnn_localization_mata.ado` fits geometry and constructs ordered
distances and ranks once; estimator commands retain signal construction,
covariance/status policy, parsing and labels. CLATE performs its paired
calculation directly from that plan.

`_rdistnn_diagnostic_mata.ado` owns matrix-input localization, NP outcome and
CATE/CLATE arm/contribution diagnostics for both fresh and restored routes.
The fresh CATE contribution-availability/reduction convention and restored-model guard
remain explicit at the shared boundary. Physical-distance diagnostics preserve
the engine recorded by localization.

`_rdistnn_model_fit.ado` coordinates current saved-model checks, numeric
evaluation and derived-result posting. `_rdistnn_model_check.ado` validates the
complete schema-4 snapshot before cached consumers, including ranks, distances,
per-target groups, selected policy, frozen geometry, stored arm inputs and
diagnostic notice facts/counts. Tie flags are checked during the existing
rank/group validation pass; weak flags retain the established concentration
criterion. These are diagnostic facts, not persisted warning-history flags.
`_rdistnn_model_decode.ado` decodes that current state;
`_rdistnn_model_mata.ado` evaluates the estimator. Posting, coefficient stripes,
descriptive metadata and uncertainty labels belong to `_rdistnn_model_post.ado`.
`_rdistnn_reuse.ado` owns the caller's estimation-state transaction, including
failures. Saved score provenance is retained; older or incomplete snapshots
require refitting.

Plot builders own rendering and delegate numeric reconstruction to the current
snapshot model. Older no-snapshot development fits require refitting.
Private modules stay flat under `src/`, are installed through
`rdistnn.pkg`, and own uniquely prefixed typed Mata functions and reload guards.

The Stata native bridge stages deletion inputs and outputs in balanced private
matrices to avoid default stripe creation on tall matrices. The dispatcher
requires packed capability. Missing capability permits automatic Mata dispatch
or a forced-native error. A failed packed execution or malformed output stops
with its reported return code; it does not trigger automatic fallback or retry
raw operations. The current directly tested raw C++ ABI operations share
the same kernels and preserve every represented payload value.
The [native protocol](../src/plugin/README.md#developer-protocol-and-verification)
describes capability and layout validation; transport does not change arithmetic
or engine requirements.

Observed-origin deletion centering also uses the actual table's maximum absolute
coefficient and compensated mass. It uses zero when the weighted supported
amplitude is subnormal, the nonzero mass is subnormal, or nonzero mass times
supported amplitude is subnormal. Exact zero mass alone retains centering, as
do ordinary tiny tails when maximum coefficient amplitude and mass are normal.
Products remain rounded before any normalization; these exceptional safeguards
do not claim exact distributivity after general signed block averaging.

### Reporting uses prepared geometry

Public geometry preparation derives exact-distance group IDs from the already
ordered ranks. Candidate weights and deletion kernels receive the same groups;
model posting retains those groups without rebuilding them. Current saved
models validate their groups against retained distances/ranks before prediction,
weights or plots use them. New target geometry prepares its own groups.
Standalone equality-feature primitives keep their generic semantics.

Nonfatal tie reporting reads the prepared groups and performs no distance
evaluation, sorting or group reconstruction. Native bridges distinguish absent
capability from a failing operation and preserve host return codes after
transport cleanup. The plugin contains every C++ exception at its ABI:
`std::bad_alloc` returns 909, direct invalid requests return 198, and unexpected
execution failures return 498.
