# Numerical stability and implementation plan

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

## Planned compensated weighted and prefix sums

This work is planned and has not been implemented. Evaluate one compensated
accumulation convention across native and Mata rank-weight mass, signed TDNN
weighted estimates and forward/backward deletion sums, together with the R
companion's shared rank-deletion and weighting layers. Reuse existing helpers
across estimators. Preserve support-aware centering, frozen geometry,
normalization and paired/cross-target covariance.

Keep correction state constant per active sum and reuse existing buffers.
Independent high-accuracy or explicit-deletion references must cover long
sums, mixed magnitudes/signs, near cancellation, constant signals and unused
extremes. Verify returned deletions, covariance, candidates and restored
prediction, and require the same discrete decisions across engines.

Benchmark representative and adversarial workloads before adoption. An
interpreted compensated loop can be slower than Mata's vectorized
`runningsum`/dot products or R's corresponding primitives. Record runtime and
peak-memory changes and justify where compensation is enabled. The
[test-suite map](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/README.md)
and [crosswalk](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/R_TEST_CROSSWALK.md)
identify the existing parity tests to extend.
