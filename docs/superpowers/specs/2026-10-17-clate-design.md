# CLATE parity with the Stata CATE workflow

The requested outcome is a usable CLATE workflow on `feature/clate`, followed
by verification, a commit, and a push to that branch. The numerical reference
is the sibling RdistNN implementation at `d655282`. The Stata baseline is
`c8957c6`. Existing CATE, NP, PLQR and plugin behavior must keep working.

## Estimator and interface

Add experimental `clate_dnn`, `clate_dnn_diagnostics`, and `clate_dnn_plot`.
The estimator accepts an outcome and externally supplied binary instrument,
binary treatment, instrument propensity, and four instrument-arm conditional
means: `dcol() acol() ecol() my0col() my1col() md0col() md1col()`.
`onesided` permits omitted `md0col()`, requires D <= A, and requires zero
supplied m_d0. A second route uses a positional reduced-form signal together
with `denominator()` and `pseudooutcome`. Never average individual ratios.

Mirror CATE's `x()/xcols()`, `x0()/x0matrix()`, if/in, exact/asymptotic
weights, general k, L2/frozen Mahalanobis geometry, scalar/paper/grid tuning,
pointwise level, and `nose` controls. Validate binary and probability inputs,
finite signals, option conflicts, and the complete estimation sample early.
Arm means condition on all adjustment variables; localization uses only the
reported X coordinates. Fit no nuisance learners.

For R in {Y,D}, construct
`phi_R = m_R1-m_R0 + A*(R-m_R1)/e - (1-A)*(R-m_R0)/(1-e)`
using one clipped instrument propensity. Localize both columns with identical
weights. Report the ratio of their weighted means. Freeze predictions, target
coordinates, metric, and s across paired rank-shift deletions. Keep numerator,
denominator and their cross-covariance together; use the full-sample ratio to
linearize centered component deletions, never delete-one ratios. Retain full
cross-target covariance in Stata, in addition to marginal paired components.

## Results and exceptional targets

Retain CATE's target, ESS, tuning and diagnostic conventions. Add canonical
`e(estimates)` (1 by m), `e(se_targets)` (1 by m), `e(components)` (m by 2),
`e(theta_vcov)` (m by m), `e(component_vcov)` (2m by 2m, interleaved Y,D),
and `e(status)` (m by 5). Status columns are ratio_status,
negative_denominator, near_zero_denominator, covariance_status, wald_available.
Ratio codes: 0=ok, 1=zero denominator, 2=numerical failure. Covariance codes:
0=ok, 1=not computed, 2=zero score variance, 3=numerical failure.
`denominatortol()` flags cancellation relative to the maximum absolute input
denominator signal; it never trims, regularizes, or changes an estimate.
Finite negative denominators remain computed and flagged.

Stata rejects missing values in `ereturn post`. Therefore `e(b)` contains
only defined target coefficients, named theta for a single target and theta1,
theta2, etc. for a batch; undefined targets remain missing in the canonical
full-target tables. Post `e(V)` only when covariance is available for every
posted coefficient. With no defined targets, retain sample and component
results without posting coefficients. Replay explicitly reports unavailable
targets and inference. Never substitute zero estimates for undefined ratios.

## Postestimation

Share the existing plotting builders with estimator-specific labels and
canonical CLATE estimates/SEs. Support estimates (all CATE styles), diagnostics,
tuning, overlap (instrument mass, arm ESS, propensity), weights (profile,
cumulative, kernel), paired influence, and dashboard; preserve selection,
saving, and estimation-result restoration. Data-dependent kernel and influence
plots verify both the original values and input row order. Influence is
`w*(phi_Y-theta*phi_D)`. Diagnostics use instrument arms and this paired score,
and retain component/status tables for first-stage assessment. Document that
these are pointwise Wald calculations without weak-IV-uniform guarantees.

## Acceptance

Add independent analytic and explicit-deletion tests, fixed R fixtures,
perfect-compliance reduction to CATE, signed/zero denominator cases, validation,
batch cross-covariance, methods/k/metrics, grid/scalar identity, replay,
estimates store/restore, plot data/graph/export behavior, and cold loading.
Register all `.do` tests. Update installed help, manifest, package checks,
command guide, release notes and a runnable CLATE example. Run licensed Stata
tests, static checks, repository maintenance tests, source archive checks and
isolated installation. Review the final branch, commit, and push only
`feature/clate`; no merge or release is requested.
