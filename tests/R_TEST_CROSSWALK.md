# R-to-Stata Test Crosswalk

This file records which high-value `RdistNN` contracts are enforced by the
Stata suite. It is a selection guide, not a claim that every R-specific API has
a Stata equivalent.

Shared contracts are checked bidirectionally and against independent oracles;
neither implementation's output alone establishes accuracy. Existing tiny-Beta
fixtures cover bounded cases, while the
[asymptotic-tail support limitation](../docs/numerical-stability.md#asymptotic-tail-limits)
remains unresolved. Passing fixed-table deletion tests does not certify all
probability construction, averaging or signed-blending boundaries.

## Shared contracts

| Contract from the R suite | Stata evidence |
|---|---|
| L2 and Mahalanobis distances, deterministic rank ties, translation/scale invariance, and singular covariance | `test-mata-distances.do`, `test-mata-rank-r-mirror.do` |
| Exact DNN weights: normalization, closed forms, brute-force subsets, support boundaries, and permutation invariance | `test-mata-weights-exact.do`, `test-weights-dnn-r-invariants.do` |
| Asymptotic DNN weights: normalization, monotonicity, exact approximation, permutation invariance, and `k = s` uniformity | `test-mata-weights-asymp.do`, `test-r-parity-priority.do` |
| Mata tiny Beta interval masses and close-order TDNN coefficients against analytic/independent references | `test-weight-numerical-stability.do`, `fixtures/r_stability_fixtures.do` |
| Extreme L2 recovery and common exceptional boundaries, actual compiled/fallback ranks, full covariance and restored prediction | `test-distance-range-stability.do`, `test-stability-r-parity.do`, `fixtures/r_stability_fixtures.do` |
| Scaled L1/L2 collectively represented subnormal norms and final rounding; R `test-distance-scaled-underflow.R` | `test-distance-range-stability.do`, `src/plugin/test_rdistnn_distance.cpp`; analytic references and actual Mata/native route assertions. |
| L1, scaled L1/L2 with R sd/MAD/type-7 IQR/range/custom scales; batched weight exports and frozen saved geometry | `test-geometry-batch.do` |
| R quantile near-tie ranks/estimates, analytic extreme fitted scales, failed-scale guards and scaled L1 recovery | `test-geometry-numerical-stability.do`, `fixtures/r_geometry_fixtures.do` |
| Actual R C++ and forced fallback L1/L2/scaled distances/ranks, large offsets, double-accumulation ties, scalar/batch and duplicate targets | `test-native-distances.do`, `fixtures/r_native_distance_fixtures.do` |
| Forced native/Mata geometry across all five estimators and both weight methods, fitted scales, grids/covariance, candidates and cold saved prediction | `test-native-distances-estimators.do` |
| Exact-distance ranks across current estimators/engines, subsets and restored predictions | `test-distance-tie-invariance.do` |
| Exact-distance block averaging, survivor-weight redistribution after actual deletion, full paired/cross-target covariance, permutation invariance, all forced engines, range guards, saved reuse and target-specific plots; R's `test-distance-block-averaging.R` | `test-duplicate-blocks.do`; independent explicit deletion and regrouping references cover every metric, both weight methods and signed TDNN. |
| Exact-distance group weights, signed/unnormalized mass, top-k boundary slots, actual deletions, target-specific reuse and unsupported-policy rejection | `test-distance-ties.do`, `test-distance-tie-r-parity.do`, `fixtures/r_distance_tie_fixtures.do`; fixtures generated from sibling R with explicit subset/deletion references. |
| Corrupted cached groups rejected before reporting-only consumers, covariance, diagnostics and plots, including empty active data/cold Mata and preserved state after rejection | `test-cached-model-guard.do`; complements R retained-state corruption/reuse tests. |
| TDNN exact/asymptotic signed rank tables, support-aware covariance, scalar/batch equivalence and saved-state reuse | `test-tdnn.do`, `test-tdnn-saved.do` |
| Mata/plugin equivalence and signed-target parsing | `test-weights-dnn-plugin.do` |
| Scalar, vector, and multi-target rank-shift jackknife covariance against explicit deletion loops | `test-mata-jackknife.do`, `test-mata-jackknife-vector.do`, `test-mata-jackknife-golden.do` |
| Shared prefix/suffix deletions, retained naive reference, native/Mata parity, missing-engine fallback, native dimensions above matsize, all-estimator covariance/grid reuse | `test-mata-jackknife-prefix.do`, `test-jackknife-engines.do` |
| Support-aware numerical centering: zero-weight outliers, active ranks after deletion, target-specific covariance, unused extremes, asymptotic tails and analytical variance 9, matching R's `test-jackknife-support-centering.R` | `test-mata-jackknife-prefix.do`, `src/plugin/test_rdistnn_core.cpp` |
| Paper defaults: estimator `k = 1`, finite-sample `lambda = 0.4`, `s = floor(0.4*n^((kx+1)/(kx+2)))`, and sensitivity grid `(0.3, 0.4, 0.5, 0.75)` | `test-r-parity-priority.do`, `test-np-dnn-parity.do`, `test-cate-dnn-tuning.do` |
| NP weighted-mean identity, exact/asymptotic methods, Mahalanobis, markout, batch covariance, tuning grids, target/outcome diagnostics, plots, and fixed R batch/candidate values | `test-np-dnn-*.do`, `test-r-parity-priority.do` |
| CATE AIPW and precomputed-score identities; propensity clipping; batch covariance; tuning grids; diagnostics; and fixed R values | `test-cate-dnn-*.do`, `test-r-parity-priority.do` |
| CLATE augmented binary-IV/paired-signal routes, common clipping, one-sided restrictions, perfect-compliance CATE reduction, signed/zero denominators, paired frozen covariance and cross-target linearization | `test-clate-mata.do`, `test-clate-dnn-core.do`, `test-clate-dnn-errors.do`, `test-clate-r-parity.do` |
| CLATE exact/asymptotic/general-k and fixed Mahalanobis weights, single/batch targets, scalar/grid equality, paired components and instrument diagnostics against fixed R values | `test-clate-dnn-tuning.do`, `test-clate-r-parity.do`, `fixtures/r_clate_fixtures.do` |
| CLATE CATE-style plot families, instrument overlap and paired influence, unavailable ratios, estimates restore, graph export and first-stage diagnostic tables | `test-clate-dnn-postestimation.do`, `tools/ci/stata_install_check.ps1` |
| Experimental PLQR supplied-theta score formula, current inputs, explicit moment covariance, method/distance behavior, markout, and fixed R moment/covariance/SE/ESS values | `test-plqr-dnn-*.do`, `test-r-parity-priority.do` |
| CATE physical 90-percent-mass radius, normalized treated/control arm ESS, and top-five absolute score-contribution share | `test-cate-dnn-diagnostics.do`, `test-r-parity-priority.do` |
| Determinism under distance ties and target reordering | `test-mata-rank-r-mirror.do`, `test-weights-dnn-r-invariants.do`, `test-r-parity-priority.do` |
| Validation, naming aliases, missing-data handling, replay, graph construction, and isolated installation | command-specific error/naming tests, `test-np-dnn-parity.do`, `test-cate-dnn-postestimation.do`, `tools/ci/stata_install_check.ps1` |

The fixed cross-language dataset and values in `test-r-parity-priority.do` are
generated by `tools/generate_r_priority_fixtures.R` from the sibling R source.
The Stata runtime suite deliberately does not require R.

CLATE reference values are separately regenerated with
`Rscript tools/generate_r_clate_fixtures.R`. RdistNN remains authoritative:
do not adjust fixtures to match a Stata discrepancy. The R API supplies
marginal component/ratio covariance and retained paired deletions; the fixture
uses those deletions to compute Stata's joint cross-target covariance. Stata's
full-target `e(estimates)`/`e(se_targets)` preserve R's missing-target behavior;
its `e(b)` omits undefined target coefficients because Stata rejects missing
coefficient posts. The CATE-compatible 25-column diagnostic schema uses
instrument arms for CLATE, and the public CLATE tables relabel them.

## Reporting, reuse and audit fixtures

| Contract | Evidence |
| --- | --- |
| Fit/replay tables, levels, missing targets, `nose` and zero variance | `test-estimator-reporting.do`; `generate_r_reporting_fixtures.R` compares actual R summary methods at 90/95 percent, mapping infinities to Stata missing values. |
| Target/candidate diagnostics, retained data/state and contributions | `test-estimator-diagnostics.do` |
| Plot options/data, dashboards, sample guards and exported graphics | `test-estimator-plots.do`; inspect rendered SVG/PNG as well as returned data. |
| Numeric serialization, candidates, cold restores and all-engine prediction | `test-saved-prediction.do`, `test-candidates.do`, `test-snapshot-observations.do`; `generate_r_reuse_fixtures.R` supplies independent frozen-deletion joint covariance. PLQR keeps moments separate from supplied theta. |
| CLATE edge cases and complete interleaved component covariance | `test-clate-edge-contracts.do`; 35 seeded cases from `generate_r_clate_audit_fixtures.R` check subset-enumeration/integrated-Beta deletion oracles. |
| Fresh/restored diagnostic equality across all five geometries, methods and score routes | `test-diagnostic-route-contracts.do` |
| Shared fresh/restored contribution recovery and concentration | `test-diagnostic-contribution-range.do`; independent rational references cover finite-operand subtraction overflow, tiny weights, zero support and overflowing totals. |
| Caller state, genuine returns, cold restores and confidence bands | `test-reporting-export-contracts.do`, `test-observation-reconstruction.do`, `test-state-helpers.do` |
| Fresh logs/markers and process exit status | `tests/ci/test_premerge_audit.py` uses mocked processes to reject stale or merely echoed markers. |

The [historical index](../markdown_notes/history.md) preserves dated audit
receipts separately from this current coverage map. The [test map](README.md)
owns execution commands and fixture regeneration. Current snapshot plots
survive data changes. R ggplot
internals are outside functional-interface parity.

## R-only contracts

The following R tests do not currently map to a Stata command and should not be
ported as superficial assertions:

- custom distance, rank-weight, score, and tuning registries;
- general local GMM models and IV simulation/replication helpers;
- R container dispatch (`matrix`, `data.frame`, tibble), S3 methods, and
  `predict()` behavior;
- ggplot object internals and `vdiffr` visual snapshots;
- PLQR density-backed coefficient inference at the supplied evaluation vector
  (the current R command also does not perform local root solving);
- PLQR multi-target tuning grids, which are not yet a Stata feature.
- optional affine-ratio confidence-set inversion and R CLATE benchmark/learner
  helpers; the Stata CLATE workflow provides CATE-style pointwise Wald reporting.

When one of these features is added to Stata, promote its corresponding R
contract into the shared table and add deterministic behavioral tests.

| Cross-cutting contract | R test owner | Stata test owner |
| --- | --- | --- |
| Host-native signaling, literal/bounded details and caller policy | `test-error-reporting.R`, `test-reporting-recovery.R`: conditions/handlers and warning promotion | `test-message-reporting.do`: standard codes and ordinary/quiet/capture channels |
| Group reuse | `test-distance-tie-localization.R`: call-local preparation and validated retained groups | `test-group-reuse.do`: public-grid preparation counters and saved prediction reuse |
| Engine selection and failure propagation | `test-compiled-fallback-parity.R`: explicit compiled/portable selection; `test-reporting-recovery.R`: unexpected callback/resource failures propagate | `test-native-failure-contracts.do`: native capability/execution/copy boundaries; standalone C++ tests provide separate kernel/protocol controls |
| Plot failure propagation | `test-reporting-recovery.R`, `test-plot-recovery-contract.R`: recognized optional recovery, unexpected failure and warning controls | `test-plot-failure-contracts.do`: bound probes and required-panel failures |
| Ratio batch reporting | `test-error-reporting.R`: bounded issue summaries with retained statuses | `test-message-reporting.do`: 100-target summaries and explicit `nose` control |
