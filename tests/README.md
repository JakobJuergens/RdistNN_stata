# Stata test-suite map

The regression suite already uses behavior-based files. Numerical references,
cross-language golden checks, command validation, and postestimation tests
remain separate. `tests.do` is the master runner; `test_framework.do` supplies
assertions, counters, and deterministic DGP programs.

## Run a suite

Change to the repository root in Stata. The master runner loads local commands,
the framework, and Mata, and resets counters on each invocation:

```stata
do "tests/tests.do"
do "tests/tests.do" "test-np-dnn-*"
do "tests/tests.do" "test-cate-dnn-*"
do "tests/tests.do" "test-clate-*"
do "tests/tests.do" "test-plqr-dnn-*"
do "tests/tests.do" "test-weights-dnn-*"
do "tests/tests.do" "test-mata-*"
do "tests/tests.do" "test-r-parity-priority"
do "tests/tests.do" "test-premerge-clate"
do "tests/tests.do" "test-saved-prediction"
do "tests/tests.do" "test-candidates"
do "tests/tests.do" "test-observation-reconstruction"
do "tests/tests.do" "test-reporting-export-contracts"
do "tests/tests.do" "test-canonical-ties"
do "tests/tests.do" "test-duplicate-blocks"
```

The optional argument is a Stata `strmatch()` wildcard matched against the
filename without `.do`. An unmatched pattern returns error 198 instead of
reporting an empty suite as passing. A focused run uses the same numerical
assertions as the full suite.

| Area | Test files |
| --- | --- |
| Translation-safe covariance/SD and local jackknife origins | `test-centering-translation.do` |
| Exact represented rank tables, CDF decisions and near-full subsets | `test-exact-rank-portability.do`, `fixtures/exact_rank_tables.csv` |
| CATE population signal behavior | `test-cate-signal-contract.do` |
| Derived coefficient stripes and verified estimation samples | `test-derived-fit-contract.do` |
| Canonical raw-covariate ordering of exact-distance ties, subset/permutation invariance and current saved prediction | `test-canonical-ties.do` |
| Exact-X block averaging, actual-observation deletion and survivor weights, naive/Mata/native parity, signed TDNN, paired covariance, extreme ranges, reuse and target-specific plots | `test-duplicate-blocks.do` |
| Original-signal block states, exact finite survivors and signed zero-mass covariance, genuine subnormals, unnormalized mass and forced engine execution | `test-block-numerical-regression.do`, native `test_rdistnn_blocks.cpp` |
| Weighted CDF boundary indices, physical radius/SD ranges and restored predictions | `test-weighted-distribution-range.do`, `fixtures/weighted_cdf_fixtures.do` |
| Distances and deterministic ranks | `test-mata-distances.do`, `test-mata-rank-r-mirror.do` |
| Exact/asymptotic weights | `test-mata-weights-exact.do`, `test-mata-weights-asymp.do`, `test-weights-dnn-*.do` |
| Native exact-table subset references and active protocol failure propagation without output writes | `src/plugin/test_rdistnn_weights.cpp`, `src/plugin/test_rdistnn_distance_plugin.cpp` |
| Scalar/vector/multi-target jackknife references | `test-mata-jackknife*.do` |
| Compiled/fast Mata/retained naive deletion parity, matrix bridge, fallback hierarchy and all-estimator grid/covariance reuse | `test-mata-jackknife-prefix.do`, `test-jackknife-engines.do` |
| Support-aware centering, inactive extremes, deletion boundary ranks, asymptotic tails and the shared R variance-9 fixture | `test-mata-jackknife-prefix.do`, `src/plugin/test_rdistnn_core.cpp` |
| NP estimation, validation, and R parity | `test-np-dnn-*.do` |
| CATE score paths, batch/grid results, diagnostics, plots, and validation | `test-cate-dnn-*.do` |
| CLATE paired signals/deletions, binary IV, statuses, batch/grid results, R parity, and all CATE plot families | `test-clate-*.do` |
| Seeded R CLATE edge fixtures, independent subset/integrated-weight deletions, singular metric and full paired covariance | `test-premerge-clate.do` |
| Adversarial reporting/diagnostic/plot returns, complete caller state, cold restores and exported bands | `test-reporting-export-contracts.do` |
| Current frozen sample geometry, signal/nuisance snapshots and independent observation kernels | `test-observation-reconstruction.do` |
| Complete e() comparison, independent caller-data snapshots and enumeration-order independence | `test-state-helpers.do`, `helpers/e-state.do`, `helpers/caller-data.do` |
| Saved-model snapshots, all-engine prediction, cold .ster restores, graph-free reports and independent R joint covariance | `test-saved-prediction.do` |
| Reusable standalone candidates, duplicate orders and parent preservation | `test-candidates.do` |
| Snapshot kernels and contributions after caller-data changes | `test-snapshot-observations.do` |
| Fresh/restored/candidate diagnostic equality across supplied/AIPW/IV routes, methods and geometries | `test-diagnostic-route-contracts.do` |
| Recoverable contribution ranges, zero support, overflowing totals and matching fresh/restored concentration | `test-diagnostic-contribution-range.do` |
| Current-result descriptive metadata and representable propensity clipping | `test-stored-result-metadata.do`, `test-clipping-validation.do` |
| Experimental PLQR moment semantics and covariance | `test-plqr-dnn-*.do` |
| L1/scaled geometry, batch weight exports, frozen new-target prediction | `test-geometry-batch.do` |
| Rank-sensitive R quantiles, finite extreme fitted SD/MAD, failed-scale validation and scaled L1 recovery | `test-geometry-numerical-stability.do` |
| Centered-design SVD, near-collinear analytic ranks, singular/extreme geometry and saved prediction | `test-mahalanobis-factor.do` |
| Shared covariance, recoverable raw-sum overflow and subnormal aggregates, unavailable variance and estimator errors | `test-covariance-range.do` |
| Actual R native/fallback distance and rank fixtures; forced distance routes, invalid protocol inputs, overflow, matrix dimensions above matsize | `test-native-distances.do` |
| Native/Mata equivalence for every estimator, metric/scaling rule and weight method, tuning and saved-fit prediction | `test-native-distances-estimators.do` |
| Two-scale signed regression/weights, full covariance, reporting and saved reuse | `test-tdnn.do`, `test-tdnn-saved.do` |
| Independent tiny Beta masses, close-order TDNN coefficients and affected estimator/prediction paths | `test-weight-numerical-stability.do`, `src/plugin/test_rdistnn_weights.cpp` |
| Extreme L2 recovery, shared arithmetic boundaries, exact ranks and genuine output overflow | `test-distance-range-stability.do`, `src/plugin/test_rdistnn_distance.cpp`, `src/plugin/test_rdistnn_distance_plugin.cpp` |
| Actual R native/fallback extreme-distance estimates, joint covariance and restored prediction | `test-stability-r-parity.do`, `fixtures/r_stability_fixtures.do` |
| Actual R native/fallback centered-design geometry, full rank decisions, extreme estimator covariance and restored prediction | `test-covariance-geometry-r-parity.do`, `fixtures/r_covariance_geometry_fixtures.do` |
| Fixed shared R values | `test-r-parity-priority.do` |
| End-to-end estimation and replay | `test-integration.do` |
| Cold-session CATE/NP plotting order | `test-plot-command-order.do` |
| Cold/repeated core loading, caller-function preservation and independent low-level references | `test-module-loading.do` |

The plot-command-order test runs first to exercise cold ado autoloading and
also runs through its exact filename filter. Shared data generators and
assertions live in `test_framework.do`.

## Fixtures and controls

Tiny-Beta and close-order fixtures are bounded numerical checks, not a
guarantee for every subnormal tail. The
[numerical reference](../docs/numerical-stability.md#asymptotic-tail-limits)
records unresolved R/Mata support differences and loss during mixtures.
Keep weight construction separate from fixed-table deletion verification.

Follow the [numerical parity contract](../CONTRIBUTING.md#numerical-parity-contract):
force and establish native/Mata routes, test fallback separately, retain
independent oracles and check discrete decisions exactly. Runtime needs no R;
regenerating R golden values requires the sibling compiled R source. Preserve
those values during editorial changes.

Run generators from the repository root only after an intentional source
contract change; each output records its reference identity where applicable:

| Fixture | Generator in `tools/` | Main consumer |
| --- | --- | --- |
| Native distances/ranks | `generate_r_native_distance_fixtures.R` | `test-native-distances*` |
| Extreme distances/Beta/TDNN | `generate_r_stability_fixtures.R` | `test-stability-r-parity` |
| Scales and quantiles | `generate_r_geometry_fixtures.R` | `test-geometry-numerical-stability` |
| Covariance/geometry | `generate_r_covariance_geometry_fixtures.R` | `test-covariance-geometry-r-parity` |
| Compensated reductions/deletions | `generate_r_accumulation_fixtures.R` | `test-accumulation-r-parity` |
| Exact rank tables | `generate_exact_rank_fixtures.R` | `test-exact-rank-portability` |
| CLATE and audit edge cases | `generate_r_clate_fixtures.R`, `generate_r_clate_audit_fixtures.R` | `test-clate-r-parity`, `test-premerge-clate` |
| Shared defaults/estimators | `generate_r_priority_fixtures.R` | `test-r-parity-priority` |
| Reporting and reusable fits | `generate_r_reporting_fixtures.R`, `generate_r_reuse_fixtures.R` | `test-estimator-reporting`, `test-saved-prediction`, `test-candidates` |
| Weighted CDF boundaries (both packages) | `generate_weighted_cdf_fixtures.py` | `test-weighted-distribution-range` |

Use `Rscript tools/GENERATOR.R` or `python3 tools/GENERATOR.py` as appropriate.
The CDF generator accepts `--check` and uses exact rational arithmetic at build
time. Stability generation additionally uses Python's exact combinations and
90-digit Decimal reference. Native/forced-fallback R paths and independent
analytic calculations establish fixture accuracy; matching seeds alone do not.
Hexadecimal accumulation fixtures preserve individual input bits and require
native sum, cross-product and rank-deletion kernels before export.

The full suite requires a matching plugin for direct parser and forced native
checks. Native absence is not a passing native comparison. The bundled binary
serves Windows x64; other hosts need a matching build. End-to-end native distance
tests hold the jackknife at Mata to isolate the distance dispatcher.

Distance selectors are `global rdistnn_distance_engine mata` or `plugin`;
jackknife selectors are `global rdistnn_jackknife_engine plugin`, `mata` or
`naive`. Auto uses the independent fallback rules in the
[plugin guide](../src/plugin/README.md#select-or-diagnose-an-engine). Reset test
globals afterward. Native tests run with `make -C src/plugin check`, or
`nmake /f Makefile.win check` from that directory under a Windows developer prompt.

[R_TEST_CROSSWALK.md](R_TEST_CROSSWALK.md) maps shared and R-only contracts,
including independent CLATE deletions and the variance-9 centering reference.

## Maintenance and documentation checks

[CONTRIBUTING.md](../CONTRIBUTING.md#verification) owns the static, licensed,
installation and documentation commands. `tests/ci/` checks package policy,
unique suite registration and source-archive contents with isolated negative
fixtures. Static documentation checks validate links, help routing and topic
coverage; source review and runnable examples establish semantic accuracy.
Native builds, licensed runtime, installed examples and archive checks provide
distinct evidence. Save their receipts in ignored working directories.

The core loading gateway is `_rdistnn_mata.ado`; source modules own definitions
and reload guards. The [numerical reference](../docs/numerical-stability.md)
records helper ownership. Tests protect caller functions: the combinatorial
reference is `rdistnn_log_binomial()` and must not occupy or drop `lncomb()`.

## Shared test ownership

`helpers/fitted-data.do` creates the deterministic sample used by independently
selectable saved-prediction, candidate and snapshot-observation suites.
`helpers/e-state.do` compares all matrix values and stripes without depending on
field enumeration order. Its full mode also compares macros, scalars and
`e(sample)`; `numeric` compares the numeric matrices for numerical checks.
Both modes preserve the active estimation state. Inert complete-e snapshots
create no estimation stores or dataset sample-marker variables. Caller-data values and variable
metadata use the separate `helpers/caller-data.do` adapter, composed with the
estimation comparator by `helpers/postestimation-state.do`.

Current metadata belongs to `test-stored-result-metadata.do`.
Each behavior suite appears exactly once in
`tests.do`, checked by `tests/ci/test_package_maintenance.py`. Run any named suite
through the master wildcard argument, including `test-saved-prediction`,
`test-candidates`, `test-observation-reconstruction` and
`test-reporting-export-contracts`. Reporting/export retains the fixed
`build/premerge/postestimation-*-profile.svg` filenames consumed by the SVG gate.

`test-packed-native-transport.do` checks balanced frame ordering and padding,
logical multi-target/component reconstruction, actual packed native execution,
and forced versus automatic behavior after a recognized packed failure. The
standalone `rdistnn_loo_plugin_test` exercises the real entry point with malformed
metadata, overflow, invalid ranks/groups, nonfinite values, padding/shape errors,
read failures, aliases and bitwise raw/packed comparisons. Legacy capability
compatibility is also covered by private audit runs against an older binary.

`test-weighted-underflow-origin.do` checks zero survivors, constant signals,
signed/unnormalized subnormal table mass, amplified rounded duplicate weights,
centered/raw paths and actual naive/Mata/native execution. It also protects
ordinary tiny-tail and exact zero-mass origin behavior.
