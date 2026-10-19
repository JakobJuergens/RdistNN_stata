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
do "tests/tests.do" "test-premerge-*"
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
| Weighted CDF boundary indices, physical radius/SD ranges and restored predictions | `test-weighted-distribution-range.do`, `fixtures/weighted_cdf_fixtures.do` |
| Distances and deterministic ranks | `test-mata-distances.do`, `test-mata-rank-r-mirror.do` |
| Exact/asymptotic weights | `test-mata-weights-exact.do`, `test-mata-weights-asymp.do`, `test-weights-dnn-*.do` |
| Native Beta convergence, endpoints/tiny tails, and numerical-error propagation without output writes | `src/plugin/test_rdistnn_weights.cpp`, `src/plugin/test_rdistnn_distance_plugin.cpp` |
| Scalar/vector/multi-target jackknife references | `test-mata-jackknife*.do` |
| Compiled/fast Mata/retained naive deletion parity, matrix bridge, fallback hierarchy and all-estimator grid/covariance reuse | `test-mata-jackknife-prefix.do`, `test-jackknife-engines.do` |
| Support-aware centering, inactive extremes, deletion boundary ranks, asymptotic tails and the shared R variance-9 fixture | `test-mata-jackknife-prefix.do`, `src/plugin/test_rdistnn_core.cpp` |
| NP estimation, validation, and R parity | `test-np-dnn-*.do` |
| CATE score paths, batch/grid results, diagnostics, plots, and validation | `test-cate-dnn-*.do` |
| CLATE paired signals/deletions, binary IV, statuses, batch/grid results, R parity, and all CATE plot families | `test-clate-*.do` |
| Seeded R CLATE edge fixtures, independent subset/integrated-weight deletions, singular metric and full paired covariance | `test-premerge-clate.do` |
| Adversarial reporting/diagnostic/plot returns, complete caller state, cold restores, sample guards and exported bands | `test-premerge-postestimation.do` |
| Saved-model snapshots, graph-free reports, all-engine prediction, reusable candidates, cold .ster restores and independent R joint covariance | `test-fitted-results.do` |
| Experimental PLQR moment semantics and covariance | `test-plqr-dnn-*.do` |
| L1/scaled geometry, batch weight exports, frozen new-target prediction | `test-geometry-batch.do` |
| Rank-sensitive R quantiles, finite extreme fitted SD/MAD, failed-scale validation and scaled L1 recovery | `test-geometry-numerical-stability.do` |
| Centered-design SVD, near-collinear analytic ranks, singular/extreme geometry and saved prediction | `test-mahalanobis-factor.do` |
| Shared covariance, recoverable raw-sum overflow and subnormal aggregates, unavailable variance and estimator errors | `test-covariance-range.do` |
| Actual R native/fallback distance and rank fixtures; forced distance routes, invalid protocol inputs, overflow, matrix dimensions above matsize | `test-native-distances.do` |
| Native/Mata equivalence for every estimator, metric/scaling rule and weight method, tuning and saved-fit prediction | `test-native-distances-estimators.do` |
| Actual older-plugin capability probing, automatic fallback and forced-mode errors (separate old-binary process) | `test-native-distances-compat.do` |
| Two-scale signed regression/weights, full covariance, reporting and saved reuse | `test-tdnn.do`, `test-tdnn-saved.do` |
| Independent tiny Beta masses, close-order TDNN coefficients and affected estimator/prediction paths | `test-weight-numerical-stability.do`, `src/plugin/test_rdistnn_weights.cpp` |
| Extreme L2 recovery, shared arithmetic boundaries, exact ranks and genuine output overflow | `test-distance-range-stability.do`, `src/plugin/test_rdistnn_distance.cpp`, `src/plugin/test_rdistnn_distance_plugin.cpp` |
| Actual R native/fallback extreme-distance estimates, joint covariance and restored prediction | `test-stability-r-parity.do`, `fixtures/r_stability_fixtures.do` |
| Actual R native/fallback centered-design geometry, full rank decisions, extreme estimator covariance and restored prediction | `test-covariance-geometry-r-parity.do`, `fixtures/r_covariance_geometry_fixtures.do` |
| Fixed shared R values | `test-r-parity-priority.do` |
| End-to-end estimation and replay | `test-integration.do` |
| Cold-session CATE/NP plotting order | `test-plot-command-order.do` |
| Cold/repeated core loading, caller-function preservation and independent low-level references | `test-module-loading.do` |

The NP parity file also tests its stored-result diagnostics and plots on the
same fixed batch/grid fixture. These related checks share explicit setup in
one modest file; further splitting would duplicate that setup. The numerical and saved-fit contracts use behavior-based filenames.
The plot-order regression runs first in the full suite to exercise cold ado
autoloading; it also runs independently through the exact filename filter.

## Fixtures and controls

The [numerical parity contract](../CONTRIBUTING.md#numerical-parity-contract)
requires within-package native/Mata comparisons as well as shared R/Stata
comparisons. Establish forced routes independently of automatic fallback;
check discrete ranks, support, statuses and missingness exactly alongside
justified tolerances for numeric values. Matching two optimized paths alone
does not replace an independent reference.

`Rscript tools/generate_r_native_distance_fixtures.R` generates distance/rank
fixtures from the sibling R package's actual C++ functions and forced R
fallbacks, and checks independent scalar references. These include large
offsets, constant coordinates, scale clamping, duplicate/reordered targets and
near ties sensitive to double accumulation. The kernel test requires a real
new distance API; it does not count unavailable native execution as a pass.
The end-to-end test fixes the jackknife to Mata to isolate distance dispatch.
The compatibility test skips in the normal suite: execute it in a separate
Stata process that has loaded an older plugin and sets
`global rdistnn_test_old_distance_plugin 1` to verify actual old-binary behavior.

`Rscript tools/generate_r_stability_fixtures.R` generates the stability fixtures
from actual R compiled and forced fallback paths. It checks analytic Beta tails,
distance/rank parity, estimates, joint covariance and serialized prediction.
The fixture header records source hashes. Regeneration additionally uses Python 3
standard-library exact combinations and 90-digit Decimal arithmetic for a signed
TDNN oracle. Licensed tests require both forced
Mata and plugin distances; unavailable native execution is a failure.

The centering regressions in `test-mata-jackknife-prefix.do` run the compiled,
Mata, naive and automatic routes. They cover all estimator commands, opposite-
sign unused observations, local variation, target-specific and joint support,
unused extremes and asymptotic tails. The variance-9 and asymptotic-tail
fixtures are independently calculated and also checked in R's
`tests/testthat/test-jackknife-support-centering.R`. Native core tests separately
exercise extreme unused values before subtraction. Run the focused Stata suite
with `do "tests/tests.do" "test-mata-jackknife-prefix"` and native tests with
`make -C src/plugin check` on Unix or `nmake /f Makefile.win check` from the
plugin directory under a Windows developer prompt.

The CLATE reference arrays in `fixtures/r_clate_fixtures.do` are generated by
`Rscript tools/generate_r_clate_fixtures.R` from the live sibling RdistNN
source. Their header records the reference commit. Tests compare clipped
binary-IV signals through estimates, paired and cross-target covariance,
SEs, geometry, instrument diagnostics, weight arrays and grid values across
exact/asymptotic, general-k and Mahalanobis cases. Analytic tests additionally
use the four-row variance-34 fixture and independent deletion enumeration.

Shared data generators are in `test_framework.do`. The R values in
`fixtures/r_priority_fixtures.do` are generated by
`tools/generate_r_priority_fixtures.R` from the sibling R source. The runtime
suite does not require R. Preserve these golden values and tolerances during
editorial or layout work. [R_TEST_CROSSWALK.md](R_TEST_CROSSWALK.md) records
shared and R-only contracts.

The full Windows suite requires the bundled compatible plugin for direct
parser and forced plugin/Mata parity checks. Use a matching native binary for
those checks on another platform. The portable estimator/Mata paths remain
available independently; do not report a full plugin suite as passed when it
has not run.

`tests/ci/` contains Python repository-policy and package-maintenance checks.
The maintenance checks require each `test-*.do` file to appear exactly once
in the master runner and validate public source-archive contents. Run them
with `python3 -m unittest discover -s tests/ci -v`.

Use [CONTRIBUTING.md](../CONTRIBUTING.md) for the runtime and isolated-install
helpers, and [AGENTS.md](../AGENTS.md) for editing rules. Licensed Stata runtime
results, static checks, native core tests, and archive checks establish
different parts of the release; retain their individual receipts.


## Documentation workflows

`python3 tools/ci/documentation_check.py` validates Markdown links/anchors,
installed help routing/SMCL anchors, public-command coverage and seeded example
presence. Optional `--r-root PATH` compares all R article/reference-family
names with the public documentation crosswalk. It is a structural check;
source review and runnable examples establish semantic accuracy.

With licensed Windows Stata, `tools/ci/stata_documentation_check.ps1` runs all
eight standalone example do-files and renders the overview help. It requires
the dedicated completion marker for every example and stops on failure.
The isolated-install helper repeats all eight examples against a temporary PLUS
payload and verifies the installed overview. This docs-only workflow complements
the numerical regression suite rather than rerunning unchanged estimator tests.

Compensated accumulation is covered by `test-compensated-accumulation.do`
(literal native/Mata/reference anchors and block boundaries) and
`test-accumulation-r-parity.do` (actual R native/interpreted fixtures,
paired deletions, estimator covariance and restored prediction). Regenerate
`tests/fixtures/r_accumulation_fixtures.do` with
`Rscript tools/generate_r_accumulation_fixtures.R` against the current compiled
sibling development package. Hexadecimal fixtures preserve input bits.
The generator requires native sum, cross-product and rank-deletion kernels
before exporting; it checks their outputs against the forced R fallback.
After regeneration, run `do "tests/tests.do" "test-accumulation-r-parity"`
to check the updated references against Mata, native and explicit deletions.

The weighted-distribution suite tests finite huge/tiny physical radii, directly
restored outcome SDs, large offsets, zero-weight extremes, saved prediction
and exact CDF boundary indices. Regenerate both companion CDF fixtures with
`python3 tools/generate_weighted_cdf_fixtures.py`; use `--check` to verify them.
Its standard-library rational arithmetic is a build-time oracle, not a runtime
package dependency. Inputs preserve binary doubles and expected selected indices
are checked against independently rounded exact totals and prefixes.

## Installation and helper ownership

`python3 tools/ci/package_hygiene.py` checks every ado/help installation input,
manifest uniqueness/local paths, explicit Stata version declarations, and the
`rdistnn_` prefix and unique source ownership of typed Mata functions. It is part of `pkg_check.sh` and
has isolated negative fixtures in `tests/ci/test_package_hygiene.py`.
When the sibling R source is available, pass `--r-root ../RdistNN` to require
byte-identical shared C++ headers; an absent or different header fails.

The core Mata entry point loads numerical modules and the separate
`_rdistnn_plot_mata.ado` adapters. Each module owns its definitions and reload
guards. The combinatorial reference helper is `rdistnn_log_binomial()` in
`_rdistnn_rank_weights.ado`; it does not occupy or drop another package's
`lncomb()` symbol. These changes preserve the existing numerical formulas.
Regenerate exact-table fixtures with `Rscript tools/generate_exact_rank_fixtures.R`;
this developer command requires the compiled sibling R source. Installed
Stata commands and the regression suite need no R runtime.
