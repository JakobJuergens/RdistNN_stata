# Native Stata plugin builds

The ado/Mata implementation is the portable package baseline. It runs on
Stata 17 or later without a native compiler. The plugin computes scalar and
batched L1/L2 distances, with optional fitted inverse scales. A single matrix
bridge serves both weight commands, NP, TDNN, CATE, PLQR, CLATE and saved-fit
prediction. Scale fitting, stable ranking and rank-weight tables stay in the
shared Mata layer. Mahalanobis factorization uses Mata; projected norms
reuse the same L2 kernels. The separate compiled
jackknife consumes fixed ranks and signals, including frozen Mahalanobis ranks.

The binary committed at `src/rdistnn.plugin` is a Windows x86-64 PE DLL.
On every host, automatic dispatch attempts substantial distance batches and falls
back to Mata if the plugin is absent, incompatible, lacks the distance API or
fails. Install a matching native Linux or macOS binary on those hosts.

Select distances independently with `global rdistnn_distance_engine plugin`,
`mata` or `auto` (the default). Explicit `plugin` requires the new API and
reports failures; explicit `mata` selects the reference implementation.
These selectors apply to L1/L2, scaled variants and projected Mahalanobis
norms. Legacy inverse snapshots retain their original Mata calculation. In `auto`, `global rdistnn_use_plugin 0` skips native calls and `1`
requires them. The explicit distance selector takes precedence over this
shared control. Drop test globals afterward to restore normal dispatch.
`rdistnn_last_distance_engine` records the last actual distance calculation;
cached candidate/same-target results may perform no new calculation.

Automatic distances stay in Mata for scalar/small jobs. The native bridge is
considered when both the covariate count and target count are at least 16 and
`n * p * m >= 20,000,000`. These conservative bounds account for matrix
transfer overhead; they are a dispatch policy, not native API limits or a
guarantee of speed. Explicit `plugin` or auto with `rdistnn_use_plugin 1`
bypasses them. In a Windows timing check, n=5,000, p=8 and m=32 favored Mata,
whereas p=32 and m=256 favored the native path. Measure representative workloads
before selecting an explicit engine for performance.

The bridge probes `distance_api` once and caches its result in
`rdistnn_distance_plugin_api`. An older plugin can still supply `rank_loo`
while distances fall back to Mata. Drop the cache or restart Stata after
replacing/reloading a plugin. All calls use the fitted training geometry;
query rows and deletion samples never refit scales. Coordinate subtraction,
inverse scaling and double accumulation follow R's compiled kernels. Near-tie
fixtures also verify the portable paths' arithmetic and stable rank ordering.
Both native L2 implementations explicitly round squared coordinates before
addition so fused multiply-add cannot change those ties.

L2 recovery uses the same exceptional boundary in R, C++ and Mata: a squared
sum or scaled raw difference reaching `2^1023`, or a nonzero coordinate whose
square is subnormal or zero. Exceptional rows use a scaled sum of squares;
huge scaled differences use binary-halved operands when the inverse scale is
below one. This keeps finite norms usable while aligning arithmetic with
Stata's missing-value boundary. Ordinary rows retain the established operation
order. See the [numerical stability guide](../../docs/numerical-stability.md)
for the shared contract and compensated accumulation.

For jackknife calculations the default order is compiled prefix/suffix sums,
Mata prefix/suffix sums, then the retained naive deletion loop if the Mata fast
function is unavailable or produces nonfinite sums. Both fast engines reuse
one weight table at sample size n-1 across all targets/deletions. The naive
engine recomputes retained weights under each deletion for parity comparisons.
All engines keep signals, geometry, k and s fixed and preserve paired and
cross-target covariance.
Covariance consumers center signal columns before dispatch and keep deletions
centered through covariance formation, avoiding spurious variance from large
common signal levels. Raw deletion helpers, including the naive reference,
remain available for direct parity comparisons.

Use `global rdistnn_jackknife_engine plugin`, `mata` or `naive` to select an
engine explicitly. Explicit `plugin` and `mata` selections surface failures;
`auto` (or dropping the global) enables the fallback hierarchy. In auto mode,
`rdistnn_use_plugin 0` skips native calls and `1` requires them. An explicit
jackknife engine takes precedence over that shared plugin control.

The native distance and `rank_loo` operations exchange temporary Stata matrices created
through Mata. Their rows and output columns may exceed `c(matsize)`;
dispatch attempts native execution without changing the user's setting.
Older plugins that implement only `weights_l2` fall back to fast
Mata automatically. Install an updated binary to enable native deletions.
The bridge removes temporary matrices after each call and leaves the dataset
unchanged.

Build a native plugin from this directory:

```text
# Linux or macOS
make clean
make

# Windows x64 Developer Command Prompt
nmake /f Makefile.win clean
nmake /f Makefile.win
```

From the repository root, `tools/ci/plugin_build_check.ps1` performs the same
build in a temporary directory, executes deterministic core, weight, distance and
entry-point protocol tests, and verifies the output format. GitHub Actions runs it on Windows
x86-64, Linux x86-64, macOS arm64, and macOS x86-64. Native binaries are
workflow artifacts; they are not interchangeable across those platforms.
The Unix Makefile selects Stata's `OPUNIX` ABI declarations on Linux and
`APPLEMAC` declarations on macOS; the Windows Makefile retains the header's
`STWIN32` default.
Published GitHub releases attach all four binaries under platform-specific
names together with `SHA256SUMS.txt`. Rename the matching asset to
`rdistnn.plugin` and place it on the Stata adopath. The default dispatch uses
it when compatible and otherwise keeps the portable Mata path active.

## Shared matrix protocol and tests

`rdistnn_distance.hpp` contains the Rcpp-independent C++17 kernel derived from
the R package's `dist_l1.cpp` and `dist_l2.cpp`. Matrices are column-major:
training X is n by p, targets X0 are m by p, and output is n by m. The plugin
operations use preallocated temporary matrices:

```text
distance_api output
distance_batch l1|l2 X X0 inverse_scales|none output
```

The capability version is 1; inverse scales have shape 1 by p. Validation
rejects malformed dimensions, missing/nonfinite values, nonpositive inverse
scales and arithmetic overflow before publishing results. Temporary output
protects the caller's data. The older `weights_l2` operation remains available
for compatibility; public weights now reuse the shared distance/rank pipeline.

The legacy operation's incomplete-Beta continued fraction has a 200-iteration
budget and requires a finite positive converged state. The pure kernel throws
`beta_evaluation_error` on failure; the plugin boundary reports numerical
error 498 and directs callers to `weights_dnn`, whose Beta weights use Mata.
No output observations are written on this failure. Ordinary successful
iterations, exact endpoints and representable tiny tails keep their arithmetic.
Native tests cover failure propagation through interval evaluation and the
actual plugin entry point, including preservation of the caller's output.

Run `make check` or `nmake /f Makefile.win check` for native tests. With a
matching plugin and licensed Stata, run
`do "tests/tests.do" "test-native-distances*"` for R-backed kernel fixtures
and all-estimator forced native/Mata comparisons. The old-binary compatibility
file requires a separate process loading an older binary and setting
`global rdistnn_test_old_distance_plugin 1`; it otherwise skips explicitly.
Fixtures regenerate with `Rscript tools/generate_r_native_distance_fixtures.R`
from actual R native and forced fallback paths. Native builds, host protocol
tests and licensed runtime parity establish distinct parts of verification.
