# Native Stata plugin builds

The ado/Mata implementation is the portable package baseline. It runs on
Stata 17 or later without a native compiler. The plugin accelerates standalone
`weights_dnn` L2 weights and frozen-score jackknife deletions for all four
estimators: NP, CATE, PLQR and CLATE. Estimator geometry and weight tables
continue to use Mata. The jackknife kernel consumes fixed ranks and signals,
so it also supports deletions under frozen Mahalanobis geometry.

The binary committed at `src/rdistnn.plugin` is a Windows x86-64 PE DLL.
On every host, `weights_dnn` attempts the plugin for eligible L2 calls and
catches load or execution failure before using Mata. Install the matching
native Linux or macOS plugin to make that fast path available. Set
`global rdistnn_use_plugin 1` to require the plugin and surface failures, or
set it to `0` to force Mata.

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

The native `rank_loo` operation exchanges temporary Stata matrices and is
therefore limited by `c(matsize)` in both rows and output columns. Automatic
dispatch uses fast Mata for larger jobs, without changing the user's matrix
limit. Older plugins that implement only `weights_l2` also fall back to fast
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
build in a temporary directory, executes deterministic exact/asymptotic core
tests, and verifies the output format. GitHub Actions runs it on Windows
x86-64, Linux x86-64, macOS arm64, and macOS x86-64. Native binaries are
workflow artifacts; they are not interchangeable across those platforms.
The Unix Makefile selects Stata's `OPUNIX` ABI declarations on Linux and
`APPLEMAC` declarations on macOS; the Windows Makefile retains the header's
`STWIN32` default.
Published GitHub releases attach all four binaries under platform-specific
names together with `SHA256SUMS.txt`. Rename the matching asset to
`rdistnn.plugin` and place it on the Stata adopath. The default dispatch uses
it when compatible and otherwise keeps the portable Mata path active.
