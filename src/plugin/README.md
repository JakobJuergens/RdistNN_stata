# Optional native acceleration

RdistNN runs on Stata 17 or later using portable ado/Mata code. The optional
plugin accelerates distances and frozen-rank jackknife deletions. These have
separate dispatchers: a small, low-dimensional fit can use Mata distances and
native jackknife calculations.

## Install a matching binary

The bundled `src/rdistnn.plugin` is a **Windows x64** binary. Linux and macOS
need their own binary. Check the assets for the matching package version on
the [release page](https://github.com/JakobJuergens/RdistNN_stata/releases).
The native release workflow produces:

| Host | Asset name |
| --- | --- |
| Windows x64 | `rdistnn-plugin-Windows-x64.plugin` |
| Linux x64 | `rdistnn-plugin-Linux-x64.plugin` |
| macOS Apple silicon | `rdistnn-plugin-macOS-arm64.plugin` |
| macOS Intel | `rdistnn-plugin-macOS-x64.plugin` |

Verify the download against that release's `SHA256SUMS.txt`, rename it to
`rdistnn.plugin`, and place it on Stata's adopath. Restart Stata after replacing
a loaded binary. If the release has no suitable binary, build from source below
or use the Mata fallback.

To locate the binary Stata will find:

```stata
findfile rdistnn.plugin
display "`r(fn)'"
```

If this fails, inspect `adopath` and the installation directory. An incompatible
or older binary normally triggers automatic fallback; an explicitly required
native engine reports an error instead.

## Select or diagnose an engine

Leave the selectors unset for automatic dispatch. For diagnosis, set them
independently before running the estimation command:

```stata
global rdistnn_distance_engine mata
global rdistnn_jackknife_engine plugin
```

| Control | Values and automatic behavior |
| --- | --- |
| `rdistnn_distance_engine` | `auto`, `mata`, `plugin`. Auto considers native batches when `p >= 16`, `m >= 16` and `n*p*m >= 20000000`; otherwise it uses Mata. |
| `rdistnn_jackknife_engine` | `auto`, `plugin`, `mata`, `naive`. Auto attempts native prefix/suffix deletions, then Mata prefix/suffix sums, then explicit deletions. The distance size thresholds do not apply. |
| `rdistnn_use_plugin` | In auto mode, `0` skips native calls and `1` requires them, bypassing distance size thresholds. An explicit engine selector takes precedence. |

Explicit `plugin` or `mata` selections surface failures rather than silently
switching engines. `rdistnn_last_distance_engine` records the last distance
calculation; a cached result may perform no new calculation. It does not report
the jackknife engine. Remove diagnostic settings afterward:

```stata
macro drop rdistnn_distance_engine rdistnn_jackknife_engine rdistnn_use_plugin
```

Distances cover L1/L2, fitted scaled variants and projected Mahalanobis norms.
Mata owns scale fitting, Mahalanobis factorization, ranking and rank weights.
Legacy inverse snapshots retain their Mata calculation. An older plugin may
support `rank_loo` without the distance API, so the two dispatchers can choose
different engines. Capability probing is cached in `rdistnn_distance_plugin_api`;
restart Stata when changing binaries. Automatic thresholds account for transfer
costs and do not guarantee speed; benchmark the workload before forcing an engine.

## Build from source

From this directory, use a C++17 compiler:

```sh
# Linux or macOS
make
make check
```

On Windows, use an x64 Visual Studio Developer Command Prompt:

```text
nmake /f Makefile.win
nmake /f Makefile.win check
```

Run the corresponding `clean` target before changing compiler settings. The
output is `rdistnn.plugin`; install it as above. Binaries are specific to their
OS and architecture. Unix builds select Stata's `OPUNIX` ABI on Linux and
`APPLEMAC` on macOS; Windows uses `STWIN32`.

From the repository root, the cross-platform build helper builds and tests in
a temporary directory and checks the output format:

```powershell
pwsh -NoProfile -File tools/ci/plugin_build_check.ps1
```

## Developer protocol and verification

`rdistnn_distance.hpp` contains the shared Rcpp-independent distance kernels.
The matrix protocol is column-major: training `X` is n by p, targets `X0`
are m by p, inverse scales are 1 by p, and the preallocated output is n by m:

```text
distance_api output
distance_batch l1|l2 X X0 inverse_scales|none output
```

Capability version 1 validates dimensions, finite inputs and positive inverse
scales before publishing output. `rank_loo` accepts fixed ranks, signals and
a deletion weight table. Temporary Mata-created Stata matrices can exceed
`c(matsize)`; the bridge removes them and leaves the dataset unchanged.

The legacy `weights_l2` operation remains for compatibility. Its incomplete-Beta
iteration has a 200-step budget and requires finite positive convergence;
`beta_evaluation_error` becomes error 498 without writing output observations.
The error directs callers to `weights_dnn`, whose Beta weights use Mata.
An older binary without `rank_loo` falls back to Mata in automatic mode.

The [numerical reference](../../docs/numerical-stability.md) owns arithmetic,
range recovery and support-aware covariance centering. Native tests cover the
core, distances, rounding, weights and actual entry-point failure propagation.
With a matching plugin and licensed Stata, run
`do "tests/tests.do" "test-native-distances*"`. Older-binary compatibility needs
a separate Stata process and `global rdistnn_test_old_distance_plugin 1`.
See the [test map](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/README.md)
for fixture generation and the distinction between native, runtime and
cross-language checks.
