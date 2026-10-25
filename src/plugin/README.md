# Optional native acceleration

RdistNN runs on Stata 17 or later using portable ado/Mata code. The optional
plugin accelerates distances and frozen-rank jackknife deletions. These have
separate dispatchers: a small, low-dimensional fit can use Mata distances and
native jackknife calculations.

Optional plugin declaration may select the portable engine when the binary is
absent or cannot be dynamically loaded on the host. Unexpected declaration
failures, Break and resource errors preserve their original Stata return code;
they are not cached as unavailable acceleration. This loader policy is separate
from the native execution and capability-probe contracts.

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

If this fails, inspect `adopath` and the installation directory. Automatic
dispatch uses Mata when native acceleration is unavailable; an explicitly
required native engine reports an error instead.

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
| `rdistnn_jackknife_engine` | `auto`, `plugin`, `mata`, `naive`. Auto uses native prefix/suffix deletions when available; otherwise it tries Mata prefix/suffix sums, then explicit deletions. Native execution failures stop. The distance size thresholds do not apply. |
| `rdistnn_use_plugin` | In auto mode, `0` skips native calls and `1` requires them, bypassing distance size thresholds. An explicit engine selector takes precedence. |

Explicit `plugin` or `mata` selections surface failures rather than silently
switching engines. Automatic fallback applies to absent or unloadable binaries
and missing native capabilities. A native execution failure stops the
command, including in `auto` mode: allocation failures retain code 909, Break
retains code 1, and unexpected plugin failures use 498. Inspect the explanation
before choosing an explicit Mata rerun; the package does not retry that failure
automatically. See the [error-handling guide](../../docs/diagnostics-and-reporting.md#warnings-and-failures)
for `capture noisily` and saving `_rc`.

`rdistnn_last_distance_engine` records the last distance
calculation; a cached result may perform no new calculation. It does not report
the jackknife engine. `rdistnn_last_jackknife_engine` separately records the
actual deletion route. Remove diagnostic settings afterward:

```stata
macro drop rdistnn_distance_engine rdistnn_jackknife_engine rdistnn_use_plugin
```

Distances cover L1/L2, fitted scaled variants and projected Mahalanobis norms.
Mata owns scale fitting, Mahalanobis factorization, ranking and rank weights.
Capability probing runs on each dispatch. Restart Stata or reload the plugin
when replacing a loaded binary. Automatic thresholds account for transfer
costs and do not guarantee speed; benchmark the workload before forcing an engine.

## Build from source

For isolated Clang builds, Bear compilation databases, clang-tidy and
ASan/UBSan checks on Linux/WSL, use the
[developer analysis profiles](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md#developer-analysis-tools)
from the repository root. They leave this directory's release build settings
and installed binaries unchanged and run all six standalone native tests.
Licensed Stata native/Mata checks remain necessary for the host integration.

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
It uses `distance_arithmetic.hpp`, mirrored byte-for-byte from RdistNN's
`inst/include/`, for common exceptional subtraction and scaled-norm arithmetic.
Both native makefiles and the portability build helper include this header;
the package hygiene check with `--r-root` detects mirror drift.
The matrix protocol is column-major: training `X` is n by p, targets `X0`
are m by p, inverse scales are 1 by p, and the preallocated output is n by m:

```text
distance_api output
distance_batch l1|l2 X X0 inverse_scales|none output
loo_packed_api output
loo_packed input output
```

Capability version 1 validates dimensions, finite inputs and positive inverse
scales before publishing output. The deletion operation accepts n by m ranks,
n by q signals and an (n-1) by 1 deletion table, returning n by (m*q) output in
target/component blocks. Grouped deletion additionally accepts n by 1 integer
group identifiers in `1:n` and computes actual-observation deletions from the
original signals. Canonical contiguous groups retain original member order and
use linear work in n; supplied ranks that interleave groups or reorder members
use native direct regrouping.
Current Mata deletion bridges use the shared `_rdistnn_native_transport.ado`
module and require packed capability. An uncached `loo_packed_api` probe writes
capability version 1 to a 1 by 1 matrix, then `loo_packed` receives a balanced
frame. Missing packed capability uses the ordinary Mata fallback in automatic
mode and errors when native execution is required. A native execution failure
or malformed packed result stops with its error instead of selecting fallback.
Only the current packed deletion transport is exposed by the plugin; obsolete
unpacked entry points are removed. The numerical kernels remain independent
of transport and are tested directly.

The packed input contains five exact doubles, `[1,n,m,q,grouped]`, followed by
column-major arrays `ranks`, `signals`, `table`, and, when grouped is 1, `groups`.
The physical frame stores this sequence in row-major order, with
`r = ceil(sqrt(length))`, `c = ceil(length/r)` and zero padding. The output frame
holds the unchanged column-major `n` by `m*q` result in the same physical order,
with zero padding. Both frames have positive balanced dimensions (`r >= c`,
`r-c <= 1`) and fewer than `r` padding cells. Native code validates metadata,
integer bounds, checked products and sums, expected frame capacities, finite
inputs and zero input padding before allocating kernel arrays or publishing
output. Conformable input/output aliases are safe because all inputs are read
first. Packed transport calls the directly tested numerical kernels.

Balanced staging avoids Stata's expensive default stripe construction for tall
matrices; the logical algorithm and its linear complexity are unchanged.
Temporary Mata-created Stata matrices can exceed `c(matsize)`; the bridge removes
them after successful transfer and before propagating a native error code.
The dataset is unchanged.

The [numerical reference](../../docs/numerical-stability.md) owns arithmetic,
range recovery and support-aware covariance centering. Native tests cover the
frozen-rank deletions, distances, rounding, packed layouts and actual entry-point
failure propagation.

The plugin consumes rank-weight tables constructed by Mata. Forcing native
deletions therefore does not address the
[known asymptotic-tail construction limit](../../docs/numerical-stability.md#asymptotic-tail-limits).
Mata constructs ranks and rank-weight tables; the plugin consumes these inputs
without a separate sorting or weight-construction implementation.
With a matching plugin and licensed Stata, run
`do "tests/tests.do" "test-native-distances*"`.
See the [test map](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/README.md)
for fixture generation and the distinction between native, runtime and
cross-language checks.
