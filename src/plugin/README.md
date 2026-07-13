# Native Stata plugin builds

The ado/Mata implementation is the portable package baseline. It runs on
Stata 17 or later without a native compiler. The plugin accelerates only the
standalone `weights_dnn` L2 path; estimator and Mahalanobis paths continue to
use Mata.

The binary committed at `src/rdistnn.plugin` is a Windows x86-64 PE DLL.
On every host, `weights_dnn` attempts the plugin for eligible L2 calls and
catches load or execution failure before using Mata. Install the matching
native Linux or macOS plugin to make that fast path available. Set
`global rdistnn_use_plugin 1` to require the plugin and surface failures, or
set it to `0` to force Mata.

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
