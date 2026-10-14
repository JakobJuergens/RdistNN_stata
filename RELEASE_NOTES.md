# Release Notes

## 0.1.0 — 2026-10-14

- Initial public release of CATE and nonparametric regression with DNN weights,
  pointwise jackknife uncertainty, tuning grids, diagnostics, and native plots.
- Experimental PLQR evaluates localized moments at supplied coefficients.
- Includes GPL-3.0-or-later licensing with explicit permission to link the
  native plugin with Stata.
- Fixes shared CATE/NP plot loading, including CATE-to-NP plotting order.
- Provides version-pinned installation, a source ZIP and checksum, and native
  plugin downloads for Windows x86-64, Linux x86-64, macOS arm64, and macOS x86-64.

## 2026-10-14: Initial public-release preparation

- Licensed the original package code under GPL-3.0-or-later, with additional
  permission to link its native plugin with Stata. The installed package and
  source ZIP include both notices; third-party copyright notices are preserved.

- Refocused onboarding and package metadata on CATE and nonparametric
  regression, retaining pointwise inference conditions and explicitly oracle
  example nuisances. Added a runnable example and detailed command guide.
- Kept PLQR supplied-coefficient moment evaluation explicitly experimental;
  it does not provide coefficient estimation or inference.
- Added contributor, release, website, and test-suite guidance. The source
  archive builder is shared with CI and includes public guides/examples while
  excluding ignored build output.
- Added focused test-runner patterns and fresh per-run counters, plus a check
  that every regression file is registered exactly once.
- Fixed shared plot-helper loading so a CATE plot followed by an NP plot works
  in a fresh session. The numerical plot-builder routines are unchanged.
- Corrected PowerShell verification paths for UNC checkouts and quoted Stata
  do-file arguments; normalized the development/test adopath for WSL shares.
  Estimator, Mata, native-plugin, and golden-fixture behavior is unchanged.

## 2026-07-12: Paper tuning alignment

- Changed the finite-sample paper default to `lambda = 0.4`.
- Changed the default sensitivity grid to `(0.3, 0.4, 0.5, 0.75)`.
- Preserved `s(lambda) = floor(lambda * n^((kx + 1)/(kx + 2)))` and all
  explicit `lambda()`, `lambdagrid()`, `s()`, `sgrid()`, and `sdefault()`
  overrides.
- Aligned NP, CATE, PLQR, and weight-command automatic orders with the
  sibling R package and the conditional Z-estimation paper.
- Regenerated cross-language fixtures from the sibling R source, including
  candidate order, resolved orders, estimates, standard errors, ESS, and
  diagnostics for NP and CATE grids.
- Updated active paper references to *Orthogonal Inference for Conditional
  Z-Estimation without Cross-Fitting*; the full citation is maintained by
  `citation("RdistNN")` in the sibling R package.
- Clarified the CATE inference contract: orthogonality controls first-order
  population nuisance drift, while same-sample reuse separately requires the
  paper's nuisance-rate and full-state stability conditions.
- Corrected the interpretation of the exact DNN kernel order: positive support
  ends at rank `n - s + k`, so larger `s` means tighter, not wider,
  localization when `n` and `k` are fixed.
- Recast `plqr_dnn` as an experimental localized-moment evaluator with
  canonical `g()` and `thetaeval()` inputs, temporary `q()`/`beta()` aliases,
  explicit moment/covariance objects, and no coefficient-style display.
- Added paper CATE diagnostics for 90-percent physical radius, normalized
  within-arm ESS, and top-five absolute score-contribution concentration, with
  deterministic R/Stata parity coverage.
- Made `weights_dnn` default to `k(1)` and accept `s=n`; completed and tested
  the `np_dnn, nose` estimate-only endpoint and replay contract.
- Removed committed CI logs and source archives. Workflows now publish
  finite-retention artifacts and attach validated archives to releases without
  bot commits to the default branch.
- Documented that the universal ado/Mata path is portable while the bundled
  plugin is Windows x86-64, and added native plugin compile checks for Windows
  x86-64, Linux x86-64, macOS arm64, and macOS x86-64. Eligible L2 calls try
  the installed plugin on every host and retain the Mata fallback. Release
  workflows publish all four native binaries with SHA-256 checksums.
- Removed duplicate root plugin sources and the generated binary under
  `src/plugin/`, made `src/plugin/` the canonical native-source location,
  removed obsolete workflow repository-write paths, and added static hygiene
  guards for these contracts.
- Fixed native Unix builds to select Stata's `OPUNIX` declarations on Linux
  and `APPLEMAC` declarations on macOS instead of inheriting the Windows-only
  `__declspec(dllexport)` default.

### Verification

Run from the indicated repository root on 2026-07-12:

| Repository | Command | Outcome |
|---|---|---|
| `RdistNN` | `Rscript -e "devtools::document(); devtools::test()"` | 3,909 passed; 5 expected skips; 0 failures; 0 warnings |
| `RdistNN_stata` | `Rscript tools/generate_r_priority_fixtures.R` | Wrote 37 deterministic matrices from current sibling source |
| `RdistNN_stata` | `powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_runtime_check.ps1` | 612 tests passed |
| `RdistNN_stata` | `bash tools/ci/pkg_check.sh` | Static package check passed |
| `RdistNN_stata` | `powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_install_check.ps1` | Fresh isolated `net install` check passed |
| `RdistNN_stata` | `bash tools/ci/archive_check.sh /tmp/rdistnn-source-validation-final.zip` | Fresh source archive passed required-file, forbidden-path, and duplicate checks |
| `RdistNN_stata` | `git diff --check` | Passed |
