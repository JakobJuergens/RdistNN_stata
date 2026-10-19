# Release Notes

## 0.1.0.9003 — unreleased development update

- Use shared ordered reductions for Mahalanobis projections and correct
  cancellation-sensitive covariance cross terms, including mean-subtraction
  residuals, to align scalar/batch and cross-platform behavior with R.
- Preserve TDNN target coefficient names after saved-fit prediction posting.
  Require every native accumulation fixture kernel and verify that the new
  rank-weight and plotting modules resolve from a clean installation.
  Refresh the shared numerical fixtures for the current arithmetic.
- Give distance/ranking, rank weights and deletion/covariance routines dedicated
  owners behind the existing shared loader, preserving commands and calculations.
- Test cold and repeated loading against independent numerical references,
  preserve unrelated caller functions, and clean up fixture state. Organize
  regression suites by behavior.
- Check unique private Mata definitions, complete installation inputs and
  canonical manifest paths, plus optional mirrored-header parity with R.
  Align documentation and the development version with the R companion.

## 0.1.0.9002 — earlier development snapshot

- Stabilize ordinary SD/covariance and mixed-sign large-offset deletion centering;
  preserve small cancellation contributions and all represented remote tails.
- Generate exact rank tables with the same bounded-workspace PMF recurrence in
  R, Mata and native code. Preserve stable signed-zero ties in native sorting.
- Preserve mu/tau/theta coefficient names after candidate/prediction posting;
  reject clipping bounds whose upper endpoint rounds to one. Clarify Robinson
  propensity sensitivity and update metric/batch-weight documentation.


- Align weighted-quantile boundary decisions with R/native compensated total
  and prefix arithmetic. Use exact binary selected-index fixtures.
- Recover finite physical CATE/CLATE radii and NP outcome SDs without squaring
  physical units outside their range; saved-fit prediction reuses the helpers.

- Align weighted signal means and rank-weight normalization with R through
  shared compensated arithmetic. Native and Mata deletion engines retain
  corrections through the prefix/suffix merge with bounded working storage;
  candidates and saved prediction use the same estimator adapters.
- Reject nonconverged Beta continued fractions in the legacy native
  `weights_l2` API. Report numerical error 498 with Mata recovery guidance
  before writing output; standard command weight tables already use Mata.
- Fit Mahalanobis whitening directly from the centered training design with
  the shared R singular-value rank cutoff. Preserve resolvable directions,
  reuse native/Mata L2 norms after projection, and keep legacy inverse snapshots.
- Centralize covariance centering, binary range normalization and symmetric
  assembly across DNN, TDNN, CATE, CLATE, PLQR moments and saved prediction.
  Recover representable extreme variances; use rescaling errors or existing
  CLATE target statuses for unrepresentable nonzero variance.
- Preserve full numeric precision when constructing reporting tables and
  returning saved-fit prediction uncertainty, including extreme scales.

- Match R quantile interpolation near rank-sensitive ties, recover finite SD
  scales after intermediate overflow/underflow, and reject failed scale fits
  before clamping. Recover scaled L1 differences consistently across Mata and
  C++, and retain signed TDNN covariance metadata after saved-fit prediction.

- Stabilize close-order TDNN coefficients, asymptotic Beta tail increments and
  extreme L2 norms consistently with R. Shared helpers serve native/Mata
  weights, estimation and saved prediction; ordinary-input distance rounding
  remains aligned. Shared compensated weighted/deletion sums are described in
  the numerical-stability guide with their precision and engine contracts.

- Add shared compiled scalar/batch L1, L2 and scaled-distance kernels for all
  estimators, weight commands and saved prediction. Keep fitted geometry,
  stable ranks and weight tables shared with Mata. Probe the distance API
  independently so older plugins retain their existing functionality and
  default distance dispatch falls back to Mata. Add forced-route controls,
  R native/fallback fixtures, overflow validation and end-to-end parity tests.
  Keep small distance workloads in Mata automatically to avoid matrix-transfer
  overhead; explicit native requirements bypass the size policy. Round L2
  squared coordinates explicitly to preserve stable ties on FMA-capable hosts.
- Add L1 and scaled L1/L2 geometry with R-compatible fitted scaling rules;
  freeze geometry across deletions and saved-fit prediction. Export batched
  weight, rank and distance matrices with original observation identifiers.
- Add `np_tdnn` and `weights_tdnn` with exact/asymptotic signed two-scale
  tables, joint delete-one covariance, numeric reporting, saved-model prediction,
  signed diagnostics and plots. All estimators share one arbitrary-rank-table
  native/Mata/reference deletion dispatcher.
- Promote supplied-prediction binary-IV CLATE and ordinary paired-signal ratio
  computation to supported workflows. Preserve pointwise inference assumptions
  and all scope limits; ratio confidence-set inversion remains R-only and
  experimental. PLQR and local IV-GMM retain their existing lifecycle status.

## 0.1.0.9001 — 2026-10-18

Development update on production `master`, aligned with the R companion.
Tagged releases and downloadable assets are published separately.

- Add `rdistnn_results`, `rdistnn_candidate` and `rdistnn_predict` with graph-free
  exports, explicit posting and named derived fits. Preserve active data and
  estimates by default and retain full CLATE target/status identities.
- Save versioned numeric model snapshots in estimation results. `.ster` files
  retain fitted numeric inputs/signals and support prediction, kernel and
  influence views without original active data. Candidate views reuse cached
  geometry; new-target prediction uses the fixed training metric.
- Align PLQR near-tie ranks with R and the shared distance-matrix convention,
  so fresh estimation and saved prediction resolve ties identically.
- Unify PLQR moment reporting, retaining supplied coefficients separately.
  Saved prediction supports batch moments; coefficient inference remains absent.
- Develop on unpublished local `development` branches. Production changes use
  reviewed pull requests and squash merges.

## 0.1.0.9000 — 2026-10-18

This post-0.1.0 development milestone follows
the R reference's CLATE/reporting milestone; the
[documentation crosswalk](docs/r-documentation-crosswalk.md) records the
remaining differences in feature coverage.

- Share compiled and Mata prefix/suffix frozen-score jackknife engines across
  NP, CATE, PLQR and CLATE. Reuse rank-weight tables across targets/deletions,
  retain explicit deletions for parity checks and as a final fallback, and
  preserve paired and cross-target covariance. Older/unavailable plugins use
  fast Mata automatically.
  Center covariance signals before accumulation to preserve zero variance
  for constant scores and accuracy under large common signal levels.
- Choose numerical centering references from each target's deletion support
  in compiled, Mata and naive jackknife paths. Exclude zero-weight outliers
  before subtraction and include ranks that can gain weight after deletion.
  Preserve the exact/asymptotic weight rules and R implementation parity.
- Allow native jackknife bridge matrices above `c(matsize)` without changing
  the user's setting; remove the artificial dispatch restriction and correct
  the related documentation.
- Reuse default-candidate covariance in NP/CATE and duplicate CLATE tuning
  fits. Short-circuit CLATE's preliminary zero-signal covariance without
  changing inference feasibility rules.

- Add experimental `clate_dnn` for augmented binary-IV conditional LATE,
  externally supplied nuisance predictions, explicit one-sided compliance,
  or precomputed paired reduced-form/first-stage signals.
- Reuse CATE's exact/asymptotic weights, general k, fixed L2/Mahalanobis
  geometry, single/batch targets and integrated tuning grids. Paired frozen
  deletions retain numerator/denominator and cross-target covariance.
- Add CLATE diagnostics and all CATE plot families with instrument-overlap
  and paired-score labels, plus installed help and an oracle-nuisance example.
- Preserve signed/zero first-stage results with component/status tables;
  undefined target coefficients are omitted from Stata's coefficient post
  while remaining visible in the full canonical target tables.
- Add analytic, validation, grid, postestimation and fixed R-reference tests.
  RdistNN remains the numerical reference implementation.

- Standardize CATE/CLATE/NP metadata headers, actual full-target `r(table)` and
  `r(summary)`, replay levels and canonical estimate/SE/covariance matrices.
- Share diagnostic selection/saving and localization returns, plot validation,
  export and selected plot data. Retain estimator-specific scientific fields.
- Check input values and observation order for all observation-level plots;
  repair confidence-band rendering across effects and regression plots.

- Expand the full Stata documentation into task-based articles aligned with
  R's article/reference coverage, add installed `help rdistnn`, strengthen all
  command references and provide six verified standalone examples. R-only
  interfaces and experimental inference boundaries are explicit.

- Align singular Mahalanobis geometry with R's SVD Moore-Penrose inverse,
  including targets outside the covariate span; retain the full-sample metric
  through deletions. This also corrects inherited CATE/regression behavior.
- Harden runtime/install/example verification against stale logs, echoed
  completion text and nonzero process exits. Repair Windows NMAKE encoding
  compatibility and expand independent numerical/postestimation audit tests.
- Avoid intermediate overflow in finite paired CLATE covariance and match
  R's unavailable-inference status while retaining finite component covariance.

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
