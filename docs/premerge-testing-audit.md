# Premerge testing audit

Historical audit of the feature snapshot below. Its branch, coverage and
experimental-status statements describe that snapshot. For the current
supported CLATE workflow and unchanged pointwise inference conditions, see
[conditional LATE](conditional-late.md) and the
[development notes](../RELEASE_NOTES.md).

Date: 2026-10-18. Starting feature snapshot: `5744599f11e61d24f21e5e1c258f1e9fb591b6ff`.
Production reference: `origin/master` at `c8957c6f3359701f84bf88727b1af40b1eafa48e`,
refreshed through the protected workspace fetch wrapper. The three feature
commits change 88 files. R reference: clean sibling RdistNN source at
`6fe48273aaa8853c5736841c5a228f8c1592ded9`.

Status: **audit complete; no unresolved merge-blocking defect in the tested
contract**. The corrected implementation is suitable for merging, retaining
CLATE's experimental inference status and the platform/version limits below.
The audit fixes and regression tests are maintained on `feature/clate`;
production has not been merged.

## Findings and corrections

| Priority | Finding | Correction and evidence |
| --- | --- | --- |
| High | Singular Mahalanobis geometry used `invsym()`, a generalized inverse that can change neighbor ranks for targets outside the data span. This predates the feature branch and is inherited by CLATE. | Match R's ordinary inverse with an SVD Moore-Penrose fallback and relative singular-value cutoff `1e-10`. The off-span R case failed five estimate/component/covariance/SE comparisons before the fix; all 228 initial new numerical assertions pass afterward. Literal projected distances, inverse, ranks and weights plus CATE/NP reductions extend the regression coverage. |
| Medium | Paired CLATE covariance could overflow in the unscaled dot product or by adding two large finite covariance matrices before halving. This affects extreme-scale inputs in the new CLATE implementation. | Scale centered deletions before the dot product and symmetrize without doubling the entries. R cases with finite component variance `7.2e307`/`5e307` first failed five comparisons. Corrected covariance preserves finite components and reproduces R's ratio-inference status, including R's effective-sample-size-scaled variance overflow check without forming that unsupported Stata value. |
| High | Runtime, installation and documentation wrappers ignored nonzero process exits; runtime `-KeepLog` could accept a previous run's log. The success parser also accepted echoed command text. These weaknesses predate CLATE. | Require zero process exit and a fresh rendered completion marker. `-KeepLog` archives the previous log as `tests.previous.log`. Anchor the positive-count success line. Behavioral mock-process tests first reproduced five false-pass scenarios and now reject all of them while accepting valid fresh runs. |
| Medium | The Windows native build failed with NMAKE `U1038: internal error: lexer` on UTF-8 Makefile comments under the installed Windows code page. | Keep `Makefile.win` comments ASCII. An isolated comment-only change made the same compiler build/core tests pass; the permanent fix passes a fresh MSVC build. The defect predates the feature branch. |
| Low | The static version-floor check rejected literal `version 16` only and accepted `15`/`16.1`. | Parse all numeric version declarations below 17. Behavioral mutation checks accept `17`/`17.1` and reject `15`/`16.1`. |
| Low | One dashboard absence assertion accepted a missing 1-by-1 matrix as evidence that `r(plot_data)` was absent. | Check actual membership in `r(matrices)`. New tests similarly establish genuine matrix membership, dimensions, values and missing masks. |

The metric correction intentionally changes affected singular-covariance
production results to agree with R. It is not a claim that previously stored
singular fits can be numerically reconstructed under the corrected geometry.
The ordinary full-rank production compatibility cases are checked separately.

## Numerical audit

The original 1,200-assertion licensed suite passed before the audit changes.
That pass did not expose the singular off-span geometry defect.

`tools/generate_r_clate_audit_fixtures.R` creates 35 seeded R reference cases
covering exact/asymptotic weights, general `k`, boundary orders, full-rank and
singular Mahalanobis geometry, ties, duplicate/reordered targets, signed/zero
and near-zero denominators, zero/singular paired covariance, extreme signal
scaling, omitted inference, supplied IV nuisances and one-sided compliance.
For supported small cases, generation separately enumerates exact subsets
and integrates beta order-statistic densities for asymptotic weights, then
forms every frozen deletion and the full interleaved component covariance.
Extreme-scale reference-only cases are identified in the generator.
Conditioning controls exercise ordinary inversion below the SVD cutoff,
fallback below machine epsilon, acceptance exactly at epsilon and a zero
covariance matrix. The exact-epsilon case exposed a strict-inequality error
in the first audit fix; changing it to `>=` passes the R boundary fixture.

`test-premerge-clate.do` compares full estimate/component/SE/status tables,
cross-target ratio covariance and all component-covariance blocks, including
their missing masks. Metamorphic checks cover affine numerator transforms,
common negative scaling, target permutations and duplicate covariance columns,
instrument relabeling, grids and invalid inputs. Fixtures retain their seed,
R version and R commit; the Stata suite needs no R installation.

The three older R fixture generators reproduce all existing numerical values.
Only the priority fixture's provenance version changes from `0.1.0` to the live
R version `0.1.0.9000`. The CLATE and reporting fixtures are byte-identical.
The new fixture generator also checks byte-identical regeneration.

Selected R source tests pass **769 expectations across 79 tests**, covering
affine-ratio estimation/methods, CLATE signals/estimation/tuning/diagnostics/
plots, Mahalanobis distance, and plot data/builder contracts. These use
`pkgload::load_all(..., compile=FALSE)`; they do not certify compiled R engines.

## Production and postestimation compatibility

An ignored `git archive origin/master` checkout supplies actual production
code without switching the working tree. Fresh production and feature fits
compare **264 stored numerical matrices across 36 fits**: NP, AIPW CATE,
precomputed CATE and PLQR moments; exact/asymptotic methods;
L2/full-rank Mahalanobis; requested/omitted uncertainty where supported.
PLQR has no `nose` option and remains a supplied-coefficient moment command.

Actual production `.ster` files are additionally checked through the new
replay, target diagnostics and estimate plots. All four historical files (NP and CATE fits from
that production snapshot) pass those checks without canonical reporting fields.

`test-premerge-postestimation.do` checks CATE/CLATE/NP selection ordering and
duplicates, candidate aliases/infeasibility, saved diagnostic contents, replay
and plot CI identities, configured confidence levels, dashboard panel metadata,
observation contributions, restricted samples and complete-case geometry,
input value/order guards, cold restores, reserved input names, undefined and
zero-variance target masks, plus success/error paths. It snapshots every `e()`
matrix and stripe, scalar, macro and `e(sample)`, and caller data values,
variable order/types/formats/labels. Positive-interval and `nose` SVG exports
have separate numerical CSV receipts for inspecting actual rendered bands.

## Executed verification

Runtime/test/tool inputs are recorded in the ignored
`build/premerge/verified-runtime-sha256.txt`: 116 files, manifest SHA-256
`e9283606005f0a4fa47c56c5ca8ec17533433d54aefc6fa81699c8f79766c6e8`.
The manifest was checked unchanged after the final full suite. Audit-receipt
edits do not change runtime inputs.

| Check | Result |
| --- | --- |
| Original licensed Stata suite | 1,200 assertions passed |
| New numerical audit after metric/covariance corrections | All 276 assertions across 35 R cases and additional identities passed |
| New postestimation audit | All 115 assertions passed |
| Final full licensed suite | All 1,591 assertions passed across 45 registered files; zero unexpected runner aborts |
| Python maintenance/adversarial harness suite | All 18 methods passed, including 10 adversarial/version-gate methods |
| Static package, documentation/live R inventory, workflow lint | Passed: 33 Markdown files, 199 links, 12 helps, 11 commands, 6 examples; 13 R articles and 14 reference families |
| Six documentation examples and installed overview | All six passed in a fresh checkout session; overview SMCL rendered successfully |
| Isolated temporary PLUS/PERSONAL installation and examples | Passed from checkout: installed command/help/private-helper discovery, estimation/reporting/plot checks and all six examples |
| Source ZIP integrity, byte equality and extracted-source installation | Passed: all 83 entries match checkout bytes; clean installation from extracted ZIP and all six examples passed |
| Linux x64 GCC native plugin build/core tests | Passed; ELF binary checked |
| Windows x64 MSVC native plugin build/core tests | Passed after ASCII fix; PE binary checked |
| Fresh MSVC plugin loaded in licensed Stata | All 58 direct parser, forced-plugin/Mata parity and sample-restriction assertions passed |
| Exported SVG geometry | Four filled bands and three `nose` profiles passed against numerical CSVs; removing a band in a copied fixture was correctly rejected |
| Linux core AddressSanitizer/UndefinedBehaviorSanitizer | Passed with leak detection disabled; LeakSanitizer is unavailable under the sandbox's ptrace environment |
| Selected R reference tests | 769 expectations / 79 tests passed |

Commands are the maintained helpers in [CONTRIBUTING](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md),
plus `python3 -m unittest discover -s tests/ci -v`,
`python3 tools/ci/documentation_check.py --r-root ../RdistNN`,
`actionlint .github/workflows/*.yml`, and
`Rscript tools/generate_r_clate_audit_fixtures.R <temporary-output>`.
The licensed suite uses `tools/ci/stata_runtime_check.ps1`; installation and
documentation use their corresponding `stata_*_check.ps1` helpers. Native
builds use `plugin_build_check.ps1` with the expected platform/architecture.
`python3 tools/ci/plot_svg_check.py` checks actual exported paths against their
numerical data; raster previews were also visually inspected. Raw logs,
compiler products, restored fits and plot receipts remain ignored under
`build/premerge/`; permanent regressions and the R generator are source files.
The final full-suite log is `full-final-tests.log`; focused numerical and
postestimation logs are `clate-final-green.log` and
`postestimation-final-green.log`. Dedicated installation/example receipts are
`install-checkout-final.log`, `install-archive-final.log` and
`documentation-examples-final.log`. Independent review found no additional
source defect; its stale fixture-count observation was corrected.

## Scope limits

Licensed runtime evidence is **Windows x64 StataNow 19.5 SE**. Source `version
17` declarations and static checks do not establish execution on Stata 17/18.
There is no local licensed Linux/macOS Stata runtime or local macOS compiler
evidence. The hosted native-build matrix includes those platforms, but its
configuration is not proof that a new hosted run passed. The checked-in
Windows plugin and a fresh Windows source build require separate runtime
evidence, recorded above when executed.

This audit verifies documented common functionality, not every R-specific API
or bitwise equality of BLAS-dependent near-ties. The R solver's reciprocal
condition estimate and Mata's explicit inverse norm can differ at a machine
precision threshold; that boundary is not certified by finite fixtures.
Custom R registries, TDNN/GMM APIs, optional confidence-set inversion and R
ggplot internals remain outside the Stata contract. PLQR is not a coefficient
estimator. CLATE's Wald inference remains experimental and conditional on the
documented first-stage, nuisance, deletion-stability, localization and bias
requirements; diagnostics do not validate those assumptions or provide
weak-instrument-uniform inference.
