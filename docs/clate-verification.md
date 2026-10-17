# CLATE implementation and verification

Accepted October 17, 2026 on `feature/clate`, from Stata baseline `c8957c6`.
RdistNN remains the numerical reference. The fixed fixtures record R commit
`6fe48273aaa8853c5736841c5a228f8c1592ded9`; its numerical R sources are
unchanged from the `d655282` implementation used by the design. The live
conditional-LATE remark and paired-covariance notes in the sibling
Cond_Z_Estimation repository supplied the scientific scope.

The implemented workflow includes supplied binary-IV nuisance predictions,
one-sided compliance, precomputed paired signals, single/batch targets,
exact/asymptotic DNN and general k, frozen Mahalanobis geometry, scalar and
grid tuning, paired component and cross-target covariance, first-stage
diagnostics, replay, and every CATE plot family. Full-target tables retain
undefined ratios; coefficient posts retain defined target identities.

| Acceptance check | Final evidence |
| --- | --- |
| Licensed Windows StataNow 19 full runtime suite | `stata_runtime_check.ps1`: all 950 assertions passed, including all six CLATE test files and the existing CATE, NP, PLQR, weights and plugin tests |
| Fixed R reference | 76 assertions compare estimates, components, paired/cross-target covariance, SEs, geometry, instrument mass, component uncertainty, weights and grid values across exact/general-k/asymptotic/Mahalanobis cases |
| Reproducible fixtures | `Rscript tools/generate_r_clate_fixtures.R` regenerates `tests/fixtures/r_clate_fixtures.do` byte-for-byte |
| Clean installation | Final `stata_install_check.ps1` passes with fresh temporary PLUS/PERSONAL directories, all public/private CLATE files resolved from PLUS, installed example and diagnostic/plot workflows executed |
| Static package integrity | `bash tools/ci/pkg_check.sh` passes, including public help pairing, manifest and CLATE inference contract |
| Maintenance and source packaging | All 8 `tests/ci` tests pass; source ZIP builder and `archive_check.sh` pass |
| Rendered plots | PNG/SVG estimate and dashboard exports inspected; positive bands, instrument labels and selected return tables checked; SVG confidence-region regression passes |
| Independent review | Both material findings corrected with failing-then-passing regressions: intermediate unused covariance overflow, and unavailable-Wald plot intervals. No other material finding remained |

The four-row analytic fixture has components (2,.5), ratio 4 and paired
variance 34. Independent finite-subset deletion tests verify general-k and
joint covariance. Perfect compliance reduces to CATE. Exceptional cases cover
signed/zero first stages, partial undefined batches with active covariance,
zero score variance, numerical overflow/underflow and `nose`.

Observation-level kernel/influence plots verify original input values and row
order. Stored estimate plots remain usable after data changes. Plotting masks
unavailable Wald intervals while keeping finite points and canonical zero SEs.
CLATE profile confidence regions retain an outline because the installed
StataNow renderer suppresses the region when `lcolor(none)` is used; a literal
standalone rendering probe and SVG regression established the workaround.

Numerical checks do not verify causal assumptions, nuisance stability/rates,
undersmoothing or weak-instrument-uniform validity. The installed help states
the same supplied-prediction and pointwise-inference scope as the R reference.
Optional R confidence-set inversion and benchmark/learner helpers extend
beyond this CATE-parity workflow; see the
[test crosswalk](https://github.com/JakobJuergens/RdistNN_stata/blob/feature/clate/tests/R_TEST_CROSSWALK.md).

Verification logs and review exports remain in ignored `build/clate/` and the
runtime log. The runnable [CLATE example](../examples/clate.do),
[design](superpowers/specs/2026-10-17-clate-design.md), and
[implementation plan](superpowers/plans/2026-10-17-clate.md) document use and scope.
