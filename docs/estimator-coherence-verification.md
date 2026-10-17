# Estimator coherence verification

This receipt covers the CATE, CLATE and NP summary/diagnostics/plot changes
following `23a0736` on `feature/clate`. Numerical estimation algorithms are
unchanged. RdistNN remains authoritative; the sibling source snapshot is
`6fe48273aaa8853c5736841c5a228f8c1592ded9`.

## Common behavior

- One metadata/reporting path produces actual full-target `r(table)` (9 by m),
  `r(summary)` (m by 6) and `r(level)` after fit and replay. Canonical stored
  estimates, SEs and covariance have the same names across estimators.
- One diagnostic selection/save path supports targets and candidates,
  returns common localization data, and preserves scientific extensions,
  original identities/order, active estimates and data on save failures.
- One plot dispatcher validates options, selects builders, exports files and
  returns selected plot data/counts or dashboard panel metadata. Unsupported
  options are rejected. Observation-level plots verify fitted values/sample
  and row order; stored-result views survive data changes and legacy fits.
- CLATE retains instrument diagnostics, components, signed/undefined ratios,
  status masks and paired influence; NP retains weighted outcome distributions.
  `nose` means missing precision. Zero CLATE variance preserves SE zero with
  no Wald bounds; mean/CATE zero variance retains point bounds. PLQR scope
  remains supplied-coefficient moment evaluation.

Full stored/summary tables retain undefined targets. Plot data contain only
finite estimate rows actually drawn; an all-undefined estimate plot returns
zero rows and no data matrix. This is distinct from legacy `r(n_plotted)`,
which continues counting selected targets/candidates/panels.

## Reference and regressions

`tools/generate_r_reporting_fixtures.R` loads sibling R source without rebuilding
compiled code. It records actual CATE and affine-ratio summary coefficients at
90 and 95 percent, R NP estimates/SEs under the same normal convention, and
actual zero-variance CLATE summary/plot output. R NP has no summary method.
Stata represents infinite z statistics as missing because its numeric matrices
cannot represent infinities; a nonzero mean estimate with zero SE retains p=0.

The new tests are `test-estimator-reporting.do`,
`test-estimator-diagnostics.do` and `test-estimator-plots.do`. They cover
actual shapes/values rather than accepting Stata's missing-scalar fallback
for absent returned matrices. Existing numerical goldens and postestimation
regressions remain registered. The clean-install check also verifies the
seven shared helpers and all three common summary/localization/plot paths.

Independent source review found and prompted regressions for overlapping
covariate names, preservation of zero SE in CLATE plot data, and invalid
styles being swallowed by dashboards. Runtime integration additionally
checked candidate-bound masking, automatic overlap axes, and missing status
rows from infeasible candidates; valid candidate intervals remain available.

## Acceptance evidence

Verified on 2026-10-17 with the installed Windows StataNow 19 SE executable,
using Stata 17 ado syntax:

| Check | Outcome |
|---|---|
| Full licensed runtime (`tools/ci/stata_runtime_check.ps1`) | 1,200 assertions passed; 0 failed |
| Isolated installation (`tools/ci/stata_install_check.ps1`) | `INSTALL_CHECK_PASSED`; all shared helpers resolved from temporary PLUS |
| R reporting fixture regeneration | Byte-identical output from the recorded sibling source |
| Static package integrity (`bash tools/ci/pkg_check.sh`) | Passed |
| Python maintenance (`python3 -m unittest discover -s tests/ci -v`) | 8 passed |
| Final source archive (`bash tools/ci/build_source_archive.sh build/coherence/source.zip`) | Built and passed required-file, forbidden-path and duplicate validation |
| Independent source review | No remaining material findings after regression fixes |

The archive is an ignored verification artifact, not a new public release.
Summary fixture SHA256:
`cd5a6e635cd683f3d4bbdd98904cd515388f98d524de4b7b04295eafaca4d0fc`.

PNG and SVG exports for CATE, CLATE and NP profiles were inspected. Their
confidence regions render as filled polygons with visible upper/lower bounds.
The NP outcome plot separately labels its weighted 10th–90th percentile band.
These are scientific equivalents of the applicable R views; R ggplot internals
and image snapshot identity are outside the comparison.
