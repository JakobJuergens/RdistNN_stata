# Website copy for the Stata companion

Ready-to-adapt text for the software section of jakobjuergens.com, describing
production source version `0.1.5`, published on 2026-10-25. Installation
and user guides are available from the repository. The current version-pinned
release is [`v0.1.5`](https://github.com/JakobJuergens/RdistNN_stata/releases/tag/v0.1.5),
with a verified source ZIP, platform-specific plugins and checksum files.

## Software entry

**RdistNN for Stata** — Distributional nearest-neighbor methods for conditional
average treatment effects and nonparametric regression. The package provides
pointwise jackknife uncertainty, tuning grids, localization diagnostics, and
Stata-native plots. Supported CLATE provides binary-IV conditional effects
with paired frozen-score covariance. The package accepts supplied nuisance
predictions or precomputed signals and complements the R package. Stata 17 or
later is required.
PLQR moment evaluation is experimental and does not provide coefficient inference.

- [Repository and getting started](https://github.com/JakobJuergens/RdistNN_stata)
- [R companion](https://github.com/JakobJuergens/RdistNN)
- [Accompanying research](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)

## 0.1.5 production-update announcement draft

RdistNN for Stata 0.1.5 strengthens validation of paired saved results,
reuses default and aliased CLATE tuning weights, and consolidates package-owned
error rendering while preserving native failures and caller state. Current
models use schema 4; earlier or incomplete saved models require refitting.
Focused tests and shorter guides retain independent numerical references and
the documented scientific limits.

Exact represented distance ties continue to receive averaged localization
weights, including after jackknife deletion. The paper assumes a continuous
distance distribution; inference with persistent distance ties remains
unestablished. Reporting respects quiet output and preserves resource failures
and interruptions. Native development tools add isolated analysis and compiler
caching while avoiding generated-output copies.

The [weight guide](https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/weights-and-tuning.md#numerical-limits)
describes remaining numerical limits, including asymptotic-tail rounding and
independently fitted Mahalanobis geometry.

Task-oriented guides, eight standalone examples and installed help support
treatment-effect and regression workflows. The paper supplies the underlying
inference conditions; numerical invariance does not establish those conditions.
Portable ado/Mata calculations and optional native acceleration remain available
on their supported platforms. Version 0.1.5 is available from production
`master` and the version-pinned `v0.1.5` release.

[Source updates](https://github.com/JakobJuergens/RdistNN_stata/blob/master/RELEASE_NOTES.md) ·
[Tagged releases](https://github.com/JakobJuergens/RdistNN_stata/releases)

Use the [verified release URL](https://github.com/JakobJuergens/RdistNN_stata/releases/tag/v0.1.5)
for version-pinned downloads and the repository README as the primary Stata
documentation entry point; the Stata companion has no separate hosted
documentation site.

## User documentation

The package has a task-oriented [article collection](https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/index.md)
and an installed `help rdistnn` overview. Its guides cover all Stata commands,
with runnable examples and an R-to-Stata topic map. Supplied-prediction CLATE
is supported under its documented pointwise inference conditions; its Wald
reporting has no weak-instrument-uniform guarantee. PLQR
evaluates moments at supplied coefficients and provides no coefficient
inference. Documentation publication to the repository does not deploy the
website or publish a tagged package release.
