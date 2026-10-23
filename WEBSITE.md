# Website copy for the Stata companion

Ready-to-adapt text for the software section of jakobjuergens.com, describing
production source version `0.1.3`, published on 2026-10-22. Installation
and user guides are available from the repository. The current version-pinned
release is [`v0.1.3`](https://github.com/JakobJuergens/RdistNN_stata/releases/tag/v0.1.3),
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

## 0.1.3 production-update announcement draft

RdistNN for Stata 0.1.3 averages localization weights over observations at
exactly equal represented distance under the selected metric. Shared Mata
helpers and the C++ plugin preserve each group's assigned mass and recompute
weights using the actual survivors after every jackknife deletion. Distinct
profiles at equal distance share weights. Refit earlier saved models to use
schema 3 and the current tie policy.

The paper assumes a continuous distance distribution, under which the old and
new rules agree almost surely. Inference with persistent distance ties remains
unestablished. Shared reporting gives bounded, actionable notices, respects
quiet output and preserves fatal resource and native-execution errors.

The [weight guide](https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/weights-and-tuning.md#numerical-limits)
describes remaining numerical limits, including asymptotic-tail rounding and
independently fitted Mahalanobis geometry.

Task-oriented guides, eight standalone examples and installed help support
treatment-effect and regression workflows. The paper supplies the underlying
inference conditions; numerical invariance does not establish those conditions.
Portable ado/Mata calculations and optional native acceleration remain available
on their supported platforms. Version 0.1.3 is available from production
`master` and the version-pinned `v0.1.3` release.

[Source updates](https://github.com/JakobJuergens/RdistNN_stata/blob/master/RELEASE_NOTES.md) ·
[Tagged releases](https://github.com/JakobJuergens/RdistNN_stata/releases)

Use the [verified release URL](https://github.com/JakobJuergens/RdistNN_stata/releases/tag/v0.1.3)
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
