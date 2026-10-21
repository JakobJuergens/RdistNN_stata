# Website copy for the Stata companion

Ready-to-adapt text for the software section of jakobjuergens.com, describing
the `0.1.0.9005` documentation and maintenance update for production `master`. Installation
and user guides are available from the repository. The latest tagged release
remains `v0.1.0`; a tagged-release announcement should link to its verified
source ZIP and platform-specific plugin assets.

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

## Source-update announcement draft

RdistNN for Stata 0.1.0.9005 makes treatment-effect and regression workflows
easier to learn, with shorter guides, eight standalone examples and updated
installed help. Stored results identify the outcome and covariance method,
including candidate and saved-fit prediction results. The accompanying paper
supplies the theory; the guides explain inputs, outputs and pointwise inference
limits. Portable ado/Mata calculations and optional native acceleration remain
available on their supported platforms.

[Source updates](https://github.com/JakobJuergens/RdistNN_stata/blob/master/RELEASE_NOTES.md) ·
[Tagged releases](https://github.com/JakobJuergens/RdistNN_stata/releases)

For a future tagged-release announcement, add the verified release URL. Use
the repository README as the primary Stata documentation entry point; the
Stata companion has no separate hosted documentation site.

## User documentation

The package has a task-oriented [article collection](https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/index.md)
and an installed `help rdistnn` overview. Its guides cover all Stata commands,
with runnable examples and an R-to-Stata topic map. Supplied-prediction CLATE
is supported under its documented pointwise inference conditions; its Wald
reporting has no weak-instrument-uniform guarantee. PLQR
evaluates moments at supplied coefficients and provides no coefficient
inference. Documentation publication to the repository does not deploy the
website or publish a tagged package release.
