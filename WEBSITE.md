# Website copy for the Stata companion

Ready-to-adapt text for the software section of jakobjuergens.com. Publish
the announcement after the GitHub release and its checked assets are live.

## Software entry

**RdistNN for Stata** — Distributional nearest-neighbor methods for conditional
average treatment effects and nonparametric regression. The package provides
pointwise jackknife uncertainty, tuning grids, localization diagnostics, and
Stata-native plots. It accepts supplied nuisance predictions or precomputed
CATE scores and complements the R package. Stata 17 or later is required.
PLQR moment evaluation is experimental and does not provide coefficient inference.

- [Repository and getting started](https://github.com/JakobJuergens/RdistNN_stata)
- [R companion](https://github.com/JakobJuergens/RdistNN)
- [Accompanying research](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)

## Release announcement

RdistNN for Stata now has a public release focused on CATE and nonparametric
regression. It combines DNN weights, pointwise uncertainty, tuning sensitivity,
diagnostics, and plots in native Stata workflows. The examples label oracle
nuisances explicitly and explain the conditions needed for inference with
estimated nuisances. Experimental PLQR evaluates localized score moments at
supplied coefficients. The ado/Mata path is portable; optional native assets
accelerate standalone L2 weight calculations on supported platforms.

Add the published release URL here once its source ZIP and plugin assets have
been verified. Use the repository README as the primary Stata documentation
entry point; the Stata companion has no separate hosted documentation site.
