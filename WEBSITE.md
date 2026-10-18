# Website copy for the Stata companion

Ready-to-adapt text for the software section of jakobjuergens.com. Publish
the announcement after the GitHub release and its checked assets are live.

## Software entry

**RdistNN for Stata** — Distributional nearest-neighbor methods for conditional
average treatment effects and nonparametric regression. The package provides
pointwise jackknife uncertainty, tuning grids, localization diagnostics, and
Stata-native plots. Experimental CLATE supports binary-IV conditional effects
with paired frozen-score covariance. The package accepts supplied nuisance
predictions or precomputed signals and complements the R package. Stata 17 or
later is required.
PLQR moment evaluation is experimental and does not provide coefficient inference.

- [Repository and getting started](https://github.com/JakobJuergens/RdistNN_stata)
- [R companion](https://github.com/JakobJuergens/RdistNN)
- [Accompanying research](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)

## Release announcement

RdistNN for Stata 0.1.0.9001 adds experimental conditional LATE and shared reporting,
diagnostics, and plotting interfaces for treatment effects and nonparametric
regression. It combines DNN weights, pointwise uncertainty, and tuning
sensitivity in native Stata workflows. The examples label oracle
nuisances explicitly and explain the conditions needed for inference with
estimated nuisances. Experimental PLQR evaluates localized score moments at
supplied coefficients. The ado/Mata path is portable; optional native assets
accelerate standalone L2 weight calculations on supported platforms.

Add the published release URL here once its source ZIP and plugin assets have
been verified. Use the repository README as the primary Stata documentation
entry point; the Stata companion has no separate hosted documentation site.


## User documentation

The package has a task-oriented [article collection](https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/index.md)
and an installed `help rdistnn` overview. Its guides cover all Stata commands,
with runnable examples and an R-to-Stata topic map. CLATE remains experimental:
its pointwise Wald reporting has no weak-instrument-uniform guarantee. PLQR
evaluates moments at supplied coefficients and provides no coefficient
inference. Documentation publication to the repository does not deploy the
website or publish a tagged package release.
