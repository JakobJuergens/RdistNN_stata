# Getting started

This walkthrough fits conditional average treatment effects (CATE) and
conditional means at selected covariate profiles. It uses simulated data with
**oracle nuisance values**. In an application, supply aligned predictions
from your own outcome and treatment models; the package does not fit them.
The [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
gives the statistical assumptions and results.

## Install and open help

Stata 17 or later is required:

```stata
net install rdistnn, from("https://raw.githubusercontent.com/JakobJuergens/RdistNN_stata/master/") replace
help rdistnn
```

See [installation alternatives](../README.md#installation) for source pins and
the optional native plugin. The [complete quick start](../examples/getting-started.do)
is a standalone do-file; the commands below are also runnable in sequence.
They clear the active dataset.

## Generate a small example

```stata
clear
set obs 300
set seed 42
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double e_oracle = invlogit(0.5*x1 - 0.5*x2)
gen byte D = runiform() < e_oracle
gen double m0_oracle = 1 + x1 - x2
gen double m1_oracle = m0_oracle + 0.5 + 0.5*sin(x1)
gen double Y = m0_oracle + D*(m1_oracle-m0_oracle) + rnormal()*0.5
```

Keep each nuisance prediction on the row of the observation it describes.
The reporting covariates in `xcols()` can be a subset of those used to fit
the nuisance models.

## Estimate CATE at selected profiles

Use `x0()` for one profile or `x0matrix()` for several. Matrix rows are
profiles; columns follow the order in `xcols()`.

```stata
matrix X0 = (0,0 \ 0.5,0 \ -0.5,0)
cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) tuningmode(grid)
matrix CATE_summary = r(summary)
matrix CATE_V = e(V)
estimates store cate_fit
matrix list CATE_summary
```

`r(summary)` has one row per profile, with estimates, standard errors and
pointwise confidence bounds. Copy `r()` results before another command
replaces them. `e(V)` retains covariance between profiles; it is not a
simultaneous confidence band. Replay can change the confidence level without
refitting:

```stata
cate_dnn, level(90)
cate_dnn_diagnostics, detail
cate_dnn_plot, type(estimates) name(cate_profile) replace
cate_dnn_plot, type(overlap) view(mass) name(cate_overlap) replace
```

Inspect treated/control support and weight concentration. Clipping a
propensity cannot create observations in a missing treatment arm.
[Reporting](diagnostics-and-reporting.md) explains the diagnostic columns;
[visualisation](visualising-fits.md) covers graph options and export.

If your pipeline already constructs a CATE pseudo-outcome, use the
[supplied-signal route](cate-and-regression.md#precomputed-pseudo-outcome).
A supplied signal is used unchanged and cannot supply treatment-overlap
information by itself.

## Fit a conditional mean

`np_dnn` estimates the conditional mean of `Y`, which differs from CATE:

```stata
np_dnn Y, xcols(x1 x2) x0matrix(X0)
matrix NP_summary = r(summary)
estimates store np_fit
matrix list NP_summary
np_dnn_diagnostics
np_dnn_plot, type(estimates) name(np_profile) replace
```

The optional `type(outcome)` plot shows weighted outcome quantiles, not
confidence bounds for the mean. Each fit replaces the active `e()` results;
restore the appropriate model before postestimation:

```stata
estimates restore cate_fit
cate_dnn
matrix list e(targets)
```

See [reusing fits](reusing-fits.md) for saved-model prediction, candidate
extraction and numeric reports. Saved fits retain observation-level data and
inherit the data's sharing restrictions.

## Read the localization and inference settings

Defaults use exact DNN, `k(1)` and L2 distance. The order `s` controls
localization: larger values put more weight close to the target. It does not
mean keeping the nearest `s` observations. `tuningmode(grid)` compares orders
while retaining the ordinary default result; it does not select a preferred
candidate. See [weights and tuning](weights-and-tuning.md) for the default
rule, explicit `s()` and candidate options. Jackknife uncertainty needs `s < n`;
`nose` skips it.

Intervals are pointwise and the sampling framework assumes independent
observations. Cluster and survey-design inference are not supplied. Row order
resolves equal-distance ties. Diagnostics describe support and concentration;
they cannot verify the paper's assumptions or remove localization bias.

For conditional LATE, continue with the [binary-IV guide](conditional-late.md).
[Notation](notation-and-conventions.md) explains shapes and result names;
the [article index](index.md) routes the remaining workflows.

## Migrating earlier CATE inputs

The former `mcol()`/`mhat()` construction is removed. Supply both arm-specific
outcome means and a propensity for AIPW, or a justified supplied signal.
Duplicating a pooled mean is not an equivalent migration. R demo columns
`m1_oracle` and `aipw_oracle` replace a different signal formerly called
`rho_oracle`.
