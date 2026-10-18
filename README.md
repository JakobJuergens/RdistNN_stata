# RdistNN for Stata

**RdistNN for Stata** estimates conditional treatment effects and regression
functions using distributional nearest-neighbor (DNN) weights. Fit at selected
covariate profiles, report pointwise estimates and delete-one jackknife
standard errors, and inspect localization, overlap and sensitivity to tuning
choices.

CATE and nonparametric regression are the main workflows. The package also
provides **experimental conditional LATE (CLATE)** with shared reporting,
diagnostics and Stata-native plots. You supply nuisance predictions for the
causal estimators; the package constructs the adjusted signals and performs
the local estimation and uncertainty calculations.

[Get started](docs/getting-started.md)
· [Command guide](docs/command-guide.md)
· [CLATE guide](docs/conditional-late.md)
· [Report an issue](https://github.com/JakobJuergens/RdistNN_stata/issues)

## Installation

Install from GitHub in Stata:

```stata
net install rdistnn, from("https://raw.githubusercontent.com/JakobJuergens/RdistNN_stata/master/") replace
help rdistnn
```

Stata **17 or later** is required. Ordinary use needs no R installation,
Python or compiler. The ado/Mata implementation supports Windows, macOS and
Linux; the optional native plugin is described below.

For a version-pinned installation, use a tag in place of `master`, for example
`v0.1.0`. [Release 0.1.0](https://github.com/JakobJuergens/RdistNN_stata/releases/tag/v0.1.0)
contains the original CATE and regression workflows and predates the recent
CLATE and shared reporting additions; open `help cate_dnn` or `help np_dnn`
with that release. Release assets also provide a source ZIP, checksums and
platform-specific plugins. See the
[installation guide](docs/getting-started.md#install-and-open-help) for alternatives.
For installation from a downloaded source ZIP, verify its checksum, extract
it, and run `net install rdistnn, from("path/to/extracted/rdistnn-source") replace`,
using the directory containing `rdistnn.pkg`.

## Choose a workflow

| Goal | Command | Inputs and scope |
| --- | --- | --- |
| Conditional average treatment effects | `cate_dnn` | Outcome, binary treatment and supplied outcome/propensity predictions, or a supplied pseudo-outcome |
| Conditional mean regression | `np_dnn` | Outcome and covariates; no nuisance predictions needed |
| Conditional effects among compliers | `clate_dnn` | Binary instrument and treatment, outcome and supplied instrument-arm nuisance predictions; **experimental** |
| A ratio of two local conditional means | `clate_dnn, pseudooutcome denominator()` | Two aligned signals with paired covariance; **experimental** |
| Inspect or reuse DNN weights | `weights_dnn` | Covariates and a target; exports exact or asymptotic rank weights |

Use `xcols()` to select localization covariates, `x0()` for one target, or
`x0matrix()` for a matrix with one target per row. CATE, CLATE and regression
support exact weights (the default), `method(asymp)`, L2/Mahalanobis distance,
and candidate tuning grids. For supplied-coefficient quantile-score moments,
see the [experimental PLQR guide](docs/plqr-moments.md); `plqr_dnn` does not
solve for coefficients or provide coefficient inference.

## Conditional average treatment effects

This self-contained example replaces the active dataset with simulated
observations. It uses **oracle nuisance values** to illustrate the interface.
In an application, replace these columns with aligned predictions of the
two conditional outcome means and treatment propensity. The estimator does
not fit those nuisance models.

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

matrix X0 = (-1,0 \ -0.5,0 \ 0,0 \ 0.5,0 \ 1,0)
cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) tuningmode(grid)
matrix cate_summary = r(summary)
matrix list cate_summary
cate_dnn_diagnostics
cate_dnn_plot, type(estimates) style(profile) replace
```

This estimates AIPW CATE at five profiles, holding `x2` at zero. For a
Robinson score, supply `ecol()` and the pooled outcome mean in `mcol()`.
For an externally constructed score, use
`cate_dnn score, xcols(x1 x2) x0(0 0) pseudooutcome`.
The [CATE and regression guide](docs/cate-and-regression.md) explains the
input routes and returned results; the
[complete quick start](examples/getting-started.do) runs each route.

## Nonparametric regression

Using the same simulated outcome and targets:

```stata
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix regression_summary = r(summary)
matrix list regression_summary
np_dnn_diagnostics
np_dnn_plot, type(estimates) replace
np_dnn_plot, type(outcome) replace
```

The estimate is the local conditional mean of `Y`. Outcome panels additionally
show weighted outcome quantiles: their 10th–90th percentile range describes
local outcome dispersion, rather than uncertainty about the mean.
Each estimator replaces Stata's active estimation results; use its matching
diagnostics and plot commands after fitting it.

## Conditional LATE (experimental)

With a binary instrument and binary treatment, `clate_dnn` estimates the
conditional effect among compliers under the IV assumptions. Supply the
instrument propensity and the conditional means of the outcome and treatment
in each **instrument arm**. The reporting coordinates in `xcols()` can be a
subset of the covariates used to construct those predictions.

The following example is self-contained and replaces the active dataset.
It uses a simulated one-sided encouragement design with labelled oracle
predictions. `A` is the instrument and `D` is the treatment.

```stata
clear
set obs 160
set seed 4711
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double e_oracle = invlogit(.3*x1-.2*x2)
gen byte A = runiform() < e_oracle
gen double my0_oracle = .5*x1-.2*x2
gen double md1_oracle = invlogit(.4+.2*x1)
gen double md0_oracle = 0
gen double tau_oracle = 1+.3*sin(x1)
gen byte D = A*(runiform() < md1_oracle)
gen double Y = my0_oracle + tau_oracle*D + rnormal()*.3
gen double my1_oracle = my0_oracle + md1_oracle*tau_oracle

matrix X0 = (-.5,0 \ 0,0 \ .5,0)
clate_dnn Y, xcols(x1 x2) x0matrix(X0) dcol(D) acol(A) ///
    ecol(e_oracle) my0col(my0_oracle) my1col(my1_oracle) ///
    md0col(md0_oracle) md1col(md1_oracle) onesided ///
    training(oracle) tuningmode(grid)
matrix clate_summary = r(summary)
matrix list clate_summary
matrix list e(components)
matrix list e(status)
```

Both component means use the same DNN weights. Ratio uncertainty retains
their paired covariance, holding the supplied signals fixed. To supply the
signals directly, use the positional numerator with `pseudooutcome` and
`denominator(variable)`. The [CLATE guide](docs/conditional-late.md) and
[runnable IV example](examples/clate.do) explain both routes and numerical
statuses. Inspect the first stage alongside the ratio: ordinary Wald
intervals can be misleading when the first stage is weak. Undefined ratios
retain their components and original target positions in the full tables.

## Shared reporting and plots

CATE, CLATE and regression follow the same reporting pattern. After the
CLATE fit above, replay at a different confidence level without refitting:

```stata
clate_dnn, level(90)
matrix coefficient_table = r(table)
matrix summary_table = r(summary)
clate_dnn_diagnostics, targets(1 3) detail
matrix diagnostic_table = r(diagnostics)
clate_dnn_plot, type(estimates) style(forest) level(90) replace
matrix figure_data = r(plot_data)
```

Use `cate_dnn` or `np_dnn` and its matching postestimation commands for the
other workflows. `r(summary)` has columns `estimate`, `std_error`, `z`,
`pvalue`, `ci_lower` and `ci_upper`. Copy `r()` results before another
r-class command replaces them. `e(estimates)`, `e(se_targets)` and
`e(estimate_vcov)` retain all targets and cross-target covariance; CLATE's
`e(b)` may contain only the subset with defined ratios.

| Plot type | What it shows |
| --- | --- |
| `type(estimates)` | Estimates and pointwise confidence bounds; profile, forest or pointrange styles |
| `type(weights)` | Rank weights, cumulative mass or observation-level kernel weights |
| `type(overlap)` | Treatment overlap for CATE; **instrument** overlap for CLATE |
| `type(influence)` | CATE score contributions or CLATE paired-score residual contributions |
| `type(tuning)` | Sensitivity across stored tuning candidates |
| `type(diagnostics)` | Localization and weight concentration |
| `type(outcome)` | Conditional means and weighted outcome quantiles; regression only |
| `type(dashboard)` | A combined view of the available panels |

Overlap and influence panels apply to CATE and CLATE; outcome panels apply
to regression.

Add `saving("figure.png")` or `saving("figure.pdf")` to export a graph;
diagnostic commands accept `saving("diagnostics.dta")`. Use `targets()` and
`candidates()` to inspect selected profiles and tuning choices. Observation-
level plots require the original estimation data and order; refit after
changing them. See [diagnostics and reporting](docs/diagnostics-and-reporting.md)
and [visualising fits](docs/visualising-fits.md) for options and returned data.

![CATE estimates with pointwise confidence bounds from the simulated diagnostics tutorial](docs/figures/cate-profile.png)

This profile comes from the [diagnostics and plots tutorial](examples/diagnostics-and-plots.do),
which uses simulated observations and oracle nuisances. The ribbon shows
pointwise confidence intervals. The [plot gallery](docs/visualising-fits.md#gallery-from-the-runnable-example)
also illustrates dashboards, CLATE forest plots and local outcome distributions.

## Localization and tuning

With default `k(1)` and no explicit order, the paper rule uses
`s = floor(.4*n^((kx+1)/(kx+2)))`, where `n` is the complete estimation-sample
size and `kx` is the number of localization covariates. `s` is the subsample
order: larger `s` tightens localization; this is
opposite to the usual bandwidth convention.
Exact top-`k` weights have positive support of size `n - s + k`; with
`k(1)`, this is `n - s + 1`.

Supply `s()` directly or use `tuningmode(grid)` to compare multipliers
`.3 .4 .5 .75`, retaining `.4` as the fitted default. `sgrid()` supplies
orders directly. Grids are sensitivity analyses; they do not automatically
select a tuning value. Plot stored candidates after a grid fit with
`clate_dnn_plot, type(tuning) allcandidates replace`, or the corresponding
CATE/regression command. See [weights and tuning](docs/weights-and-tuning.md).

Jackknife uncertainty requires `s < n`. `nose` skips uncertainty calculations.
The manuscript's DNN derivations use exact weights and `k=1`; larger fixed
`k` and asymptotic weights need their corresponding statistical conditions.

## Inference and localization

Reported confidence intervals are **pointwise**, rather than simultaneous
bands across a covariate profile. Valid inference requires the relevant
design, overlap, nuisance accuracy, deletion stability and localization-bias
conditions. Orthogonality alone does not justify arbitrary same-sample
nuisance reuse; the paper's nuisance-rate and full-state stability conditions
also matter. Oracle, out-of-fold and full-sample predictions each require
the assumptions appropriate to their training procedure.

Supplied propensity predictions must be finite and lie in `[0,1]`; nuisance-
based CATE/CLATE scores clip them to `[clip,1-clip]`, with default `clip=1e-8`.
Supplied pseudo-outcomes are used unchanged. Clipping and favorable effective
sample size or overlap diagnostics do not establish valid inference, and
jackknife uncertainty does not remove localization bias.

CLATE additionally needs consistency, exclusion, conditional instrument
independence, monotonicity and positive local compliance. Its target generally
differs from population CATE. The Stata interface provides pointwise Wald
reporting, without a weak-instrument-uniform guarantee. See the
estimator-specific inference sections in the [CATE/regression](docs/cate-and-regression.md)
and [CLATE](docs/conditional-late.md) guides before using these intervals.

## Optional native acceleration

The optional plugin accelerates standalone `weights_dnn` calls with L2
distance and frozen-score jackknife deletions for NP, CATE, PLQR and CLATE.
Estimator geometry and weight tables use Mata, including Mahalanobis geometry.
Both compiled and Mata jackknife engines use prefix/suffix sums with one
rank-weight table per order, avoiding repeated delete-one weight calculations.
The bundled plugin is a Windows x86-64 PE binary. Eligible calls attempt the plugin on every platform;
execution falls back to Mata if it is missing or incompatible.

Jackknife dispatch tries compiled prefix/suffix sums, then the fast Mata
implementation, then retained naive deletions if the fast Mata engine is
unavailable or its sums are nonfinite. Set `global rdistnn_jackknife_engine`
to `plugin`, `mata` or `naive` for explicit route/parity checks; drop the global
to restore automatic dispatch. Native jackknife bridge matrices are created
through Mata and can exceed `c(matsize)` without changing that setting.
See the plugin guide for controls.

Release assets include plugins for Windows x86-64, Linux x86-64, macOS arm64
and macOS x86-64, with checksums. To use native acceleration, install the asset
matching your platform as `rdistnn.plugin` on Stata's adopath. See the
[plugin guide](src/plugin/README.md) for installation, build and fallback controls.

## Further documentation

- [Getting started](docs/getting-started.md) covers installation and a first result.
- [Notation and conventions](docs/notation-and-conventions.md) explains arguments and returned quantities.
- [Diagnostics and reporting](docs/diagnostics-and-reporting.md) covers local support, summaries and saved tables.
- [The article collection](docs/index.md) links all guides and six standalone examples; `help rdistnn` opens the installed command references.
- [The R documentation crosswalk](docs/r-documentation-crosswalk.md) maps shared workflows and R-only features. [RdistNN](https://github.com/JakobJuergens/RdistNN) is the numerical reference for shared implementations.

The package accompanies
[Orthogonal Inference for Conditional Z-Estimation without Cross-Fitting](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
See the [R companion's citation information](https://github.com/JakobJuergens/RdistNN/blob/master/inst/CITATION)
for the methodological references. When reporting results, identify the Stata
package version or commit you used. [Release notes](RELEASE_NOTES.md) record
changes and experimental features.

## Help and contributions

Report problems or ask usage questions through the
[issue tracker](https://github.com/JakobJuergens/RdistNN_stata/issues), including
a minimal reproducible example, your Stata version, operating system and the
installed release or commit. Contributions are welcome;
[CONTRIBUTING.md](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md)
links development and verification procedures.

The original code is licensed under
[GNU GPL version 3 or later](LICENSE.md), with
[additional permission to link the native plugin with Stata](LINKING_EXCEPTION.md).
Original third-party notices are preserved. Installation also places the
license notices in Stata's PLUS directory.
