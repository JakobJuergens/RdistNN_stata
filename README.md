# RdistNN for Stata

`RdistNN_stata` provides distributional nearest-neighbor (DNN) estimation of
conditional average treatment effects and nonparametric regression, with
pointwise uncertainty, localization diagnostics, tuning grids, and Stata-native
plots. It is the Stata companion to
[RdistNN](https://github.com/JakobJuergens/RdistNN) and
[Orthogonal Inference for Conditional Z-Estimation without Cross-Fitting](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).

The initial public release focuses on CATE, regression, and their supporting
weights and postestimation tools. PLQR functionality is experimental: it
evaluates score moments at supplied coefficients and provides moment
covariance, without solving for coefficients or providing coefficient inference.

The package's original code is licensed under the
[GNU General Public License version 3 or later](LICENSE.md)
(`GPL-3.0-or-later`), with [additional permission to link the native plugin
with Stata](LINKING_EXCEPTION.md). Original third-party copyright notices are
preserved.

`net install` also installs the notices as `rdistnn_license.md` and
`rdistnn_linking_exception.md` in Stata's PLUS directory.

| Command | Purpose |
| --- | --- |
| `cate_dnn` | CATE from AIPW, Robinson, or precomputed scores |
| `cate_dnn_diagnostics`, `cate_dnn_plot` | CATE diagnostics and plots |
| `np_dnn` | Pointwise nonparametric regression |
| `np_dnn_diagnostics`, `np_dnn_plot` | Regression diagnostics and plots |
| `weights_dnn` | Compute and export exact or asymptotic DNN rank weights |
| `plqr_dnn` | Experimental PLQR score moments at supplied coefficients |

## Installation

Stata 17 or later is required. Install release **0.1.0**:

```stata
net install rdistnn, from("https://raw.githubusercontent.com/JakobJuergens/RdistNN_stata/v0.1.0/") replace
help cate_dnn
help np_dnn
```

Alternatively, download [rdistnn-source.zip](https://github.com/JakobJuergens/RdistNN_stata/releases/download/v0.1.0/rdistnn-source.zip)
and `SHA256SUMS-source.txt` from the [v0.1.0 release](https://github.com/JakobJuergens/RdistNN_stata/releases/tag/v0.1.0),
verify the checksum, extract the ZIP, and install from the directory containing
`rdistnn.pkg`:

```stata
net install rdistnn, from("path/to/extracted/rdistnn-source") replace
```

For the current development version, use the default branch:

```stata
net install rdistnn, from("https://raw.githubusercontent.com/JakobJuergens/RdistNN_stata/master/") replace
```

Stata's [net manual](https://www.stata.com/manuals/rnet.pdf) describes local
and remote installation. The ado/Mata implementation runs on Windows,
macOS, and Linux without a compiler after installation. Plugin support is
described below.

## Getting started

This simulated example uses **oracle nuisances**: the true propensity and
conditional outcome means that generate the data. In applied work, supply
appropriate fitted nuisance predictions as variables; the package does not
fit them. The complete runnable example is
[examples/getting-started.do](examples/getting-started.do).

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

* One target or a matrix whose rows are target profiles.
matrix X0 = (0,0 \ 0.5,0 \ -0.5,0)
cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) tuningmode(grid)
matrix list e(b)
matrix list e(V)
cate_dnn_diagnostics
cate_dnn_plot, type(estimates) replace
cate_dnn_plot, type(tuning) allcandidates replace
cate_dnn_plot, type(overlap) view(mass) replace

* NP regression uses the same target and tuning conventions.
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
np_dnn_diagnostics
np_dnn_plot, type(estimates) replace
np_dnn_plot, type(outcome) replace
```

Each estimator replaces Stata's active estimation results. Run its
postestimation commands after the corresponding fit. Batch `e(V)` retains
cross-target covariance; the displayed intervals are pointwise and do not
provide simultaneous coverage over a target grid. `nose` requests estimates
without jackknife uncertainty.

For a precomputed CATE score, use the positional signal with `pseudooutcome`:

```stata
gen double eclip = min(max(e_oracle,1e-8),1-1e-8)
gen double rho = m1_oracle-m0_oracle + D*(Y-m1_oracle)/eclip ///
    - (1-D)*(Y-m0_oracle)/(1-eclip)
cate_dnn rho, xcols(x1 x2) x0(0 0) pseudooutcome nose
```

For a Robinson score, provide `ecol()` and the pooled outcome mean in
`mcol()`; see `help cate_dnn`. The
[command guide](docs/command-guide.md) covers score routes, tuning,
diagnostics, plots, compatibility aliases, and experimental PLQR examples.

## Inference and localization

Neyman orthogonality controls first-order population nuisance drift. Valid
same-sample nuisance reuse additionally requires the paper's nuisance-rate
and full-state stability conditions, together with its overlap, localization,
and target-bias requirements. `cate_dnn` cannot establish those conditions
for an arbitrary learner. Supply oracle quantities, cross-fitted or out-of-fold
predictions, or full-sample predictions from learners for which the required
conditions are credible. Jackknife uncertainty does not itself remove
localization bias. The README example illustrates mechanics using known
nuisances; it is not empirical evidence about nuisance estimation.

The default `k(1)` selects one nearest observation within each size-`s`
subsample. The paper rule sets
`s = floor(.4*n^((kx+1)/(kx+2)))`, where `n` is the complete estimation-sample
size and `kx` is the covariate dimension. `tuningmode(grid)` evaluates
`lambda = (.3, .4, .5, .75)` while keeping the ordinary fit at the default
`.4` candidate. Use `lambda()`, `s()`, `lambdagrid()`, or `sgrid()` for
documented overrides.

For exact top-`k` weights, positive support has size `n - s + k`; under
`k(1)`, it is `n - s + 1`. The direction of `s` is opposite to the usual bandwidth convention:
larger s tightens localization, while smaller s spreads weights over more
observations. Distinguish effective sample size, weight concentration,
Gaussian inference, and target bias when interpreting diagnostics. Full
contracts and the fixed-metric Mahalanobis jackknife are in the command guide
and installed help.

## Optional native plugin

The bundled `src/rdistnn.plugin` is a Windows x86-64 PE binary. It accelerates
only standalone `weights_dnn` L2 calls. Estimators and Mahalanobis weights
use Mata. Eligible L2 calls attempt the plugin on every platform; an
incompatible or missing plugin falls back to Mata.

```stata
global rdistnn_use_plugin 0   // force the Mata reference path
global rdistnn_use_plugin 1   // require the plugin and surface failures
macro drop rdistnn_use_plugin // restore default dispatch
```

Published releases receive plugin assets for Windows x86-64, Linux x86-64,
macOS arm64, and macOS x86-64, plus `SHA256SUMS.txt`. Use the asset matching
your platform, rename it to `rdistnn.plugin`, and place it on the Stata
adopath. Native binaries are platform-specific. Build and parity instructions
are in [src/plugin/README.md](src/plugin/README.md).

## Development and release preparation

[CONTRIBUTING.md](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md) describes local setup and verification.
[The test-suite map](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/README.md) lists numerical, cross-language,
validation, and postestimation coverage. [RELEASING.md](https://github.com/JakobJuergens/RdistNN_stata/blob/master/RELEASING.md) records
the source archive, installation, plugin, and publication checks;
[WEBSITE.md](https://github.com/JakobJuergens/RdistNN_stata/blob/master/WEBSITE.md) supplies announcement copy. Release preparation does
not itself publish a release. See [RELEASE_NOTES.md](RELEASE_NOTES.md) for
the change history.
