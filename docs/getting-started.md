# Getting started

RdistNN for Stata estimates conditional average treatment effects (CATE) and
conditional means with distributional nearest-neighbor (DNN) weights. You
supply the observations, localization covariates and target profiles. For
CATE, you also supply nuisance predictions or a previously constructed
pseudo-outcome. The package computes local estimates, delete-1 jackknife
uncertainty, diagnostics and Stata graphs.

These articles document package **0.1.0.9000** on `master`, including experimental
conditional LATE and shared reporting. Start with the
[article index](index.md) to find a particular workflow.

## Install and open help

Stata 17 or later is required. Install from GitHub:

```stata
net install rdistnn, from("https://raw.githubusercontent.com/JakobJuergens/RdistNN_stata/master/") replace
help rdistnn
help cate_dnn
help np_dnn
```

For a version-pinned installation, replace `master` in the address with a
published tag. The historical `v0.1.0` release includes CATE, regression,
weights and experimental PLQR moments; it predates CLATE and the shared
reporting updates. The [README](../README.md#installation)
also describes installation from a downloaded source archive. Installation
does not require R or a compiler. Estimation uses the portable ado/Mata
implementation; the optional native plugin accelerates eligible standalone
`weights_dnn` calls. See [weights and tuning](weights-and-tuning.md) for the
plugin and fallback controls.

The examples below assume the commands are installed. From a source checkout,
`do dev-setup.do` at the repository root adds `src/` to the session's adopath.
The complete onboarding script is
[examples/getting-started.do](../examples/getting-started.do); a second
[CATE and regression script](../examples/cate-and-regression.do) covers all
three CATE score routes and result extraction. Both clear the active dataset.

## Generate a small example

This simulation supplies **oracle nuisances**: the propensity and conditional
outcome means used to generate the observations. It demonstrates the command
interface without making claims about a fitted nuisance learner.

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

In an application, replace the oracle columns with nuisance predictions from
your own estimation procedure. Keep each prediction on the row of the
observation it describes. The localization variables in `xcols()` may be a
subset of the adjustment variables used to fit those predictions.

## Estimate CATE at selected profiles

Use `x0()` for one profile, or `x0matrix()` for a matrix whose rows are
profiles. Target columns always follow the order in `xcols()`; matrix column
names do not rearrange the supplied values.

```stata
matrix X0 = (0,0 \ 0.5,0 \ -0.5,0)
cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) tuningmode(grid)
matrix CATE_summary = r(summary)
matrix CATE_b = e(b)
matrix CATE_V = e(V)
estimates store cate_fit
matrix list CATE_summary
```

This call constructs the augmented inverse-probability weighted (AIPW)
pseudo-outcome, then averages it with DNN weights around each profile.
`e(b)` is a 1-by-3 estimate vector. `e(V)` is the 3-by-3 joint jackknife
covariance, retaining dependence across targets. The displayed normal
intervals have **pointwise** coverage; a profile graph does not turn them into
a simultaneous confidence band.

`r(summary)` contains one row per target and columns `estimate std_error z
pvalue ci_lower ci_upper`. It is temporary: copy it before another command
overwrites `r()`. Replay can change the confidence level without fitting
again:

```stata
cate_dnn, level(90)
matrix CATE_summary90 = r(summary)
cate_dnn_diagnostics, detail
cate_dnn_plot, type(estimates) name(cate_profile) replace
cate_dnn_plot, type(overlap) view(mass) name(cate_overlap) replace
cate_dnn_plot, type(tuning) allcandidates name(cate_tuning) replace
```

Diagnostics report concentration, distance and local treatment support.
Clipping an extreme propensity does not create treated or control
observations near a target. Diagnostic tables describe the realized sample;
they do not establish the assumptions required for causal interpretation or
inference. Details on retaining, selecting and saving results are in
[diagnostics and reporting](diagnostics-and-reporting.md), and graph options
are in [visualising fits](visualising-fits.md).

## Choose the CATE score input

The AIPW call above is one of three supported routes. If your workflow
supplies a pooled conditional outcome mean `m(X) = E[Y|X]`, use the Robinson
route:

```stata
gen double m_oracle = e_oracle*m1_oracle + (1-e_oracle)*m0_oracle
cate_dnn Y, tcol(D) xcols(x1 x2) x0(0 0) ///
    ecol(e_oracle) mcol(m_oracle)
```

Its signal is `(D-eclip)*(Y-m)/(eclip*(1-eclip))`; the command takes a
DNN-weighted mean of that signal. It does not solve a residual least-squares
R-learner problem. With a precomputed CATE pseudo-outcome, supply the signal
as the positional variable and add `pseudooutcome`:

```stata
gen double eclip = min(max(e_oracle,1e-8),1-1e-8)
gen double rho = m1_oracle-m0_oracle + D*(Y-m1_oracle)/eclip ///
    - (1-D)*(Y-m0_oracle)/(1-eclip)
cate_dnn rho, xcols(x1 x2) x0(0 0) pseudooutcome nose
```

`nose` requests an estimate without jackknife uncertainty. A supplied score
is used unchanged; the package cannot recover its nuisance inputs or show
treatment-overlap diagnostics from it. The full input and missingness
contracts are in [CATE and regression](cate-and-regression.md).

## Fit a conditional mean

`np_dnn` applies the same localization to the observed outcome. It estimates
`E[Y|X=x0]`, which differs from CATE.

```stata
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix NP_summary = r(summary)
estimates store np_fit
np_dnn_diagnostics, detail
np_dnn_plot, type(estimates) name(np_profile) replace
np_dnn_plot, type(outcome) name(np_outcomes) replace
```

The outcome view describes weighted outcome spread and quantiles. Those
quantiles are not confidence bounds for the conditional mean. The estimate
view uses jackknife uncertainty for the mean and still needs the regression
assumptions and adequate localization-bias control.

Every fit replaces the active `e()` results. Restore the appropriate fit
before using estimator-specific postestimation:

```stata
estimates restore cate_fit
cate_dnn
matrix list e(targets)
matrix list e(se_targets)
```

To evaluate new target profiles, issue another estimation call with those
targets. There is no separate Stata `predict` command for these estimators.
Stored-result summaries and graphs can use saved estimation results;
observation-level kernel and influence graphs additionally require the
original estimation data, values and observation order.

## Read the localization and inference settings

Defaults are `k(1)`, `method(exact)` and `distance(l2)`. The paper rule sets
`s = floor(.4*n^((kx+1)/(kx+2)))`, using the complete estimation-sample size
`n` and number of localization covariates `kx`. `k` counts nearest
observations inside each size-`s` subsample. Exact weights can extend through
full-sample rank `n-s+k`; `s` is not a request to keep the nearest `s` rows.
Larger `s` tightens localization. The paper's stated theory uses `k=1`.

`tuningmode(grid)` evaluates the default `.3 .4 .5 .75` multiplier grid
while preserving the ordinary `.4` fit as the main result. It is a
sensitivity calculation, not automatic selection of a preferred candidate.
Use `lambda()`, explicit `s()`, or a custom grid as explained in
[weights and tuning](weights-and-tuning.md).

For CATE, causal interpretation requires the identifying assumptions,
including appropriate treatment adjustment and overlap. Neyman
orthogonality controls first-order population nuisance drift; valid
same-sample nuisance reuse additionally requires the paper's nuisance-rate
and full-state stability conditions. Oracle, out-of-fold and full-sample
predictions are all accepted inputs. Their acceptance does not certify the
conditions for local inference. The jackknife holds the supplied scores
fixed and does not fit or refit nuisance learners, remove localization bias,
or check those assumptions.

Continue with [CATE and regression](cate-and-regression.md) for full score,
sample and covariance conventions. The [notation guide](notation-and-conventions.md)
connects the commands to the paper, while [conditional LATE](conditional-late.md)
and [PLQR moments](plqr-moments.md) describe the experimental workflows.
The [R documentation crosswalk](r-documentation-crosswalk.md) identifies
features available only in the R package.
