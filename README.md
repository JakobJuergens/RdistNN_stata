# RdistNN_stata

`RdistNN_stata` is a Stata package for distributional nearest-neighbor (DNN)
local estimation and pointwise conditional Z-estimation. It is the Stata
companion to the R package
[RdistNN](https://github.com/JakobJuergens/RdistNN) and the paper
[Orthogonal Inference for Conditional Z-Estimation without Cross-Fitting](https://github.com/JakobJuergens/Cond_Z_Estimation).

The current package exposes four estimation/weight commands and paired
postestimation commands for NP regression and CATE:

| Command | Purpose |
|---|---|
| `np_dnn` | Pointwise nonparametric regression at a target covariate profile |
| `np_dnn_diagnostics` | Inspect or save NP target-level and tuning-candidate diagnostics |
| `np_dnn_plot` | Plot NP estimates, diagnostics, tuning, weights, outcomes, and dashboards |
| `cate_dnn` | Pointwise conditional average treatment effect from AIPW, Robinson, or precomputed scores |
| `plqr_dnn` | Experimental localized PLQR score moments at a supplied theta vector |
| `weights_dnn` | Utility for computing and exporting exact DNN rank weights |
| `cate_dnn_diagnostics` | Inspect or save target-level and tuning-candidate diagnostics |
| `cate_dnn_plot` | Plot estimates, diagnostics, tuning, overlap, weights, influence, and dashboards |

Stata 17 or later is required.

## Installation

Install the current default-branch package directly from GitHub:

```stata
net install rdistnn, from("https://raw.githubusercontent.com/JakobJuergens/RdistNN_stata/master/") replace
```

For local development, clone or copy this repository, then run
`dev-setup.do` once per Stata session to add `src/` to the adopath:

```stata
do "path\to\RdistNN_stata\dev-setup.do"
```

To make the adopath addition permanent across sessions, add the `adopath ++`
line to your `profile.do` (find it with `di c(sysdir_personal)`).

## Platform compatibility

The ado and Mata implementation is the portable baseline and is intended to
run under Stata 17 or later on Windows, macOS, and Linux. It does not require a
native compiler after installation.

The bundled `src/rdistnn.plugin` is currently a Windows x86-64 PE binary. It
accelerates only the standalone `weights_dnn` L2 path. On another operating
system, install the matching native binary to make the fast path available.
Eligible L2 calls attempt the plugin on every platform; load or execution
failure falls back to Mata. Mahalanobis weights and the estimator commands use
Mata regardless. Setting `global rdistnn_use_plugin 1` requires a working
native plugin and therefore fails instead of falling back. Setting it to `0`
forces Mata.

Native plugin sources and build instructions are under `src/plugin/`. The
plugin portability workflow compiles and inspects Windows x86-64, Linux
x86-64, macOS arm64, and macOS x86-64 binaries and publishes them as
finite-retention workflow artifacts. Native plugin binaries are
platform-specific and cannot be copied between those targets. Published
GitHub releases receive all four named plugin assets plus a SHA-256 checksum
file. After placing the matching asset on the Stata adopath as
`rdistnn.plugin`, default dispatch uses it when compatible; without a
compatible plugin, Mata remains available.

`src/plugin/` is the canonical native source directory. Build products from
that directory are ignored; the only plugin binary tracked for ordinary
`net install` is `src/rdistnn.plugin`.

## Quick start

```stata
* Simulate data
clear
set obs 1000
set seed 42
gen x1 = rnormal()
gen x2 = rnormal()
gen e  = invlogit(0.5*x1 - 0.5*x2)
gen D  = (runiform() < e)
gen m0 = 1 + x1 - x2
gen m1 = m0 + 0.5 + 0.5*sin(x1)
gen Y  = m0 + D*(m1 - m0) + rnormal()*0.5

* Nonparametric regression at one or several targets
np_dnn Y, xcols(x1 x2) x0(0 0)
matrix X0 = (0,0 \ 0.5,-0.25 \ -0.5,0.75)
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix list e(V)
matrix list e(tuning_estimates)
np_dnn_diagnostics
np_dnn_plot, type(dashboard)

* Estimate CATE at the origin with oracle nuisances.
* Orthogonality alone does not validate arbitrary same-sample nuisances;
* see the inference contract below.
cate_dnn Y, tcol(D) xcols(x1 x2) x0(0 0) ecol(e) m0col(m0) m1col(m1)

* Equivalent precomputed-score path; nose skips jackknife inference
gen eclip = min(max(e,1e-8),1-1e-8)
gen rho = D*(Y-m1)/eclip - (1-D)*(Y-m0)/(1-eclip) + m1-m0
cate_dnn rho, xcols(x1 x2) x0(0 0) pseudooutcome nose

* Robinson/R-learner score from a pooled outcome regression
gen mhat = e*m1 + (1-e)*m0
cate_dnn Y, tcol(D) xcols(x1 x2) x0(0 0) ecol(e) mcol(mhat)

* Evaluate several targets jointly; e(V) retains cross-target covariance
cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0) ///
    ecol(e) m0col(m0) m1col(m1)
matrix list e(targets)
matrix list e(V)

* Evaluate the R paper tuning grid; e(b)/e(V) remain at lambda(.4)
cate_dnn rho, xcols(x1 x2) x0matrix(X0) pseudooutcome tuningmode(grid)
matrix list e(tuning_candidates)
matrix list e(tuning_estimates)
matrix list e(tuning_diagnostics)
matrix list e(diagnostics)

* Postestimation diagnostics and Stata-native R-analogous plots
cate_dnn_diagnostics
cate_dnn_diagnostics, allcandidates
cate_dnn_plot, type(estimates)
cate_dnn_plot, type(diagnostics)
cate_dnn_plot, type(tuning) allcandidates
cate_dnn_plot, type(overlap) view(mass)
cate_dnn_plot, type(weights) view(kernel) targets(1 3)
cate_dnn_plot, type(influence) targets(1 3)
cate_dnn_plot, type(dashboard)

* Inspect DNN weights
weights_dnn x1 x2, x0(0 0) k(5) generate(w)
summarize w
```

## CATE inference contract

Neyman orthogonality controls first-order population nuisance drift. It does
not, by itself, validate arbitrary nuisance predictions estimated on the same
observations used by `cate_dnn`. Valid full-sample reuse additionally requires
the nuisance-rate and full-state stability conditions in the accompanying
paper.

`cate_dnn` computes the localized estimator and rank-shift variance procedure
from the supplied score or nuisance variables. It cannot diagnose whether an
arbitrary learner satisfies the paper's theoretical rate and stability
conditions. Users may supply:

- oracle quantities;
- cross-fitted or out-of-fold predictions; or
- full-sample predictions from learners for which the paper's nuisance-rate
  and full-state stability conditions are credible.

The Stata package deliberately does not fit, select, or prescribe nuisance
models. Nuisance construction remains modular and external to `cate_dnn`.

## Experimental PLQR moment evaluation

`plqr_dnn` is an experimental localized-moment evaluator, not a complete PLQR
coefficient estimator. It evaluates the orthogonal score at a supplied
`thetaeval()` using supplied `g()` and `md()` nuisance variables. The command
posts `e(moment)`, `e(moment_vcov)`, and `e(theta_eval)`; compatibility copies
remain in `e(b)`, `e(V)`, and `e(beta)` for Stata eclass workflows.

```stata
gen md1 = 0.3*x1 - 0.2*x2
gen md2 = -0.1*x1 + 0.4*x2
gen d1  = md1 + rnormal()*0.4
gen d2  = md2 + rnormal()*0.4
gen g   = 0.25*x1 - 0.15*x2
gen Yq  = d1 - 0.5*d2 + g + rnormal()*0.5
plqr_dnn Yq, x(x1 x2) x0(0 0) d(d1 d2) g(g) md(md1 md2) ///
    thetaeval(1 -0.5) tau(0.5)
```

The displayed result states that no local root solving was performed. The
reported covariance belongs to the localized score moment, not to theta.
Temporary `q()` and `beta()` aliases remain available but cannot be combined
with `g()` or `thetaeval()`. Local root solving and density-backed coefficient
inference are future functionality.

## Verification

Static package checks:

```powershell
bash tools/ci/pkg_check.sh
bash tools/ci/archive_check.sh "$env:TEMP/rdistnn-source.zip"
git diff --check
```

Compile and inspect the native plugin for the current host platform:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/plugin_build_check.ps1
```

Full Stata runtime tests on a machine with Stata 17 or later:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_runtime_check.ps1
```

The shared-contract crosswalk is in `tests/R_TEST_CROSSWALK.md`. Regenerate the
fixed R golden values used by `test-r-parity-priority.do` with:

```powershell
Rscript tools/generate_r_priority_fixtures.R
```

Fresh install smoke test, using isolated temporary `PLUS` and `PERSONAL`
directories:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_install_check.ps1
```

Both Stata helpers accept `-StataExe` to override executable discovery.  On the
current Windows development machine, StataNow 19 is discovered as
`C:\Program Files\StataNow19\StataSE-64.exe`.

The package-snapshot workflow builds and validates `source.zip`, then publishes
it as a GitHub Actions artifact with finite retention. Published releases also
receive the validated archive as `rdistnn-source.zip`. Generated archives are
not committed to the default branch. For local release preparation, build a
ZIP under `_dist/` or a temporary directory and pass that path to
`tools/ci/archive_check.sh`.

Failure logs and debug source copies are uploaded as finite-retention GitHub
Actions artifacts on the failed workflow run. They are no longer harvested
into committed `GH_logs/` directories. Open the workflow run's **Artifacts**
section to download them.

## Repository layout

```text
RdistNN_stata/
|-- src/             Core ado, Mata, help, and plugin files
|-- tests/           Stata-based regression and smoke tests
|-- tools/           Development helpers
|-- dev-setup.do     Adds src/ to the adopath for local development
|-- rdistnn.pkg      Package metadata for net install workflows
`-- stata.toc        Table of contents for Stata package distribution
```

## Notes

- `np_dnn` now follows the same multi-target and integrated-tuning contract as
  `cate_dnn`. Use `xcols()` and either `x0()` or `x0matrix()`. Batch `e(V)` is
  the joint rank-shift jackknife covariance, and `tuningmode(grid)` evaluates
  the R paper grid while reusing target ranks. `nose` skips inference and
  permits `s=n`.
- Every NP fit posts `e(diagnostics)` with the R package's weight and distance
  fields plus localized outcome mean, standard deviation, and weighted 10th,
  50th, and 90th percentiles. `np_dnn_diagnostics` displays or saves target and
  candidate tables; `np_dnn_plot` provides estimate, diagnostics, tuning,
  weight, outcome, and dashboard views.
- Deterministic fixtures enforce R/Stata NP parity for batch estimates,
  standard errors, ESS, all target diagnostics, and every paper-grid candidate.
- `cate_dnn` supports the same three core score paths as the R estimator:
  a precomputed positional signal with `pseudooutcome`, full AIPW nuisances
  with `ecol()`, `m0col()`, and `m1col()`, or a Robinson score with
  `ecol()` and `mcol()`. See the CATE inference contract above before using
  same-sample nuisance predictions; orthogonality alone is not sufficient.
  The default `k(1)` matches the paper theory and the R package. Specify
  `nose` for estimate-only fits.
  Use `x0()` for one target or `x0matrix()` with a named target matrix for
  batch evaluation. Batch `e(V)` is the joint rank-shift jackknife covariance;
  its diagonal agrees with separate R/Stata target fits.
- `cate_dnn, tuningmode(grid)` evaluates the paper sensitivity grid
  `lambda = (.3, .4, .5, .75)` in one call and reuses the target ranks.
  `lambdagrid()` customizes the multiplier grid. The Stata-native `sgrid()`
  supplies integer orders directly, with `sdefault()` selecting the default.
  Candidate estimates, marginal SEs, and diagnostics are returned in
  `e(tuning_estimates)`, `e(tuning_se)`, and `e(tuning_diagnostics)`;
  `e(tuning_candidates)` records feasibility and duplicate-computation IDs.
  The deterministic test fixture asserts all paper-grid values against the R
  package, not only qualitative invariants.
- Every CATE fit posts `e(diagnostics)` with the R diagnostic contract's
  numerical fields: fixed top-1/top-5/top-10 weight shares, nonzero support,
  nearest-rank radii, weighted distance summaries, and, for nuisance-based
  score paths, treatment-mass and propensity-overlap summaries. Deterministic
  R fixtures cover both distance and overlap diagnostics. Paper diagnostics
  additionally report the ordinary-Euclidean radius containing 90 percent of
  DNN weight, normalized within-arm treated/control ESS, and the top-five share
  of absolute weighted score contributions. These remain distinct from the
  rank-`s` radius, overall ESS, and top-five raw-weight share. Arm ESS requires
  nuisance-based treatment inputs; the physical radius and score-contribution
  fields are available for precomputed scores as well.
- Nuisance-based CATE paths first require every raw propensity to be finite and
  in `[0,1]`, then use `e_used = min(max(e, clip), 1-clip)` consistently in
  AIPW or Robinson scores, propensity diagnostics, and postestimation score
  reconstruction. `clip()` defaults to `1e-8`, must be finite and strictly
  between 0 and 0.5, and does not alter a precomputed pseudo-outcome.
- Because the estimator and tuning rules are deterministic, the regression
  suite enforces R/Stata parity for target estimates, covariance matrices,
  candidate grids, and the shared diagnostic fields.
- `cate_dnn_diagnostics` turns the stored target or tuning-candidate matrices
  into compact postestimation tables, returns the full table in
  `r(diagnostics)`, and can save it as a `.dta` file.
- `cate_dnn_plot` provides Stata-native analogues of the R package's CATE
  estimate, diagnostics, tuning, overlap, weight, influence, and dashboard
  plots. Stored-result views do not refit. Effective-kernel and influence
  views reconstruct weights only after validating the current estimation
  sample and referenced variables. All graph builders support named graphs,
  `nodraw`, and export through `saving()`.
- Stata option names cannot contain underscores, so the CATE interface follows
  the R names by removing underscores only: `x_cols` -> `xcols()`, `t_col` ->
  `tcol()`, `e_col` -> `ecol()`, `m0_col` -> `m0col()`, `m1_col` -> `m1col()`,
  `m_col` -> `mcol()`, and `pseudo_outcome` -> `pseudooutcome`. R's `X0`
  matrix is passed with `x0matrix()`, `tuning_mode` with `tuningmode()`, and
  `lambda_grid` with `lambdagrid()`. Earlier Stata names (`x()`, `treat()`,
  `propensity()`, `m0()`, `m1()`, and `mhat()`) remain compatibility aliases.
- The same naming rule applies to NP regression: R's `x_cols`, `X0`,
  `tuning_mode`, and `lambda_grid` are Stata's `xcols()`, `x0matrix()`,
  `tuningmode()`, and `lambdagrid()`. The earlier NP `x()` option remains a
  compatibility alias.
- Experimental `plqr_dnn` requires user-supplied `g()` and `md()` nuisance
  variables and a supplied `thetaeval()` vector. It reports localized score
  moments and their jackknife covariance; it does not solve for theta, estimate
  nuisances, perform GMM, or provide coefficient inference. `q()` and `beta()`
  are temporary compatibility aliases.
- `weights_dnn` uses the installed L2 plugin path by default when the plugin is
  available. Set global `rdistnn_use_plugin` to `0` to force the Mata reference
  path, or to `1` to require the plugin path and fail on plugin errors.
- `weights_dnn` defaults to `k(1)`. Because it constructs weights without a
  jackknife, it accepts the full mathematical endpoint `k <= s <= n`, including
  `s(n)`. Exact/asymptotic methods, L2/Mahalanobis distance, generated weights,
  plugin fallback, and returned `r()` metadata use the same interface.
- `np_dnn, nose` skips the rank-shift jackknife, permits `s=n`, omits `e(V)`,
  and records `e(se_method)="none"` and `e(jackknife)="none"`. Initial display
  and replay state that inference was not requested. Without `nose`, `s<n`
  remains required.
- `k` is the number of nearest observations selected within each size-`s`
  subsample; it is not the number of observations with positive weight in the
  full-sample rank-weight vector. For exact top-`k` DNN weights, positive
  support ends at global rank `n - s + k`, so the support size is
  `n - s + k`. Under the paper-default `k(1)`, this becomes `n - s + 1`.
- `s` is the DNN subset/kernel order and a localization or bandwidth analogue,
  but its direction is opposite to the usual bandwidth convention. Holding
  `n` and `k` fixed, larger `s` shrinks positive-weight support and produces
  tighter localization; smaller `s` spreads weight across more observations
  and increases the effective local sample size. The finite-sample default is
  `lambda(.4)`, resolving `s = floor(.4*n^((kx+1)/(kx+2)))`; `lambda(1)` is a
  more aggressive localization benchmark because it yields a larger `s`.
- Use `distance(mahalanobis)` when covariates are strongly correlated. The
  jackknife keeps the full-sample covariance metric fixed across leave-one-out
  deletions. All public commands require at least two complete estimation
  observations for Mahalanobis distance; rank-deficient covariance matrices
  with larger samples retain generalized-inverse behavior.
- Standard errors are computed with the fast rank-shift jackknife, exploiting
  the L-statistic structure of the DNN estimator.

## Related repositories

- [RdistNN](https://github.com/JakobJuergens/RdistNN): R package implementing
  the broader DNN weighting and local-estimation framework.
- [Cond_Z_Estimation](https://github.com/JakobJuergens/Cond_Z_Estimation):
  paper repository for the conditional Z-estimation theory motivating these
  commands.
- [Cond_GMM_with_DNN](https://github.com/JakobJuergens/Cond_GMM_with_DNN):
  companion conditional GMM paper that builds on the same localization layer.
- [JK_GenUStat](https://github.com/JakobJuergens/JK_GenUStat): supporting
  generalized U-statistic and jackknife variance theory.
