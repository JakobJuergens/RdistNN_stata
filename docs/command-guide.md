# Command guide

Start with [the README](../README.md) and the runnable
[getting-started example](../examples/getting-started.do). Stata's
`help COMMAND` is the reference for syntax, stored results, and validation.

## Workflow support

| Workflow | Estimation and results | Diagnostics and plots |
| --- | --- | --- |
| CATE | `cate_dnn`, one or several targets, joint jackknife covariance, tuning grids, `nose` | `cate_dnn_diagnostics`, `cate_dnn_plot` |
| Experimental CLATE (feature branch) | `clate_dnn`, augmented binary-IV or paired signals, paired joint covariance, batch/grid controls | `clate_dnn_diagnostics`, `clate_dnn_plot` |
| Nonparametric regression | `np_dnn`, one or several targets, joint jackknife covariance, tuning grids, `nose` | `np_dnn_diagnostics`, `np_dnn_plot` |
| Standalone weights | `weights_dnn`, exact/asymptotic L2 or Mahalanobis weights, generated variable and `r()` metadata | Inspect or export the generated weights |
| Experimental PLQR | `plqr_dnn`, one supplied target and coefficient vector, score moments and moment covariance | Explicit moment display and stored matrices |

The Stata companion has no TDNN or general local GMM command, automatic
nuisance fitting, or PLQR coefficient inference. The R package has a broader
interface; shared contracts and differences are recorded in the
[R-to-Stata test crosswalk](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/R_TEST_CROSSWALK.md).

## Conditional LATE (development branch)

Install `feature/clate` explicitly; CLATE is not in the v0.1.0 release.
The runnable [CLATE example](../examples/clate.do) supplies oracle nuisance
predictions for a binary instrument A and binary treatment D:

```stata
clate_dnn Y, xcols(x1 x2) x0matrix(X0) dcol(D) acol(A) ecol(e) ///
    my0col(my0) my1col(my1) md1col(md1) onesided tuningmode(grid)
clate_dnn_diagnostics, detail
clate_dnn_plot, type(overlap) view(mass) replace
clate_dnn_plot, type(influence) targets(1 3) replace
clate_dnn_plot, type(tuning) allcandidates replace
```

`onesided` requires D<=A and treats an omitted zero-arm treatment prediction
as zero. Without it, supply `md0col()`. The four arm means condition on
instrument assignment; their predictions may use all adjustment covariates
even when localization reports effects using fewer X coordinates.
The package fits no nuisance learners. `training()` is descriptive.

Alternatively, use `clate_dnn phi_y, denominator(phi_d) pseudooutcome ...`.
Both routes localize the two signals with identical weights and divide their
averages. Do not average individual ratios. All CATE weighting/target/tuning
controls and plot types are available. Overlap panels describe the instrument;
influence is `w*(phi_y-theta*phi_d)`. A paired-signal fit has no instrument
overlap information. Diagnostics retain paired components, first-stage
uncertainty and computational status.

`e(estimates)` and `e(se_targets)` contain every target, including missing
values. `e(components)` retains reduced-form and first-stage means;
`e(theta_vcov)` retains full cross-target ratio covariance and
`e(component_vcov)` retains interleaved Y,D component covariance. `e(status)`
flags zero, negative, near-zero and numerically unavailable cases. Negative
denominators are computed; the diagnostic tolerance never regularizes a ratio.
Since Stata rejects missing posted coefficients, `e(b)` contains only defined
targets, with their original target numbers; it is absent if none are defined.
`e(V)` is posted when every posted coefficient has available covariance.

Causal interpretation additionally requires instrument independence given the
full adjustment variables, consistency, exclusion, monotonicity, overlap and
positive local compliance. Full-sample inference requires the affine theorem's
nuisance accuracy, deletion stability, design, oracle and undersmoothing
conditions. Paired frozen-score deletions do not refit nuisances, eliminate
localization bias or provide weak-instrument-uniform confidence guarantees.
The [implementation design](superpowers/specs/2026-10-17-clate-design.md)
records the numeric contract and acceptance scope. RdistNN is the reference;
shared fixed fixtures test estimates, covariance and diagnostics.

## Plotting after a nuisance-based CATE grid fit

On the simulated data from the getting-started example, fit the nuisance
route before asking for treatment-overlap views:

```stata
cate_dnn Y, tcol(D) xcols(x1 x2) x0matrix(X0) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) tuningmode(grid)
cate_dnn_diagnostics, allcandidates
cate_dnn_plot, type(tuning) allcandidates replace
cate_dnn_plot, type(overlap) view(mass) replace
cate_dnn_plot, type(weights) view(kernel) targets(1 3) replace
cate_dnn_plot, type(influence) targets(1 3) replace
```

Overlap panels need the treatment and propensity inputs from an AIPW or
Robinson fit. A precomputed-score fit does not supply those inputs.
NP plots additionally provide `type(outcome)`; its weighted outcome
quantiles describe the local outcome distribution, rather than confidence
bounds for the conditional mean. NP and CATE graph commands accept
`name()`, `replace`, `nodraw`, and `saving()`. See their installed help for
view-specific options and data requirements.

## Localization and command contracts

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
  `ecol()` and `mcol()`. See [the CATE inference contract](../README.md#inference-and-localization) before using
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

## Experimental PLQR moment evaluation

The example below continues from the getting-started dataset. Its nuisance
variables and coefficient vector are oracle inputs used to illustrate the
score formula.

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
