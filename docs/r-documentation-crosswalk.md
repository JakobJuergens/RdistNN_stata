# R documentation coverage and Stata equivalents

This map compares Stata package **0.1.0.9000** on `master` with the RdistNN
source at `6fe48273aaa8853c5736841c5a228f8c1592ded9`. The reference inventory is
R's [article source directory](https://github.com/JakobJuergens/RdistNN/tree/6fe48273aaa8853c5736841c5a228f8c1592ded9/vignettes)
and [pkgdown reference index](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/_pkgdown.yml).
It covers all 13 articles and 14 reference families in that snapshot. R remains
the numerical reference for shared functionality. Documentation parity here
means usable coverage of Stata's full public surface and clear routing for
R-only functionality; it does not add R APIs to Stata.

## Articles

| R article source | Stata coverage | Boundary |
|---|---|---|
| `getting-started.Rmd` | [Getting started](getting-started.md), [CATE and regression](cate-and-regression.md) | All three CATE routes and NP; two-scale regression remains R-only |
| `notation-and-conventions.Rmd` | [Notation and conventions](notation-and-conventions.md) | Dataset columns/target matrices and e()/r() replace R input containers/objects |
| `diagnostics-reporting.Rmd` | [Diagnostics and reporting](diagnostics-and-reporting.md) | Shared localization and estimator-specific tables; PLQR conditioning diagnostics remain R-only |
| `visualising-fits.Rmd` | [Visualising fits](visualising-fits.md) | Applicable Stata-native views and returned plot data; ggplot composition remains R-only |
| `kernel-order-tuning.Rmd` | [Weights and tuning](weights-and-tuning.md#integrated-candidate-grids) | Paper/default/direct-order grids; arbitrary tuning specifications/registries remain R-only |
| `dnn-weights.Rmd` | [Weights and tuning](weights-and-tuning.md) | Exact/asymptotic rules and one-target weight variables; public batch weight/rank/order/distance objects remain R-only |
| `conditional-z-estimation.Rmd` | [CATE and regression](cate-and-regression.md), [PLQR moments](plqr-moments.md), [jackknife discussion](weights-and-tuning.md#jackknife-and-the-localized-target) | Implemented moments and inference conditions, not a generic Z/GMM root solver |
| `jackknife-gen-ustats.Rmd` | [Jackknife and localized targets](weights-and-tuning.md#jackknife-and-the-localized-target), [numerical tests](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/R_TEST_CROSSWALK.md) | Stata uses fast rank-shift covariance internally; the generic compact-support L-statistic utilities and selectable fast/slow/LOO APIs remain R-only |
| `extending.Rmd` | [Contribution workflow](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md), [native plugin](../src/plugin/README.md) | No public Stata score, distance, rank-weight or tuning registry; private Mata helpers are not a supported extension API |
| `conditional-late.Rmd` | [Conditional LATE](conditional-late.md) | Binary-IV/supplied-signal routes, reporting, grids and applicable plots; R demo/learner helpers remain R-only |
| `conditional-late-inference.Rmd` | [CLATE covariance and inference](conditional-late.md), [reporting](diagnostics-and-reporting.md) | Paired pointwise Wald inference/status; ratio confidence-set inversion remains R-only |
| `local-iv-estimation.Rmd` | [CLATE scope](conditional-late.md), this map | Experimental local IV-GMM and its simulation tools have no Stata command; CLATE is not a replacement for general GMM |
| `paper-affine-replication.Rmd` | This map and [inference scope](../README.md#inference-and-localization) | R's installed affine manuscript evidence/simulation/SIPP bundle is R-only; Stata examples are simulated tutorials, not paper replication |

## Reference families

Names in the first column match R's reference groups. The Stata package
reference is `help rdistnn`; each public command has its own installed help.

| R reference family | Applicable Stata interface/documentation | R-only or different surface |
|---|---|---|
| Package overview | `help rdistnn`, [article index](index.md), [README](../README.md) | No separate Stata pkgdown site |
| CATE and regression — estimators and weights | `cate_dnn`, `np_dnn`, `weights_dnn`; [core guide](cate-and-regression.md), [weights](weights-and-tuning.md) | `np_tdnn`, `weights_matrix`, batch weight objects and separate exported exact/asymptotic function names; Stata chooses the method by an option |
| S3 methods | Fit/replay, `r(table)`, `r(summary)`, `e(estimates)`, `e(estimate_vcov)`, diagnostic commands, stored grid matrices | No S3 objects, `candidate_fit()`/`candidate_fits()` or `predict()` method. Refit with a new x0/X0 to evaluate new targets; select stored grid rows to report existing candidates |
| Plotting | `cate_dnn_plot`, `clate_dnn_plot`, `np_dnn_plot`; [plot guide](visualising-fits.md) | No ggplot/autoplot object editing or independent plot-data generic. Stata single-view commands return `r(plot_data)` after building a graph |
| Kernel-order tuning | `lambda()`, `tuningmode(grid)`, `lambdagrid()`, `sgrid()`, `sdefault()` in CATE/CLATE/NP; [tuning](weights-and-tuning.md) | No `new_tuning`, `as_tuning`, registry or arbitrary evaluator. PLQR/standalone weights expose explicit s rather than integrated grids |
| Demo data & utilities | Deterministic [examples](index.md#start-with-a-task) generate their own data | No packaged `cate_demo_data`, `np_demo_data` or demo-subset helpers |
| General CATE simulation tools | Tutorial examples and [test map](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/README.md) | `sim_cate_*`, method registries and Monte Carlo summary utilities remain R-only; test success is not simulation evidence |
| Experimental — conditional LATE | `clate_dnn`, CLATE diagnostics/plots; [article](conditional-late.md) | R signal constructors and generic affine-ratio object methods map to inline supplied columns and Stata matrices. R `ratio_confidence_set`, learner/benchmark helpers and extra jackknife diagnostic exports have no Stata equivalents |
| Experimental — PLQR | `plqr_dnn`; [moment article](plqr-moments.md) | Single target and supplied theta moment covariance only. R's density-backed Jacobian/coefficient reporting at supplied theta, conditioning/calibration plots, reusable result methods and batch candidate interfaces are absent |
| Experimental — local IV-GMM (companion paper) | No Stata public command | `gmm_dnn`, `new_gmm_model`, DR-IV models, GMM result methods, IV demo/simulation helpers remain R-only |
| Developer API — scores (U-stat kernels) | Built-in CATE/NP/CLATE/PLQR score paths described in their guides | Score specifications, registration, generic evaluation and parameter-aware derivative paths remain R-only |
| Developer API — distances (rankers) | `distance(l2|mahalanobis)`; [weight guide](weights-and-tuning.md#exact-and-asymptotic-rules) | R L1/scaled-L1/scaled-L2 and custom distance specifications/registries are unavailable as Stata options |
| Developer API — rank → weight mappings | `method(exact|asymp)`, k/s controls | Custom rank-weight specifications, registries and two-scale mappings remain R-only |
| Advanced — TDNN and jackknife-paper infrastructure | [Weight/jackknife explanation](weights-and-tuning.md), internal algorithms tested against deletion references | No public TDNN estimator/weight/scales plot or generic `jackknife_fast_rank`, `jackknife_rank_full`, `jackknife_slow_rank` interface |

## Translating a shared workflow

| R convention | Stata convention |
|---|---|
| `data=` / matrix and vector inputs | Active dataset variables plus a named matrix for multiple targets |
| `x_cols`, `X0`, `tuning_mode`, `lambda_grid` | `xcols()`, `x0matrix()`, `tuningmode()`, `lambdagrid()` in CATE/CLATE/NP |
| `pseudo_outcome` | Positional signal with `pseudooutcome`; CLATE additionally requires `denominator()` |
| `level=.90` | `level(90)`; normally the configured Stata confidence level defaults to 95 |
| `coef()`, `summary()`, `confint()` | Full canonical e() matrices and r() tables from fit/replay; copy r() results promptly |
| `diagnostics(fit)` | Estimator-specific `_diagnostics` command; `.dta` saving and common localization return |
| `candidate_fit(fit, id)` | Inspect row id in stored grid matrices, or refit at its resolved s if a new active fit is required |
| `predict(fit, X0=...)` | Run the original estimator call on the required data with a new target matrix; replay does not predict |
| `plot_data_*()` | Build the applicable `_plot` view and copy `r(plot_data)`; a dashboard returns panel metadata |
| `se=FALSE` / `jackknife="none"` | `nose` in CATE/CLATE/NP; standalone weights have no SE, PLQR always computes moment covariance |

Old Stata aliases remain documented for compatibility. Removing underscores
is a useful naming guide for CATE/NP options, not a universal translator:
PLQR uses `x()`, `d()`, `g()`, `md()` and `thetaeval()`, while CLATE has its
own binary-IV and paired-signal grammar.

Stata's joint target covariance is useful for contrasts. It still does not
make default pointwise intervals simultaneous, and a joint CLATE covariance
matrix is different from R's marginal covariance-list presentation. Use
full canonical target matrices for CLATE because undefined coefficients are
omitted from the ordinary Stata coefficient post.

The [R-to-Stata test crosswalk](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/R_TEST_CROSSWALK.md) records numerical
and behavioral evidence. That test map and this documentation map answer
different questions: one shows which contracts are checked; the other tells
a reader where to learn and how to translate the available workflows.
