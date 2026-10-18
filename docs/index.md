# RdistNN for Stata documentation

These guides describe package **0.1.0.9000** on the local `development` branch, including experimental
conditional LATE and the shared reporting and postestimation interfaces.
See [installation](../README.md#installation). In Stata, `help rdistnn` opens
the package overview and links to every installed command reference.

## Start with a task

| Task | Guide | Runnable example |
|---|---|---|
| Install, fit and extract a first result | [Getting started](getting-started.md) | [Original quick start](../examples/getting-started.do) |
| Choose variables and understand shapes/names | [Notation and conventions](notation-and-conventions.md) | [CATE and regression](../examples/cate-and-regression.do) |
| Supply CATE scores or nuisances; estimate conditional means | [CATE and regression](cate-and-regression.md) | [CATE and regression](../examples/cate-and-regression.do) |
| Inspect/reuse weights; choose kernel orders and candidate grids | [Weights and tuning](weights-and-tuning.md) | [Weights and tuning](../examples/weights-and-tuning.do) |
| Save fits, extract tables and evaluate new targets | [Reusing fits](reusing-fits.md) | [Saved results](../examples/reusing-fits.do) |
| Read tables, assess localization, save diagnostic data | [Diagnostics and reporting](diagnostics-and-reporting.md) | [Diagnostics and plots](../examples/diagnostics-and-plots.do) |
| Choose panels, compare targets/candidates, export graphs | [Visualising fits](visualising-fits.md) | [Diagnostics and plots](../examples/diagnostics-and-plots.do) |
| Estimate experimental conditional LATE with supplied IV predictions | [Conditional LATE](conditional-late.md) | [Oracle binary-IV example](../examples/clate.do) |
| Evaluate experimental PLQR moments at a supplied coefficient | [PLQR moments](plqr-moments.md) | [PLQR moments](../examples/plqr-moments.do) |
| Translate R workflows and identify unavailable features | [R documentation crosswalk](r-documentation-crosswalk.md) | [Numerical test crosswalk](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/R_TEST_CROSSWALK.md) |

Each example is a standalone do-file that clears the active dataset and sets
its seed. From a checkout, change to the repository root and run
`do "dev-setup.do"` once, then run an example. With the package
installed, download an example and run it directly. The diagnostics example
creates its optional saved data/graphs under `build/documentation/` relative
to the current directory. Read its opening comments before running it.

## Command references

| Commands | Installed help | Role |
|---|---|---|
| `cate_dnn` | `help cate_dnn` | AIPW, Robinson-form or supplied-score CATE |
| `np_dnn` | `help np_dnn` | Nonparametric conditional means |
| `weights_dnn` | `help weights_dnn` | One target's exact/asymptotic DNN weights |
| `clate_dnn` | `help clate_dnn` | Experimental augmented binary-IV or paired-signal ratios |
| `rdistnn_results`, `rdistnn_candidate`, `rdistnn_predict` | `help COMMAND` | Graph-free extraction and reusable saved fits |
| `plqr_dnn` | `help plqr_dnn` | Experimental supplied-coefficient score moments |
| `cate_dnn_diagnostics`, `clate_dnn_diagnostics`, `np_dnn_diagnostics` | `help COMMAND` | Stored target/candidate reports and `.dta` export |
| `cate_dnn_plot`, `clate_dnn_plot`, `np_dnn_plot` | `help COMMAND` | Stored-result and observation-level graph views |

The [command guide](command-guide.md) is a compact workflow and option map.
The installed help is the full syntax/result reference; the articles explain
how to make choices and interpret the outputs. RdistNN is the numerical
reference for shared implementations, rather than a promise that every R
object method or experimental workflow is available in Stata.

## Inference and implementation details

Start with the estimator-specific inference sections in
[CATE and regression](cate-and-regression.md), [CLATE](conditional-late.md),
and [PLQR](plqr-moments.md). The [weight guide](weights-and-tuning.md#jackknife-and-the-localized-target)
connects rank weights, deletion covariance and localization bias. Valid
pointwise inference requires the relevant design, nuisance, stability,
variance and bias conditions; a successful command is not a check of those
conditions. Weighted outcome quantiles are distribution summaries, and
neither ESS nor a candidate grid supplies an inference guarantee.

For contributing, use [CONTRIBUTING.md](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md), the
[test-suite map](https://github.com/JakobJuergens/RdistNN_stata/blob/master/tests/README.md), and the
[native plugin guide](../src/plugin/README.md). [RELEASING.md](https://github.com/JakobJuergens/RdistNN_stata/blob/master/RELEASING.md)
describes release validation. Dated [CLATE verification](clate-verification.md)
and [estimator-coherence verification](estimator-coherence-verification.md)
record the earlier implementation checks; they are separate from current
user instructions. The [documentation verification receipt](documentation-verification.md)
records the earlier article/help/example update. [Release notes](../RELEASE_NOTES.md)
record package changes and experimental features.
