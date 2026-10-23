# RdistNN for Stata documentation

Start with [Getting started](getting-started.md) for a first fit. These guides
cover package **0.1.3**. In Stata, `help rdistnn` opens the command
references; [release notes](../RELEASE_NOTES.md) record changes.

## Start with a task

| Task | Guide | Runnable example |
|---|---|---|
| Install, fit and extract a first result | [Getting started](getting-started.md) | [First fit](../examples/getting-started.do) |
| Choose variables and understand shapes/names | [Notation and conventions](notation-and-conventions.md) | [CATE and regression](../examples/cate-and-regression.do) |
| Supply CATE scores or nuisances; estimate conditional means | [CATE and regression](cate-and-regression.md) | [CATE and regression](../examples/cate-and-regression.do) |
| Inspect/reuse weights; choose kernel orders and candidate grids | [Weights and tuning](weights-and-tuning.md) | [Weights and tuning](../examples/weights-and-tuning.do) |
| Save fits, extract tables and evaluate new targets | [Reusing fits](reusing-fits.md) | [Saved results](../examples/reusing-fits.do) |
| Read tables, assess localization, save diagnostic data | [Diagnostics and reporting](diagnostics-and-reporting.md) | [Diagnostics and plots](../examples/diagnostics-and-plots.do) |
| Choose panels, compare targets/candidates, export graphs | [Visualising fits](visualising-fits.md) | [Diagnostics and plots](../examples/diagnostics-and-plots.do) |
| Estimate conditional LATE with supplied IV predictions | [Conditional LATE](conditional-late.md) | [Oracle binary-IV example](../examples/clate.do) |
| Evaluate experimental PLQR moments at a supplied coefficient | [PLQR moments](plqr-moments.md) | [PLQR moments](../examples/plqr-moments.do) |
| Fit two-scale regression with signed weights | [Two-scale regression](two-scale-regression.md) | [Two-scale example](../examples/two-scale-regression.do) |
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
| `cate_dnn` | `help cate_dnn` | AIPW or supplied-score CATE |
| `np_dnn` | `help np_dnn` | Nonparametric conditional means |
| `weights_dnn` | `help weights_dnn` | Single/batch exact/asymptotic DNN weights and geometry |
| `np_tdnn`, `weights_tdnn` | `help COMMAND` | Signed two-scale regression and weights |
| `np_tdnn_diagnostics`, `np_tdnn_plot` | `help COMMAND` | Signed-weight diagnostics and two-scale plots |
| `clate_dnn` | `help clate_dnn` | Augmented binary-IV or paired-signal ratios |
| `rdistnn_results`, `rdistnn_candidate`, `rdistnn_predict` | `help COMMAND` | Graph-free extraction and reusable saved fits |
| `plqr_dnn` | `help plqr_dnn` | Experimental supplied-coefficient score moments |
| `cate_dnn_diagnostics`, `clate_dnn_diagnostics`, `np_dnn_diagnostics` | `help COMMAND` | Stored target/candidate reports and `.dta` export |
| `cate_dnn_plot`, `clate_dnn_plot`, `np_dnn_plot` | `help COMMAND` | Stored-result and observation-level graph views |

The [command guide](command-guide.md) is a compact workflow and option map.
The installed help is the full syntax/result reference; the articles explain
how to make choices and interpret the outputs. Shared implementations follow
the RdistNN contract, with bidirectional engine comparisons and independent
references. The [crosswalk](r-documentation-crosswalk.md) identifies R-only
methods, and the [numerical reference](numerical-stability.md#asymptotic-tail-limits)
records the unresolved asymptotic-tail support limitation.

## Theory and development

The [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
gives the statistical assumptions and results. The estimator guides explain
how to supply inputs, inspect local support and interpret uncertainty outputs.

For development and verification, see
[CONTRIBUTING.md](https://github.com/JakobJuergens/RdistNN_stata/blob/master/CONTRIBUTING.md).
The [plugin guide](../src/plugin/README.md) covers optional acceleration;
[numerical implementation](numerical-stability.md) is a specialist reference.
