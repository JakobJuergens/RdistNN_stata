# Experimental PLQR moment evaluation

`plqr_dnn` evaluates a localized partially linear quantile-regression score
at a **supplied coefficient vector**. It returns the moment and its full
rank-shift jackknife covariance. It does not find a coefficient, estimate a
threshold density, or report coefficient standard errors, confidence
intervals, p-values or significance tests. The PLQR interface is experimental
and can change. The initial release focuses on CATE and regression.

Use [examples/plqr-moments.do](../examples/plqr-moments.do) for complete seeded
data generation and extraction, or `help plqr_dnn` in an installed package.
The [documentation index](index.md) covers installation and other workflows.

## Specify the score and its units

Let `D` contain `q` partially linear regressors, `U` the full adjustment
state, `X` the reporting covariates, and `tau` a quantile index in `(0,1)`.
For observation `i`, component `j` is

```text
threshold_i = sum_j D_ij*theta_eval_j + g_i
psi_ij     = (tau - 1{Y_i <= threshold_i})*(D_ij - mD_ij)
moment_j(x) = sum_i w_i(x)*psi_ij
```

The threshold uses **less than or equal to**: observations exactly on it
receive indicator one. `g()` supplies the additive shape nuisance `g_0(U)`,
not a complete conditional quantile by itself. The complete threshold is
`D' theta_eval + g`. `md()` supplies `m_D(U)=E[D|U]` in exactly the same
variable order as `d()`. Neither nuisance is fitted by this command; their
predictions may depend on adjustment variables beyond `x()`.

`thetaeval()` supplies one coefficient per `d()` regressor in that order.
`Y` and `g` use outcome units; each coefficient has outcome units per unit of
its associated regressor. The indicator residual is dimensionless, so moment
component `j` has regressor-`j` units. A moment covariance entry `(j,h)` has
the product of those regressor units. It is not on a coefficient scale.

In particular, neither `e(b)` nor a small localized moment is a solved
quantile-regression coefficient. A finite sample need not give a zero moment
even when the supplied value is the generating coefficient. Assessing a
population moment condition or a coefficient's statistical validity requires
its own identification, nuisance, localization, bias and sampling conditions.

## One target, explicit variable order

The following fragment assumes the named data variables already exist:

```stata
plqr_dnn Y, x(x1 x2) x0(0 0) d(d1 d2) g(g_tau) ///
    md(md1 md2) thetaeval(1 -.5) tau(.5) k(1) s(12)
matrix moment = e(moment)
matrix moment_covariance = e(moment_vcov)
matrix supplied_theta = e(theta_eval)
```

Here `theta_eval=(1,-.5)` and `md1,md2` correspond respectively to `d1,d2`.
`x0()` supplies one coordinate per `x()` variable, in the same order.
Only one target is supported in each call: this command has neither
`x0matrix()` nor the `xcols()` alias used by the CATE/CLATE/NP front ends.
For several targets, make separate calls and store each moment/covariance
before the next call replaces active results. These calls do not provide a
joint covariance across targets.

The complete-case sample includes the outcome, all `x()` and `d()`
variables, `g()` and all `md()` variables within the requested `if`/`in`
range. Check `e(sample)` and `e(N)` after fitting. External predictions must
be correctly aligned with observations before the call.

Preferred input names are `g()` and `thetaeval()`. `q()` is a temporary
compatibility alias for `g()`; it does not change the shape nuisance into a
complete quantile. `beta()` is a temporary compatibility alias for
`thetaeval()`; it does not request coefficient estimation. Each preferred
name conflicts with its own alias.

## Localization and frozen moment covariance

Defaults are `tau(.5)`, `k(1)`, `distance(l2)`, `method(exact)`, and
`s=floor(.4*n^((p+1)/(p+2)))`, with complete-case sample size `n` and `p`
reporting covariates. `s()` is an integer subsample order, not a neighbor
count. `k()` counts observations averaged within each
subsample; all score components receive the same reporting weights.
`1<=k<=s<n` is required because PLQR always computes the jackknife.
There is no PLQR `nose` option or tuning-grid interface.

Exact weights can be positive through rank `n-s+k` (`n-s+1` under `k(1)`).
Larger `s` tightens localization, opposite to an ordinary bandwidth; exact
`k=s` gives uniform weights. `method(asymp)` uses normalized Beta-rank
weights. `distance(mahalanobis)` requires at least two reporting covariates
and uses the covariance geometry of the full estimation sample.

For each deletion the command removes one row, shifts the other rows' ranks,
and recomputes the weights at the same target, `s`, `k` and method. It freezes
the supplied coefficient, nuisance predictions, observation-level score and
distance metric, including full-sample Mahalanobis geometry. The covariance
retains off-diagonal dependence among the `q` score components:

```text
L_i = localized vector moment after deleting row i
Lbar = average_i L_i
moment_vcov = (n-1)/n * sum_i (L_i-Lbar)'*(L_i-Lbar)
```

`e(moment_vcov)` is already the covariance of the weighted local moment.
Do not divide it by `e(neff)` again. `e(neff)=1/sum_i w_i^2` describes weight
concentration, not an independent sample size that establishes a CLT.
The [weights and tuning article](weights-and-tuning.md) explains the shared
geometry and weight algorithms; its CATE/CLATE/NP grid options are not PLQR
options.

## Extract the returned objects

Stata uses a **row vector** for this single target. This differs from R's
moment-component-by-target matrix convention.

| Result | Dimensions/content |
| --- | --- |
| `e(moment)` | 1 by `q`, columns in `d()` order, row named `moment` |
| `e(moment_vcov)` | `q` by `q`, rows and columns named by `d()` |
| `e(theta_eval)` | 1 by `q`, supplied coefficients in `d()` order |
| `e(b)` | Compatibility copy of the moment, **not** the coefficient |
| `e(V)` | Compatibility copy of moment covariance |
| `e(beta)` | Temporary compatibility copy of the supplied evaluation vector |
| `e(N)`, `e(neff)` | Complete-case sample size and DNN effective sample size |
| `e(k)`, `e(s)`, `e(kx)`, `e(kd)`, `e(tau)` | Localization, reporting dimension, score dimension and quantile index |

The macros `e(theta_source)="supplied"`, `e(inference_target)="moment"`,
and `e(root_solver)="none"` make the scale explicit. `e(gvar)`, `e(mdvars)`,
`e(dvars)`, `e(xvars)`, `e(yvar)` and `e(x0)` record inputs;
`e(distance)` and `e(method)` record geometry/weights. `e(g_option)` and
`e(theta_option)` record preferred names or aliases, and `e(qvar)` retains
the resolved shape-variable compatibility name.

```stata
matrix m = e(moment)
matrix C = e(moment_vcov)
matrix theta_supplied = e(theta_eval)
scalar first_moment_sd = sqrt(C[1,1])
matrix list m
matrix list C
matrix list theta_supplied
display "First moment's jackknife SD: " first_moment_sd
```

This SD describes one score component. It is not the SE of the first
coefficient. Generic coefficient commands that interpret `e(b)`/`e(V)` as
coefficient estimates/covariance would use the wrong statistical object.
The custom display and replay print moments and covariance without
coefficient tests. `level()` is accepted but produces no interval or test.
There are no PLQR diagnostics or plot commands in Stata; use explicit
moment/covariance extraction with accurate labels.

## What remains R-only

The R companion also accepts positive supplied threshold-density values and
constructs a density-backed coefficient sandwich with analytic or full
jackknife moment meat. That operation needs a coefficient Jacobian in
addition to moment covariance. Stata `plqr_dnn` has no `density()` input,
Jacobian, sandwich route or coefficient-inference fields. Taking the square
root of a moment-variance diagonal cannot substitute for that missing bread.

R's multiple-target support, tuning grids, retained score/weight extras and
PLQR-specific reporting/prediction/plot methods likewise do not establish
those interfaces in Stata. Neither front end solves a local PLQR root.
For implemented Stata scope, use this article and installed help; for
comparison, see the source-pinned R [PLQR reference](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/man/plqr_dnn.Rd),
[conditional Z-estimation article](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/vignettes/conditional-z-estimation.Rmd),
and [notation article](https://github.com/JakobJuergens/RdistNN/blob/6fe48273aaa8853c5736841c5a228f8c1592ded9/vignettes/notation-and-conventions.Rmd).
The [documentation crosswalk](r-documentation-crosswalk.md) records the full
package boundary.

## Saved moment reporting and prediction

`rdistnn_results` returns moment tables with supplied `r(theta_eval)` separately.
`rdistnn_predict, x0matrix(NEW) sefit` evaluates the saved frozen scores at a
batch of new targets and returns joint **moment** covariance. `post` or
`store(name)` creates a reusable derived moment fit. These commands work after
`estimates save`/`estimates use` without the original active dataset. They do
not solve for coefficients or supply coefficient uncertainty. See [reusing
fits](reusing-fits.md) and installed `help rdistnn_predict`.
