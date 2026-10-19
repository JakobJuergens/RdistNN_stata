# Save, report and reuse fitted results

Version 0.1.0.9002 provides three common commands for DNN/TDNN regression, CATE, CLATE and
PLQR. They use the R result methods as the reference for shared behavior.
PLQR retains its supplied-coefficient **moment-only** scope.

## Extract numeric reports

```stata
rdistnn_results, targets(3 1) level(90) saving("estimates.dta") replace
matrix estimates = r(fit)
matrix joint_covariance = r(vcov)
rdistnn_results, type(weights) targets(1) saving("weights.dta") replace
```

No graph is drawn. `r(fit)` and `r(se_fit)` have one row per component and
one column per selected target; NP/CATE/CLATE have one component. `r(table)`
and `r(summary)` follow the common reporting contract. `r(index)` identifies
the original target and component for every flattened covariance column.
The covariance is already on the estimator scale and retains cross-target
terms. Intervals remain pointwise. Undefined CLATE ratios retain their rows;
`nose` retains missing uncertainty.

`type(vcov)` exports a long covariance table with row/column target and
component IDs. `type(diagnostics)`, CLATE `type(components|status)`, and
CATE/CLATE `type(influence)` provide independent numeric extraction.
`r(data)` is the table written by `saving()`. The original active data and
every active estimation result are preserved, including on I/O errors.
See `help rdistnn_results` for exact column schemas. Estimator-specific
`_diagnostics` commands continue to provide their richer scientific reports
and all-candidate exports.

## Reuse a candidate fit

```stata
np_dnn Y, xcols(x1 x2) x0(0 0) tuningmode(grid) sgrid(5 10 15) sdefault(10)
rdistnn_candidate, candidate(3) store(order15)
estimates restore order15
np_dnn, level(90)
np_dnn_diagnostics
```

Candidate IDs are original rows of `e(tuning_candidates)`. Extraction reuses
cached ranks and distances and evaluates the candidate's resolved order with
the same frozen signals. Duplicate orders remain identical. The parent fit
stays active unless `post` is explicit; `store(name)` creates a reusable
named fit. Derived fits omit the parent's tuning grid and retain their own
model snapshot. Infeasible candidate IDs error without changing the parent.
Stata PLQR has no integrated candidate grid.

## Predict using saved observations and signals

```stata
estimates save fitted_model, replace
clear
estimates use fitted_model
matrix NEW = (-1,0 \ 0,0 \ 1,0)
rdistnn_predict, x0matrix(NEW) sefit level(90) store(new_targets)
matrix predicted = r(fit)
matrix covariance = r(vcov)
```

Targets have one row per profile and one column per fitted covariate, in
`e(xvars)` order. `x0()` supplies one profile. The default omits covariance;
`sefit` or `interval(confidence)` requests frozen-score jackknife covariance
and requires `s<n`. Predictions use the saved training sample, k, s, metric
and frozen signals. Mahalanobis geometry always uses that training sample.
New-target distances are computed once. A candidate at existing targets
uses the saved geometry. No nuisance predictions are refitted.

`rdistnn_predict` preserves active data and estimates by default. `post`
explicitly activates the derived fit; `store(name)` saves it independently.
A posted fit can be replayed, diagnosed, plotted and predicted again.
`e(sample)` is retained only when the original rows and their values/order
can still be verified; otherwise its sample is empty. Saved observation IDs
refer to the original fitted dataset, not row numbers of unrelated active data.

## Storage and experimental boundaries

Current fits store a versioned numeric snapshot in `e(model_*)`: fitted
inputs and covariates, frozen signals, targets, ranks, distances and original
observation IDs. `estimates save` writes this state into the `.ster` file.
No live Mata globals, original dataset, sidecar files or R installation are
needed after restoration. Storage grows with sample size and target count;
`.ster` files now contain the fitted numeric observations and signals.

Kernel and influence plots use the snapshot after data edits, reordering or
`clear`. Legacy fits without a snapshot retain their original-data/order
guards and require a refit for prediction or candidate extraction.

PLQR reporting and prediction concern localized score moments at the supplied
`theta_eval`. `r(quantity)` is `moment`, and `r(theta_eval)` reports the supplied
coefficients separately. Coefficient uncertainty, density-backed Jacobians,
root solving and R's PLQR conditioning/calibration plots remain outside this
upgrade. R-only estimators, custom metrics and ratio confidence sets also
remain outside the shared surface; see the [crosswalk](r-documentation-crosswalk.md).

## Two-scale regression

`np_tdnn` also supports `rdistnn_results` and `rdistnn_predict`, including
full covariance, saved-model restore and explicit posting. Both orders, the
localization dimension, normalization and fitted geometry remain fixed. Its
diagnostics describe signed weights, and no candidate grid is retained.
See [two-scale regression](two-scale-regression.md).
