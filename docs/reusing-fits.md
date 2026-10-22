# Save, report and reuse fitted results

Version 0.1.2 provides three common commands for DNN/TDNN regression, CATE, CLATE and
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

**Saved-fit confidentiality:** `estimates save` writes observation-level numeric
inputs and frozen signals into the `.ster` file. Apply the underlying data's
sharing restrictions to saved fits; clearing the active dataset does not remove
those observations from the saved model. Do not distribute confidential fits.

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

Fits use exact distance followed by lexicographic original covariates, then
observation order for identical covariate rows as bookkeeping. Final weights
are averaged within the identical-covariate blocks, including every deletion's
remaining block. Cached ranks and retained raw covariates jointly determine
these weights. The command warns that this modification is outside the paper's
theory. Saved-model schema version 2 records the policy; refit older saved objects
before reusing them with the current package. Current snapshots are also
validated for consistent dimensions, controls, observation/target identities,
frozen geometry and cached ordering before reconstruction. An invalid snapshot
requests refitting before prediction or observation-level extraction and
preserves the active data and fitted results. Replaying already stored numeric
estimates does not reconstruct the old model under the new weight policy.

`rdistnn_predict` preserves active data and estimates by default. `post`
explicitly activates the derived fit; `store(name)` saves it independently.
A posted fit can be replayed, diagnosed, plotted and predicted again.
`e(sample)` is retained only when the original rows and their values/order
can still be verified; otherwise its sample is empty. Saved observation IDs
refer to the original fitted dataset, not row numbers of unrelated active data.

For NP, TDNN, CATE and CLATE, derived fits preserve `e(depvar)` from the
current fitted state. Default reuse leaves the parent's complete `e()` state
unchanged. `e(vce)` and `e(vcetype)` describe the requested
covariance calculation, with `none` and an empty label for point-only results.
For CLATE, `e(status)` still determines whether ratio inference is available.

## Migrating saved CATE fits after Robinson removal

Fresh nuisance-based CATE uses AIPW and requires a propensity and both
arm-specific outcome means. A pooled outcome mean is not a substitute.
Alternatively, `pseudooutcome` localizes a supplied signal unchanged; its
construction, causal interpretation and inferential justification remain
the caller's responsibility. The former `mcol()` and `mhat()` options remain
rejected. The paper's residualized-slope estimating equation is distinct from
the removed normalized Robinson pseudo-outcome; this migration adds no
residualized-slope estimator.

Refit objects created by older development versions. Postestimation uses the
current numeric snapshot and does not reconstruct retired score methods from
live inputs. Current supplied signals are reused unchanged and keep their
recorded provenance.

The R companion additionally reports descriptive signal provenance and the
original fit's creation version. In either interface, provenance does not certify same-sample
nuisance accuracy or stability, identification, or localization conditions.
Localization bias and the pointwise scope of inference remain relevant; see
the [CATE uncertainty guide](cate-and-regression.md#what-the-uncertainty-measures).

## Storage and experimental boundaries

Current fits store a versioned numeric snapshot in `e(model_*)`: fitted
inputs and covariates, frozen signals, targets, ranks, distances and original
observation IDs. `estimates save` writes this state into the `.ster` file.
No live Mata globals, original dataset, sidecar files or R installation are
needed after restoration. Storage grows with sample size and target count;
`.ster` files now contain the fitted numeric observations and signals.

Kernel and influence plots use the current snapshot after data edits,
reordering or `clear`. Fits without the current required snapshot need refitting.

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
posted coefficient names remain `target1`, `target2`, and so on, including
`target1` for a single target, so existing coefficient contrasts remain valid.
Its diagnostics describe signed weights, and no candidate grid is retained.
See [two-scale regression](two-scale-regression.md).
