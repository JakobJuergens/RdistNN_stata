# Two-scale nonparametric regression

`np_tdnn` implements the R companion's `np_tdnn()` for conditional mean
regression. It combines two DNN rank-weight tables at orders `s1<s2`, using
the localization dimension `dim`. The combination is designed to cancel a
leading localization-bias term under the corresponding smoothness conditions.
The weights can be negative. This is a separate regression workflow; CATE
and CLATE continue to use single-scale DNN.

The [complete example](../examples/two-scale-regression.do) generates seeded
data and demonstrates fitting, diagnostics, plots and saved-fit prediction.
Installed `help np_tdnn` and `help weights_tdnn` document every option.

```stata
matrix X0 = (-1,0 \ 0,0 \ 1,0)
np_tdnn Y, xcols(x1 x2) x0matrix(X0) s1(4) s2(12)
matrix estimates = e(estimates)
matrix joint_covariance = e(estimate_vcov)
np_tdnn_diagnostics
np_tdnn_plot, type(estimates) replace
np_tdnn_plot, type(weights) replace
```

`s1()` and `s2()` default to 1 and 2. `dim()` defaults to the number of
localization covariates. Exact and asymptotic rules use the same two-scale
combination as R. `nonormalize` preserves R's unnormalized asymptotic tables;
it has no effect on the exact tables. Targets use either `x0()` or
`x0matrix()`. Distance and scaling controls follow the
[shared geometry contract](weights-and-tuning.md#distance-geometry).

The full-sample weights use the fitted training geometry. Each delete-one
calculation removes an observation, shifts retained ranks, and rebuilds both
rank tables at the same `s1`, `s2` and `dim`. Geometry remains fixed.
`e(estimate_vcov)` retains covariance across targets and is already on the
estimator scale. Covariance requires `s2<n`; `nose` permits `s2=n` and leaves
uncertainty unavailable. Pointwise normal intervals require the relevant
sampling, variance, smoothness and bias conditions. They are not simultaneous
bands, and the two-scale construction does not establish those conditions.

## Interpret signed weights

The signed diagnostic table reports positive, negative and absolute weight
mass, the signed sum, squared-weight concentration, maximum absolute weight
and nonzero support. An inverse squared-weight concentration is an algebraic
concentration measure; negative weights prevent an ordinary probability-mass
interpretation. The distance summary uses absolute weights and is labelled
accordingly. Probability-weighted outcome quantiles and cumulative probability
shares are not supplied for these signed weights.

`weights_tdnn` exports the observation-by-target weight, rank and distance
matrices together with original observation identifiers. Its `generate()`
option writes a weight variable for a single target. For example:

```stata
weights_tdnn x1 x2, x0matrix(X0) s1(4) s2(12)
matrix W = r(weights)
matrix ranks = r(ranks)
matrix distances = r(distances)
```

## Reuse the fitted model

The shared numeric commands also accept two-scale fits:

```stata
rdistnn_results, type(vcov)
rdistnn_results, type(weights)
estimates store regression
matrix NEW = (-.5,0 \ .5,0)
rdistnn_predict, x0matrix(NEW) sefit store(new_profiles)
```

Saved numeric state includes both orders, the localization dimension, the
normalization choice, outcomes and fitted geometry. Prediction and reporting
work after `estimates save`/`estimates use` without the original dataset.
Prediction preserves active estimates unless `post` is requested. This
workflow has no integrated tuning grid; `rdistnn_candidate` reports that
no grid is available. See [reusing fits](reusing-fits.md).
