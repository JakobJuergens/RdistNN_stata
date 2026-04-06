# Updated audit note on the Stata/Mata layer in `RdistNN_stata`

I reviewed the Stata/Mata sources in the Stata archive and then compared them against the newly provided **R package source**, which appears to be the reference implementation the Stata package is based on. This is still a **static code audit** rather than an executed test run, but the R reference substantially sharpens several of the earlier findings.

## Updated bottom line

After comparing against the R package, my view is stronger than before:

1. **The Stata port is not just incomplete in a few places; in several spots it has drifted away from the reference implementation.**
2. The most important Stata-side problems are now even clearer:
   - **`distance(mahalanobis)` is advertised but not implemented** in the Stata/Mata layer.
   - **The jackknife standard errors in Stata are not the same algorithm as in the R reference.**
   - **`cate_dnn` in Stata does not match the reference package’s CATE score interface and seems to collapse distinct score constructions into a single Robinson-style path.**
3. Some of the earlier bugs remain fully valid:
   - reused global helper name `Display`
   - `weights_dnn` fallback bug with `generate()`
   - `weights_dnn` incorrectly requiring at least two covariates
   - missing `markout` in `weights_dnn`

So the main concern is no longer only “there are coding bugs,” but also **“the Stata implementation is not faithfully reproducing the reference package’s estimator stack.”**

---

## What the R reference clarified

The R package is helpful because it spells out the intended architecture much more explicitly.

In the R code:

- distances are handled through a proper registry/API (`R/api-distances.R`),
- Mahalanobis distance is genuinely implemented (`R/dist_mahal.R`),
- DNN weights are separated into exact and asymptotic rank-weight rules (`R/wrank-dnn_exact.R`, `R/wrank-dnn_approx.R`, `R/wrank-dnn-wrapper.R`),
- jackknife SEs are implemented through actual leave-one-out logic, with both a literal `loo` mode and a faster but still exact rank-update path (`R/user-dnn-core.R` and related files),
- and these behaviors are backed by tests, especially in
  `tests/testthat/test-distance-mahal.R` and
  `tests/testthat/test-user-np-dnn.R`.

That makes it much easier to distinguish true Stata bugs from deliberate design choices.

---

## Most important updated findings

### 1. `distance(mahalanobis)` is definitely a real implementation gap, not just a missing optimization

My earlier note said that `distance(mahalanobis)` was accepted by the Stata interface but not actually used. The R reference makes that conclusion much stronger.

In the **R package**, Mahalanobis distance is fully implemented:

- `R/dist_mahal.R` defines `distance_mahalanobis()` with a `setup()` step that computes `cov(X)^{-1}` once, then reuses it in `eval()` / `eval_batch()`.
- It explicitly handles singular or near-singular covariance matrices through a pseudoinverse fallback.
- `tests/testthat/test-distance-mahal.R` contains a fairly serious test suite for correctness, setup behavior, singular covariance, alias registration, and integration with `cate_dnn()`.

By contrast, in the **Stata/Mata** code:

- the `.ado` interfaces accept `distance(mahalanobis)`,
- the help files describe it as available,
- but the actual computation always goes through L2:
  - plugin calls are always `weights_l2 ...`
  - `_rdistnn_mata.ado` only defines `rdistnn_l2dist()`
  - `rdistnn_weights()` always calls the L2 distance path.

So this is not just “Mahalanobis support hasn’t been optimized yet.” It is a **clear behavioral mismatch with the reference package**. A user calling `distance(mahalanobis)` in Stata is currently not getting the estimator they asked for.

### 2. The Stata jackknife is not faithful to the reference implementation

This is the biggest substantive finding from the comparison.

In the **R package**, the jackknife logic is explicit and careful:

- `R/user-dnn-core.R` implements two modes:
  - `jackknife = "loo"`: literal recomputation on the leave-one-out sample,
  - `jackknife = "fast"`: reuse full-sample ranks, update ranks after deletion, then recompute the DNN weights on the reduced sample.
- The fast path is still a genuine delete-1 calculation, just done efficiently.
- The tests in `tests/testthat/test-user-np-dnn.R` check that the fast and loo paths agree.
- The code also explicitly requires **`s < n`** for delete-1 jackknife validity with fixed order `s`.

In the **Stata** code, `np_dnn.ado` and `cate_dnn.ado` compute leave-one-out values as

```text
(theta_hat - w_i * signal_i) / (1 - w_i)
```

using the full-sample weights `w_i`.

That is the correct deletion formula for a **fixed weighted average with renormalization**, but it is **not** the same as recomputing a DNN rank-weight estimator after deletion. DNN weights depend on:

- the ranks,
- the sample size,
- and the rank shifts induced by dropping an observation.

The R reference confirms that the intended jackknife should recompute those weights on the reduced sample, either directly or through a rank-update shortcut. The Stata code does not do that.

So relative to the reference package, the Stata SEs are not just “fast”; they are using a **different jackknife algorithm**.

That matters because the help files say “fast leave-one-out jackknife,” which strongly suggests fidelity to delete-1 recomputation. Right now that description is misleading.

### 3. `cate_dnn` in Stata appears to collapse distinct CATE score constructions that the R package treats separately

This was not fully obvious before seeing the R package.

In the **R reference**, `cate_dnn()` clearly supports two distinct pseudo-outcome constructions:

1. **Full AIPW score** using `m0`, `m1`, and `e`:

```text
psi = (D (Y - m1) - (1-D)(Y - m0)) / (e(1-e)) + (m1 - m0)
```

2. **Robinson-style score** using `m = E[Y|X]` and `e`:

```text
psi = (Y - m)(D - e) / (e(1-e))
```

The R interface is explicit about this distinction in `R/user-dnn-cate.R`.

In the **Stata** implementation, `cate_dnn.ado` only constructs

```text
rho_i = (D_i - e_i)(Y_i - m_i) / (e_i (1-e_i))
```

where `m_i` is built as

```text
m = m0 + e (m1 - m0)
```

if both `m0()` and `m1()` are supplied.

So the Stata code is effectively forcing everything through a Robinson-style score based on `m = E[Y|X]`, even when the user supplies `m0` and `m1` separately.

That is a **real estimator-level deviation** from the R reference, which treats the full AIPW score as its own supported path.

Even more concerning, the Stata code allows cases where only one of `m0()` or `m1()` is supplied and then substitutes:

- if only `m0()` is present, it sets `m = m0`
- if only `m1()` is present, it sets `m = m1`

This fallback is not mirrored in the R reference and is statistically hard to justify as a default behavior. If the intended score is Robinson-style, the natural nuisance is `m(X) = E[Y|X]`, not one arm-specific mean substituted ad hoc.

So for `cate_dnn`, I think the Stata implementation is currently mixing together:

- the public interface of the richer R function,
- the Robinson score construction,
- and a fallback policy for missing nuisance pieces that the reference code does not endorse.

### 4. The Stata help/examples for `cate_dnn` reinforce the score confusion

The help file says the estimator uses Robinson pseudo-outcomes, which is fine in itself. But the example then calls

```stata
propensity(e_oracle) m0(m0_true) m1(m0_true)
```

in a design where treatment effects are nonzero.

If the command is using a Robinson score with `m = E[Y|X]`, then passing `m1 = m0_true` is generally wrong, because it makes

```text
m = m0 + e(m1 - m0) = m0
```

instead of the actual unconditional regression function.

If the command were using the full AIPW score, then passing `m1 = m0_true` would also generally be wrong, because `m1` should be `E[Y|X,D=1]`.

So the example is not just a little loose; it suggests the Stata-side treatment of nuisance functions has not been cleanly aligned with either score formulation.

### 5. `weights_dnn` incorrectly blocks one-dimensional DNN even though the reference package clearly supports it

In Stata, `weights_dnn.ado` declares

```stata
syntax varlist(min=2)
```

which forbids one-dimensional covariates.

The R package clearly does **not** impose this restriction for ordinary DNN weighting. One-dimensional input is treated as perfectly valid throughout the API.

The only place where dimension matters specially is Mahalanobis distance, where `R/dist_mahal.R` explicitly requires at least two covariates.

So the Stata restriction is almost certainly accidental and too strong. It should not be there for L2-based DNN weights.

### 6. Missing-value handling in `weights_dnn` is weaker than in both the other Stata commands and the R reference

This issue remains valid and the R package makes the gap more visible.

- `np_dnn.ado` and `cate_dnn.ado` call `markout`.
- `weights_dnn.ado` does not.
- The R package consistently validates shapes and missingness before distance/ranking operations.

So `weights_dnn` is the odd command out in the Stata port. Missing covariates can leak into ranking/weight calculations there in a way that is inconsistent with the rest of the port and with the reference package’s more defensive design.

### 7. The Stata port is missing the reference package’s explicit `s < n` jackknife guard

In the R code, delete-1 jackknife logic explicitly refuses `s >= n`:

```text
Delete-1 jackknife requires s < n
```

That makes sense because the leave-one-out sample has size `n-1`, so keeping the same order `s` requires `s <= n-1`.

The Stata code does not impose this guard in `np_dnn` or `cate_dnn`. Given that its current SE formula is already not a literal leave-one-out recomputation, this omission is another sign that the Stata jackknife path is operating under a different, looser logic than the reference implementation.

---

## Earlier findings that still stand unchanged

These issues are still present and still worth fixing even apart from the R comparison.

### 8. Reused helper program name `Display`

Both `np_dnn.ado` and `cate_dnn.ado` define a global Stata program named `Display`. That is fragile because Stata program names live in a shared namespace. These should be renamed to command-specific helpers such as:

- `np_dnn_display`
- `cate_dnn_display`

### 9. `weights_dnn` fallback bug with `generate()`

If the plugin path is attempted, pre-creates the target variable, and then fails, the Mata fallback tries to `gen` the same variable again. That can cause a secondary failure exactly when the fallback is supposed to rescue the call.

### 10. `np_dnn` still does not properly validate or use `distance()`

Compared with the R package, this now looks even more clearly like an unfinished port. `np_dnn.ado` exposes `distance()` but does not really implement a distance dispatch layer at all.

### 11. `np_dnn` does not explicitly validate `k >= 1`

The other Stata entry points do. `np_dnn` should too.

### 12. The package metadata still appears incomplete for deployment

The Stata package metadata does not appear to list the plugin binary for installation, and the README refers to setup material that is not in the archive. These are secondary issues, but they do affect whether the Stata package is reproducibly installable.

---

## How I would reinterpret the Stata port after seeing the R package

Before seeing the R code, one could read the Stata package as a rough early version with some bugs.

After seeing the R code, I think the more accurate diagnosis is:

- the **weight formulas** in Mata look broadly aligned with the intended DNN exact/asymptotic rules,
- but the **distance layer**, **jackknife layer**, and **CATE score layer** are only partially ported,
- and some user-facing interfaces were copied over before the underlying semantics were fully implemented.

That is why the package feels inconsistent: some parts are reasonably faithful, while others still look like placeholders or compressed approximations of richer R functionality.

---

## What I would fix first

If the goal is to make the Stata implementation genuinely faithful to the reference package, I would prioritize fixes in this order:

1. **Disable `distance(mahalanobis)` immediately unless and until it is truly implemented.**
2. **Decide what `cate_dnn` is supposed to estimate.**
   - If it is Robinson-only, simplify the interface and require a single nuisance `m(X)`.
   - If it is supposed to mirror the R package, add distinct support for both Robinson and full AIPW pseudo-outcomes.
3. **Replace the current SE computation with a genuine delete-1 jackknife path.**
   - literal recomputation is fine as a first correct implementation,
   - then add the rank-update shortcut later if desired.
4. **Rename the `Display` helpers.**
5. **Fix the `weights_dnn generate()` fallback bug.**
6. **Allow one-dimensional covariates in `weights_dnn`.**
7. **Add `markout` and stronger front-end validation (`k >= 1`, `s <= n`, and for jackknife specifically `s < n`).**

---

## Short verdict

Yes, I still think there are real problems in the Stata/Mata layer, and the R reference makes that conclusion stronger.

The most important updated takeaways are:

- **Mahalanobis support is genuinely missing in Stata while being fully implemented in R.**
- **The Stata jackknife is not the same as the reference package’s delete-1 jackknife.**
- **`cate_dnn` in Stata does not cleanly match the reference package’s CATE score construction and likely implements a different estimator than the interface suggests.**

So the central issue is not only code robustness. It is also **implementation fidelity**: the Stata port currently does not look like a fully faithful translation of the R package it is based on.
