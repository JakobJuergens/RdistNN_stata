# Stata implementation: what needs to be fixed next

Status update, 2026-04-11: the first implementation pass has addressed the
highest-priority estimator issues described below: `cate_dnn` now uses
arm-specific AIPW denominators, requires `propensity()`, `m0()`, and `m1()`,
rejects non-binary treatment and out-of-range propensity values, and the
package manifest/help files have been updated.  The Mahalanobis jackknife is
now documented as a full-sample fixed-metric rank-shift jackknife rather than
as a changing-covariance leave-one-out recomputation.  The remaining strategic
question is whether to later implement a slower true leave-one-out
Mahalanobis jackknife.  Status update, 2026-05-30: the Stata-version
hygiene pass standardizes the active ado/test entry files and user-facing
documentation to Stata 17 or later.  Status update, 2026-05-30: the next
non-GMM Z-estimation layer is `plqr_dnn`, a fixed-beta localized PLQR
moment evaluator backed by a vector rank-shift jackknife covariance helper.
This keeps `gmm_dnn`, TDNN, automatic nuisance fitting, diagnostics, and true
changing-metric Mahalanobis leave-one-out jackknife as later strategic layers.
Status update, 2026-06-01: the active next layer is install and release
hardening: repeatable local Stata runtime wrappers, fresh `net install`
smoke checks outside the development adopath, and source-archive validation.
This remains a non-estimator layer; `gmm_dnn`, TDNN, automatic nuisance
fitting, diagnostics, and true changing-metric Mahalanobis leave-one-out
jackknife remain later strategic layers.  During this pass, the installed
plugin was made discoverable by name, but computation continues to default to
the Mata path; plugin fast-path parity for all `x0()` values is a separate
future hardening item.

The rest of this note is retained as historical context for the review that
motivated the 2026-04-11 implementation pass.

## 1. Fix the `cate_dnn` pseudo-outcome first

### Why this is the top priority

The current `cate_dnn` implementation still appears to construct the wrong signal. That is not a cosmetic issue. It changes the estimator itself.

At present, the full-nuisance branch computes a pseudo-outcome of the form

```stata
rho = D*(Y-m1)/[e(1-e)] - (1-D)*(Y-m0)/[e(1-e)] + (m1-m0)
```

But the usual AIPW / doubly robust signal for a CATE-style target is

```stata
rho = D*(Y-m1)/e - (1-D)*(Y-m0)/(1-e) + (m1-m0)
```

Equivalently,

```stata
rho = (D-e)*(Y-m_D)/[e(1-e)] + (m1-m0),
```

where `m_D = D*m1 + (1-D)*m0`.

The current code divides both residual terms by the common product `e(1-e)`, which over-scales both arms relative to the usual score.

### What should be changed

In `src/cate_dnn.ado`, replace the current common-denominator construction with separate denominators for the treated and control residual terms.

The safest implementation is to clip the propensity once, then use `e` and `1-e` separately.

### Recommended implementation pattern

Use something like:

```stata
tempvar eclip rho_var

quietly gen double `eclip' = min(max(`propensity', 1e-8), 1 - 1e-8) if `touse'

quietly gen double `rho_var' = ///
    `treat' * (`outcome' - `m1') / `eclip' ///
    - (1 - `treat') * (`outcome' - `m0') / (1 - `eclip') ///
    + (`m1' - `m0') ///
    if `touse'
```

This is much cleaner than building a single `den_var = e(1-e)` and then using it everywhere.

### What to do about the current fallback branch

The current one-nuisance fallback also looks conceptually wrong:

```stata
rho = (D-e)*(Y-m_use)/[e(1-e)]
```

with `m_use` equal to either `m0` or `m1`.

That expression only makes theoretical sense if `m_use` is a single pooled regression function such as `m(X) = E[Y|X]`, not a treatment-arm-specific regression like `m0(X)` or `m1(X)`.

So one of the following needs to happen:

#### Preferred immediate fix
Require **both** `m0()` and `m1()` and remove the current fallback altogether.

That is the safest choice if the goal is to keep the Stata package aligned with the intended estimator.

#### Alternative, if single-nuisance support is genuinely desired
Redesign the API and the theory around a pooled nuisance, for example with an explicit option such as

```stata
m(varname)
```

rather than overloading `m0()` or `m1()` into a role they do not actually play.

### Practical recommendation

For the next revision, I would do the conservative thing:

- make `propensity()` mandatory,
- make **both** `m0()` and `m1()` mandatory,
- delete the current one-nuisance branch,
- reintroduce a fallback only after the target parameter and score are written down clearly.

---

## 2. Exclude missing outcomes explicitly in both `cate_dnn` and `np_dnn`

### Why this matters

Both estimators currently mark out the covariates and nuisance variables, but not the outcome itself.

That means an observation with missing `Y` can survive into `touse`, then get passed into Mata through `st_data()`, and contaminate the weighted sum or jackknife.

This is a real implementation bug, not just a stylistic issue.

### What should be changed

#### In `src/np_dnn.ado`

Current pattern:

```stata
marksample touse
markout `touse' `x'
```

This should become:

```stata
marksample touse
markout `touse' `varlist' `x'
```

or, equivalently,

```stata
local outcome `varlist'
marksample touse
markout `touse' `outcome' `x'
```

#### In `src/cate_dnn.ado`

Current pattern:

```stata
marksample touse
markout `touse' `treat' `x' `propensity' `m0' `m1'
```

This should include the outcome:

```stata
marksample touse
markout `touse' `varlist' `treat' `x' `propensity' `m0' `m1'
```

### Why this should be fixed immediately

It is a low-effort change with high payoff. It protects the entire estimation path and the jackknife from silently ingesting missing values.

---

## 3. Tighten input validation in `cate_dnn`

### Why this matters

Right now `cate_dnn` says that `treat()` is binary and `propensity()` is a propensity score, but those assumptions are not enforced tightly enough.

If invalid values slip through, the code may still run and return a number, which is worse than throwing an early error.

### What should be validated

#### Treatment indicator
Check that, on the estimation sample, `treat` only takes values `0` and `1`.

A standard pattern is:

```stata
quietly count if `touse' & !inlist(`treat', 0, 1)
if r(N) > 0 {
    di as error "treat() must be coded 0/1."
    exit 198
}
```

#### Propensity score range
Check that `propensity` is in `[0,1]` on the estimation sample:

```stata
quietly count if `touse' & (`propensity' < 0 | `propensity' > 1)
if r(N) > 0 {
    di as error "propensity() must lie in [0,1]."
    exit 198
}
```

#### Optional but useful
If you want to be stricter, warn or error when the propensity is numerically too close to 0 or 1 on many observations, since clipping can otherwise hide serious overlap problems.

### Recommended approach

At minimum, add hard validation for:

- treatment is binary,
- propensity is bounded in `[0,1]`.

That should happen before any pseudo-outcome is generated.

---

## 4. Clarify or redesign the Mahalanobis jackknife

### Why this is still a substantive issue

The current jackknife routine in Mata is now much better than before, but it is only truly exact if deletion affects the sample solely through rank removal and rank shifting.

That is true for a fixed metric like L2.

It is **not** automatically true for Mahalanobis distance, because the Mahalanobis metric depends on the sample covariance matrix. If one observation is deleted, the covariance matrix changes, so the metric changes too. In that case, the leave-one-out ranks are not generally obtained by a simple rank shift.

### Current mismatch

The code currently:

- computes Mahalanobis ranks from the full-sample covariance,
- then uses the rank-shift jackknife as if those ranks were stable under deletion.

That means the current jackknife is not a literal delete-1 recomputation under Mahalanobis distance unless the covariance matrix is intentionally held fixed.

### There are three viable ways to fix this

#### Option A: Document a fixed-metric jackknife
Define the Mahalanobis metric using the full-sample covariance matrix and hold it fixed across jackknife deletions.

If that is the intended design, then the current rank-shift logic is much more defensible, but this must be stated clearly in the documentation and comments.

This is the easiest path if speed is important.

#### Option B: Implement a true leave-one-out Mahalanobis jackknife
For each deleted observation:

- rebuild `X_{-i}`,
- recompute the leave-one-out covariance inverse,
- recompute Mahalanobis distances,
- recompute ranks,
- recompute weights,
- then recompute the leave-one-out estimate.

This is slower but conceptually exact.

#### Option C: Disallow jackknife SE for Mahalanobis until it is implemented correctly
This is the safest short-run move if you do not want to commit yet to either fixed-metric or true-LOO semantics.

### My recommendation

For the next revision, choose one of these explicitly. My preference would be:

- either **document the Mahalanobis jackknife as fixed-metric**, and say so clearly,
- or **disable Mahalanobis jackknife SE temporarily** until the exact delete-1 version is implemented.

The current middle position is not ideal because it risks claiming an exactness property that is only correct under L2 or under a fixed metric.

---

## 5. Finish packaging the plugin properly

### Why this matters

The plugin may work in a development environment, but the package manifest still does not include the plugin binary. That means a user installing the package through the normal packaging route will not necessarily receive the file needed for the fast path.

### What should be changed

In `rdistnn.pkg`, add the plugin binary to the file list.

At a minimum, if the installable binary is intended to be

```text
src/rdistnn.plugin
```

then it needs to appear in the package manifest.

If platform-specific binaries are part of the intended release process, that should also be reflected consistently in the package structure and installation instructions.

### Recommended follow-up

After updating the manifest, test a clean install in a fresh adopath-style environment rather than in the development tree.

That is the only reliable way to verify that the plugin path is actually complete.

---

## 6. Clean up documentation mismatches that can confuse users

These are lower priority than the estimator issues above, but they still should be fixed.

### `cate_dnn` help-file example typo

The help example currently uses:

```stata
m0(m0_true) m1(m0_true)
```

The second argument should presumably be:

```stata
m1(m1_true)
```

This is a small fix, but it matters because help examples are where users learn the intended API.

### README reference to missing setup file

The README refers to `dev-setup.do`, but that file is not present in the archive.

Either:

- add the missing file,
- or remove / replace the instruction.

Broken setup references make the development workflow harder to reproduce.

---

## 7. Supported Stata version is now encoded consistently

Status update, 2026-05-30: the active ado programs and test harness now target
Stata 17 via `version 17`, and the README/help files document Stata 17 or
later as the supported floor.

### Current situation

The intended support window, from prior discussion, is Stata 17, 18, and 19.
The code now uses the oldest version in that window as the version target.

### Remaining release validation

Before a release, test under Stata 17, 18, and 19 when those runtimes are
available.

### Why this matters

Using the oldest genuinely supported Stata version as the version target makes the package behavior more predictable and avoids accidental reliance on newer features that were never tested in older releases.

---

## Suggested order of work

If the goal is to make the next revision substantively correct as quickly as possible, I would proceed in this order:

1. **Correct the `cate_dnn` signal construction**
2. **Remove or redesign the single-nuisance fallback**
3. **Mark out missing outcomes in `cate_dnn` and `np_dnn`**
4. **Add hard validation for binary treatment and bounded propensity**
5. **Resolve the Mahalanobis jackknife semantics**
6. **Add the plugin to the package manifest**
7. **Fix help-file / README inconsistencies**
8. **Standardize version targeting and test across Stata 17/18/19** (source and docs now target Stata 17+; multi-version runtime validation remains a release check)

---

## Bottom line

The next revision should focus on three core items above all else:

- **make the `cate_dnn` score correct,**
- **protect the estimation sample from missing outcomes and invalid inputs,**
- **be honest and explicit about what the Mahalanobis jackknife is computing.**

Once those are fixed, the remaining tasks are packaging and documentation cleanup rather than estimator correctness.
