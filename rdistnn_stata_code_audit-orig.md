# Audit note on the Stata/Mata layer in `RdistNN_stata`

I reviewed the Stata and Mata sources in the archive, focusing on the `.ado` files and `_rdistnn_mata.ado`. I did **not** execute the code in Stata in this environment, so the note below is a **static code audit**. The issues flagged here are grounded in the source and, in several cases, look like clear execution bugs rather than stylistic concerns.

## Bottom line

There are several real problems on the Stata/Mata side.

The most important ones are:

1. **`distance(mahalanobis)` is advertised but not actually implemented.** In practice the code always goes through the L2 path.
2. **The helper program name `Display` is reused across multiple `.ado` files**, which is dangerous in Stata because program names live in a global namespace.
3. **`weights_dnn` has a fallback bug when `generate()` is used and the plugin path fails**: it can try to `gen` the same variable twice.
4. **`weights_dnn` incorrectly requires at least two covariates** via `varlist(min=2)`.
5. **`weights_dnn` does not `markout` missing covariates**, so missing `x` values can leak into the distance/rank calculation.
6. **The reported “LOO jackknife” does not appear to be a literal delete-1 recomputation of DNN weights**, but rather a fixed-weight self-normalization formula. That may be intentional as an approximation or shortcut, but it is not the ordinary delete-1 jackknife for a rank-weight estimator.

I would treat items 1, 2, 3, and 6 as the most important.

---

## Detailed findings

### 1. `distance(mahalanobis)` is not implemented, despite being documented and accepted by the interface

**Where this shows up**

- `cate_dnn.ado` accepts and validates `distance(l2|mahalanobis)` at lines 57–60.
- `weights_dnn.ado` accepts and validates it at lines 47–50.
- `np_dnn.ado` exposes the option in `syntax` at line 19 and sets a default at line 41.
- The help files and README explicitly describe Mahalanobis support.

**Why this is a problem**

The computation path never uses a Mahalanobis distance:

- In `cate_dnn.ado`, the plugin call is always
  `weights_l2 ...` at lines 127–128.
- In `np_dnn.ado`, the plugin call is always
  `weights_l2 ...` at lines 67–68.
- In `weights_dnn.ado`, the plugin call is always
  `weights_l2 ...` at lines 79–80 and 85–86.
- In `_rdistnn_mata.ado`, the only distance routine is `rdistnn_l2dist()` at lines 61–66.
- The top-level Mata wrapper `rdistnn_weights()` always calls `rdistnn_l2dist()` at lines 217 and 221.

So at the moment, **`distance(mahalanobis)` is silently ignored**. The user may believe they are getting Mahalanobis weighting while the code still computes Euclidean/L2 weights.

**Severity**

High. This is not a cosmetic issue. It can materially change the estimator the user thinks they ran.

**Suggested fix**

Either:

- **remove/disable the option** until it is truly implemented, or
- add a genuine Mahalanobis path in both Mata and the plugin, and only dispatch to the plugin fast path when `distance(l2)` is requested.

If you want a clean interim patch, I would strongly prefer:

- reject `distance(mahalanobis)` with an error for now, rather than accept it silently.

---

### 2. The helper program name `Display` is reused across `.ado` files

**Where this shows up**

- `cate_dnn.ado`: `program define Display` at lines 177–188
- `np_dnn.ado`: `program define Display` at lines 114–123

Both main programs also call `Display` directly:

- `cate_dnn.ado` lines 10 and 173
- `np_dnn.ado` lines 9 and 111

**Why this is a problem**

Stata program names are global in the session. Reusing a generic helper name like `Display` across package files is risky:

- one file can overwrite or collide with another,
- loading one command can affect replay/output behavior of another,
- in the worst case, a second `program define Display` can fail because the name already exists.

Even if it happens to work in some load orders, it is fragile and unnecessary.

**Severity**

High. This is a classic namespace collision bug in Stata package code.

**Suggested fix**

Rename the helpers to unique internal names, e.g.

- `cate_dnn_display`
- `np_dnn_display`

and call those explicitly.

---

### 3. `weights_dnn` can fail in the fallback path when `generate()` is used

**Where this shows up**

In `weights_dnn.ado`:

- plugin branch with `generate()` pre-creates the variable at line 78,
- then attempts the plugin call at lines 79–80,
- if the plugin call fails, `__plugin_ok` becomes 0 at line 88,
- fallback branch then tries to `gen double \`generate\' = .` again at line 112.

**Why this is a problem**

If the plugin path is attempted and fails, and the user specified `generate()`, the fallback code tries to create a variable that already exists.

That means a plugin failure can cascade into a secondary Stata error like “variable already defined”, even though the pure-Mata fallback should have rescued the call.

This is especially relevant because the plugin path is exactly where unsupported options or deployment issues can surface.

**Severity**

High. This is a concrete control-flow bug.

**Suggested fix**

Do one of the following:

- only create `generate()` **after** the plugin has succeeded, or
- if the variable was already created before the plugin attempt, use `replace`/`st_store()` in the fallback instead of `gen`, or
- always create a tempvar for the plugin attempt and copy into `generate()` only after success.

That last option is the cleanest.

---

### 4. `weights_dnn` incorrectly requires at least two covariates

**Where this shows up**

- `weights_dnn.ado`, line 7: `syntax varlist(min=2)`

**Why this is a problem**

This blocks one-dimensional use of `weights_dnn`, even though 1D DNN weights are perfectly meaningful and your other commands do not require `kx >= 2`.

This looks accidental rather than intentional.

**Severity**

Medium-high.

**Suggested fix**

Change

```stata
syntax varlist(min=2) ...
```

to

```stata
syntax varlist(min=1 numeric) ...
```

or similar.

---

### 5. `weights_dnn` does not remove observations with missing covariates

**Where this shows up**

- `weights_dnn.ado` uses `marksample touse` at line 18,
- but unlike the other estimators it does **not** call `markout` on the covariates.

Compare with:

- `cate_dnn.ado` line 32: `markout \`touse\' \`treat\' \`x\' \`propensity\' \`m0\' \`m1\'`
- `np_dnn.ado` line 25: `markout \`touse\' \`x\'`

**Why this is a problem**

As written, `weights_dnn` can leave missing covariate values inside the estimation sample. Those missings then flow into Mata or the plugin distance computation.

At best, missing-distance observations get dumped somewhere at the end of the rank ordering; at worst, they distort or break the intended logic. In any case, that behavior is not what users expect.

**Severity**

High for robustness, especially in practical use.

**Suggested fix**

Add something like

```stata
markout `touse' `xvars'
```

after parsing the varlist.

---

### 6. The current “LOO jackknife” is not a literal delete-1 recomputation of the DNN estimator

**Where this shows up**

- `cate_dnn.ado`, lines 136–146
- `np_dnn.ado`, lines 76–86

The leave-one-out value is computed as

```stata
(tau_hat - w_i * rho_i) / (1 - w_i)
```

or analogously for `np_dnn`.

**Why this is a problem**

That formula is exact for a weighted mean when the weights are treated as fixed and one observation is removed with renormalization.

But DNN weights are **not fixed** in the sample. They depend on:

- the full rank ordering,
- the sample size `n`,
- and, after deletion, the ranks of the remaining observations.

So the true delete-1 estimator would ordinarily require recomputing the DNN weights on the sample with observation `i` removed.

The current code does **not** do that. It keeps the original weights and applies a self-normalized deletion formula.

That means one of two things must be true:

1. this is only an approximation / shortcut, or
2. there is a nontrivial theorem showing that for this particular L-statistic representation the shortcut equals the intended delete-1 quantity.

From the code alone, I do **not** see such a theorem being encoded. So as a code audit finding, I would say the current implementation is **not the ordinary delete-1 jackknife**.

**Severity**

High, but different in character from the syntax bugs. This is more of a statistical/algorithmic concern than a parser/runtime bug.

**Suggested fix**

At minimum, clarify the documentation:

- if this is an approximation, say so explicitly;
- if it is meant to be exact, document the algebra showing why.

If you want the literal delete-1 jackknife, you would need to recompute the weights on the leave-one-out sample.

---

### 7. `np_dnn` exposes `distance()` but does not really validate or use it

**Where this shows up**

- `np_dnn.ado` includes `Distance(string)` in `syntax` at line 19.
- It sets the default at line 41.
- But there is no `inlist()` validation analogous to the other two commands.
- More importantly, the option is never used in the actual computation path.

**Why this matters**

This makes `np_dnn` especially brittle:

- arbitrary strings can slip through,
- and even valid strings do nothing.

Also, unlike `cate_dnn`, it does not store `e(distance)`.

**Severity**

Medium.

**Suggested fix**

Either fully implement the distance dispatch or reject unsupported values now. Also store the chosen distance in `e()` if the option remains part of the public interface.

---

### 8. `np_dnn` does not explicitly check `k >= 1`

**Where this shows up**

- `np_dnn.ado`: there is no `k` validation comparable to `cate_dnn` line 52 or `weights_dnn` line 24.

**Why this matters**

A user can in principle pass `k(0)` or something nonsensical. The downstream code will then fail later and less transparently.

**Severity**

Medium.

**Suggested fix**

Add the same user-facing validation used in the other commands.

---

### 9. No explicit `s <= n` validation on the Stata side

**Where this shows up**

- `cate_dnn.ado`, `np_dnn.ado`, and `weights_dnn.ado` only check `s < k`, not `s > n`.

**Why this matters**

Conceptually the DNN subset size should satisfy `k <= s <= n`.

In the Mata exact path, `_rdistnn_mata.ado` does check this at lines 120–123. But:

- the Stata side does not catch it early,
- the plugin side does not appear to guard against `s > n` before computing,
- and the asymptotic path may still run with an impossible `s` value.

**Severity**

Medium.

**Suggested fix**

Add explicit Stata-side validation in all three entry points:

```stata
if `s' > `n' {
    di as error "s(`s') exceeds the sample size n(`n')."
    exit 198
}
```

---

## Secondary/documentation/package issues

### 10. Documentation currently overclaims what the code does

This is mostly downstream of item 1 and item 6, but it is worth stating explicitly.

- The help files and README present `distance(mahalanobis)` as available.
- The help files and README present the SEs as coming from a “fast leave-one-out jackknife”.

Given the current code, both statements are at least misleading:

- Mahalanobis is not implemented.
- The jackknife is not a literal recomputation-based delete-1 jackknife.

I would adjust the documentation immediately, even before changing the code.

---

### 11. The package metadata does not include the plugin binary

**Where this shows up**

- `rdistnn.pkg` lists the `.ado`/`.sthlp`/Mata files,
- but it does **not** list `src/rdistnn.plugin`.

**Why this matters**

If this archive is what you intend to distribute via Stata’s package mechanism, the plugin will not get installed from the package metadata as written.

That is not a C++ logic bug, but it does affect the Stata-side deployment story.

---

### 12. `README.md` refers to `dev-setup.do`, but it is not in the archive

The README says to run `dev-setup.do`, but I do not see that file in the archive listing.

Again, not a core estimator bug, but it is a concrete development/setup inconsistency.

---

## What I would fix first

If I were patching this package in order of priority, I would do:

1. **Remove or hard-disable `distance(mahalanobis)` until it is actually implemented.**
2. **Rename the helper `Display` programs to unique names.**
3. **Fix `weights_dnn` fallback logic when `generate()` is used.**
4. **Add `markout` in `weights_dnn`.**
5. **Relax `weights_dnn` to allow one covariate.**
6. **Add uniform validation for `k >= 1` and `s <= n`.**
7. **Revisit the jackknife implementation and either document it carefully or replace it with a true delete-1 recomputation if that is what the theory requires.**

---

## Short verdict

Yes, I do think there are real problems in the Stata/Mata layer.

The biggest concrete implementation bug is the **silent non-implementation of `distance(mahalanobis)`**. The biggest code-structure bug is the **shared `Display` program name**. The biggest substantive concern is the **jackknife shortcut**, which does not look like a true leave-one-out recomputation for DNN rank weights.

