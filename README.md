# RdistNN_stata

Stata implementation of the DNN (distributional nearest neighbors) estimator
for conditional Z-estimation, companion to the R package
[RdistNN](../RdistNN) and the paper:

> **Robust Inference for Conditional Z-Estimation via Distributional Nearest
> Neighbors**  
> Jakob Juergens

---

## Commands

| Command | Description |
|---|---|
| `cate_dnn` | Conditional average treatment effect at a point x0 |
| `np_dnn` | Nonparametric regression E[Y\|X=x0] at a point |
| `weights_dnn` | Compute DNN rank weights (utility) |

All three commands require Stata 16 or later.

---

## Installation (development)

Clone or copy this repository, then run `dev-setup.do` once per Stata
session to add `src/` to the adopath:

```stata
do "path\to\RdistNN_stata\dev-setup.do"
```

To make the adopath addition permanent across sessions, add the
`adopath ++` line to your `profile.do`
(find it with `di c(sysdir_personal)`).

---

## Quick start

```stata
* Simulate data
clear
set obs 1000
set seed 42
gen x1 = rnormal()
gen x2 = rnormal()
gen e  = invlogit(0.5*x1 - 0.5*x2)
gen D  = (runiform() < e)
gen m0 = 1 + x1 - x2
gen m1 = m0 + 0.5 + 0.5*sin(x1)
gen Y  = m0 + D*(m1 - m0) + rnormal()*0.5

* Estimate CATE at the origin with oracle nuisances
cate_dnn Y, treat(D) x(x1 x2) x0(0 0) propensity(e) m0(m0) m1(m1)

* Inspect DNN weights
weights_dnn x1 x2, x0(0 0) k(5) generate(w)
summarize w
```

---

## Package structure

```
RdistNN_stata/
├── .gitignore
├── README.md
├── rdistnn.pkg          ← package metadata (for net install)
├── stata.toc            ← table of contents (for net install)
├── dev-setup.do         ← adds src/ to adopath for development
└── src/
    ├── _rdistnn_mata.ado    ← Mata core: distances, ranks, weights
    ├── cate_dnn.ado         ← CATE estimator
    ├── cate_dnn.sthlp       ← help for cate_dnn
    ├── np_dnn.ado           ← nonparametric regression
    ├── np_dnn.sthlp         ← help for np_dnn
    ├── weights_dnn.ado      ← weight utility
    └── weights_dnn.sthlp    ← help for weights_dnn
```

---

## Notes

- **Nuisances**: `cate_dnn` requires user-supplied propensity scores and
  outcome models.  These may be oracle values or first-stage estimates;
  Neyman orthogonality ensures robustness to first-stage bias at the local
  effective sample size scale.

- **Tuning**: `k` controls the number of nearest neighbors; `s` is the
  kernel order (bandwidth analog).  The default `s(0)` uses the
  undersmoothing rule `floor(n^((kx+1)/(kx+2)))` from the paper.

- **Mahalanobis distance**: use `distance(mahalanobis)` when covariates
  are strongly correlated (corresponds to `dnn_oracle_mahal` in the R
  simulation study).

- **Standard errors**: computed via the fast leave-one-out jackknife,
  which exploits the L-statistic structure of the DNN estimator.
