# Coherent estimator reporting and postestimation

The user requests standardization of summaries, diagnostics and plots across
CATE, CLATE and nonparametric regression, followed by a commit and push to the
existing `feature/clate` branch. Baseline: `23a0736`. RdistNN remains the
numerical reference; this is a Stata presentation/interface change, preserving
the estimators and their scientific distinctions.

## Common contract

1. CATE, CLATE and NP use one summary adapter and common metadata header:
   sample size, target/covariate count, covariates, method, metric, k/s,
   ESS or ESS range, tuning and covariance method. Fit and replay return a
   conventional 9-by-m `r(table)` (b,se,z,pvalue,ll,ul,df,crit,eform) plus
   `r(summary)` with m rows and six columns estimate,std_error,z,pvalue,
   ci_lower,ci_upper. The full target order is preserved even if CLATE omits
   undefined coefficients from e(b). Replay level changes bounds without
   changing e() or the fit. PLQR remains a supplied-coefficient moment
   evaluator; it can share the metadata header without claiming coefficient
   inference.
2. Add canonical `e(estimates)`, `e(se_targets)` and `e(estimate_vcov)` to
   all three scalar estimators; keep existing coefficient names, matrices and
   aliases. `nose` means missing uncertainty, never zero precision. Preserve
   R's CLATE status masks and its distinction between zero ratio variance and
   ordinary mean/CATE zero variance; no fabricated weak-IV inference.
3. Every diagnostics command supports the same targets(), candidates(),
   allcandidates, detail, saving(), replace grammar; target and candidate
   selection are mutually exclusive. Return common `r(localization)` fields
   target_id, target coordinates, and the first eleven shared weight/distance
   diagnostics, preserving existing r(diagnostics)/r(compact) and scientific
   extensions. Preserve identities/order, e() and data. Shared helpers own
   selection, candidate metadata and saving; adapters own field names and extras.
4. Plot commands share validation, naming, exporting and return conventions.
   All single-view plots expose selected `r(plot_data)`, with consistently
   named fields for comparable views. Return `r(n_rows)` for data rows while
   retaining legacy `r(n_plotted)` meanings; dashboard composition instead
   returns panel metadata. Estimate/tuning tables and bounds use the common
   reporting data. Support target/candidate selections wherever meaningful;
   reject unsupported combinations instead of silently ignoring options.
5. CATE, CLATE and NP observation-level kernel plots, and effect influence
   plots, verify original values/sample/order with one shared signature helper.
   Stored-result plots continue working after data changes. Confidence bands
   must render for all applicable estimators, including regression outcomes.
6. Keep CLATE instrument overlap, paired score influence, components and
   denominator statuses; keep CATE treatment overlap; keep NP outcome
   distribution views. Public family entry points and installed help remain
   independently usable, including cold loading and estimates store/restore.

## Acceptance

Cross-estimator regressions cover positive and unavailable uncertainty,
zero variance, undefined/partial CLATE targets, batch identities, replay level,
legacy result restoration, ordinary r(table), diagnostics/plot selections and
saved data, incompatible options, changed/reordered samples and cold loading.
Add actual R reporting fixtures where applicable. Retain existing numerical
goldens and run the full licensed Stata suite, static checks, Python
maintenance, isolated install and source archive checks. Inspect rendered
exports and obtain an independent final review before commit/push.
