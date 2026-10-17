# Estimator coherence implementation plan

> **For agentic workers:** Use superpowers:executing-plans for the parent implementation; independent diagnostics and plotting tasks have separate ownership and may use superpowers:dispatching-parallel-agents. Track progress in ignored build/coherence/progress.md.

**Goal:** Standardize applicable summary, diagnostics and plotting interfaces and push the verified changes to feature/clate.

**Architecture:** Common reporting/state helpers plus thin scientific adapters; shared selection/dispatch/export paths replace independent policy decisions. Preserve legacy fields and estimator-specific interpretation.

**Tech Stack:** Stata 17 ado/Mata, licensed Windows StataNow19, R reference fixtures, Python/shell checks.

**Spec:** `docs/superpowers/specs/2026-10-17-estimator-coherence-design.md`

## Global Constraints

- Remain on feature/clate; no default-branch merge or release.
- R remains the numerical reference; no sibling edits or estimator changes.
- Preserve legacy coefficient names and diagnostic schemas, with additive common outputs.
- Keep partial/undefined targets visible in full stored/summary tables (plot data describe finite rows actually drawn); never replace unavailable estimates/uncertainty by zero.
- PLQR remains a supplied-coefficient moment evaluator.

## Review Focus

- r(table) must contain actual full-target statistics after fit/replay, not a missing scalar accepted as a matrix.
- nose, zero SE and CLATE status masking must remain correct without stale r() results.
- Selected targets/candidates must match displayed, returned and saved data.
- Data guards must reject reordered/mutated observation-level inputs across all three estimators.
- Shared helper loading must work from clean installation and restored estimates.

### Task 1: Parent-owned reporting and state

**Files:** New `_rdistnn_report.ado`, `_rdistnn_estimation_state.ado`, `_rdistnn_sample_signature.ado`; modify three scalar estimator ado display paths and existing CLATE signature alias; add `tests/test-estimator-reporting.do` and R reporting fixtures/generator. Optionally share the header with PLQR while retaining moment-only display.

**Interfaces:** `_rdistnn_report [, level(#)]` rclass returns table/summary; `_rdistnn_estimation_state varlist` eclass attaches canonical vectors/covariance and ordinary/ordered signatures. `_rdistnn_sample_signature varlist [if] [in]` rclass returns signature; original CLATE helper remains an alias. Canonical fields and e(signature_vars)/e(ordered_signature)/e(datasignaturevars) become common.

- [x] Write and observe failing tests for actual r(table), six-column summary, replay level, metadata and exceptional targets.
- [x] Implement shared report/state helpers and adapters, preserving scientific masks and legacy results.
- [x] Compare reported numerical fields to R reference fixtures; run focused licensed tests.

### Task 2: Diagnostics worker

**Files:** Three diagnostics ado/help files, new shared diagnostic helper(s), `tests/test-estimator-diagnostics.do`.

**Interfaces:** Common targets/candidate grammar and r(localization) as in spec. Existing full/compact matrices and CLATE component/status tables remain available. Parent handles manifests/runner/docs.

- [x] Observe failing cross-family target/return/saving tests.
- [x] Share selection, candidate assembly and saving policy; add common localization output.
- [x] Run focused tests for selection identity, incompatible options, e()/data preservation and family-specific extras.

### Task 3: Plotting worker

**Files:** Three plot ado/help files, shared plot helpers/dispatcher, `tests/test-estimator-plots.do`.

**Interfaces:** Read canonical reporting fields when present, preserve legacy fallback. Use `_rdistnn_sample_signature e(signature_vars) if e(sample)` plus ordinary signature for observation-level views. Expose selected r(plot_data), r(n_rows), applicable target/candidate and dashboard metadata, with shared export and validation.

- [x] Observe failing cross-family plot-data/selection/guard/render tests.
- [x] Share dispatcher/export policy and builders; align compatible options and reject no-op combinations.
- [x] Render profile confidence bands for CATE/CLATE/NP and NP outcomes; run focused and cold-loading regressions.

### Task 4: Parent-owned acceptance and publication

**Files:** Manifest, runner, install/static checks, docs/help for estimator reporting, README, release notes and verification receipt.

- [x] Document the common interfaces and scientific exceptions; install all new helpers.
- [x] Run full licensed runtime, fixed R references, static/Python/archive checks and isolated install.
- [x] Obtain final independent review, fix material findings with regression tests, inspect staged diff/check.
- [x] Commit and push feature/clate.


Implementation and acceptance evidence: [verification receipt](../../estimator-coherence-verification.md).
Publication uses the existing feature branch as requested; no default-branch merge.
