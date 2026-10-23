# Error reporting and testing guidance

This note records the reporting and testing contract for implementation and
maintenance. Shared messages, native boundaries and supported public notice
and reuse paths implement this contract. One-off validators can retain their
existing messages; changes to them must follow this guidance.
Use [AGENTS.md](../AGENTS.md) for editing rules,
[CONTRIBUTING.md](../CONTRIBUTING.md) for commands and checks, and
[the test map](../tests/README.md) for suite selection. The
[numerical reference](../docs/numerical-stability.md) owns arithmetic and module
responsibilities. This note keeps durable requirements rather than session
plans, execution receipts or historical-fit migration procedures.

## Shared ownership and parity

Keep input validation, saved-state validation, optional-extra access and
scientific statuses at their existing shared owners. A small common message
utility may centralize rendering and termination; do not create an independent
message framework for each estimator. Kernels return numerical facts/statuses;
public orchestration owns user context, warning frequency and summaries.

Shared R/Stata methods should agree on causes, consequences, recovery and
affected target identities. They use native argument syntax and host controls:
R condition objects/handlers and Stata return codes/quiet output are distinct
mechanisms. No new runtime dependency or parallel numeric error-code space is
needed. NP, TDNN, CATE and CLATE lead shared changes; PLQR remains experimental.

## Severity and actionable messages

| Situation | Required policy |
| --- | --- |
| Invalid input or infeasible main order | Early error naming the option, expected constraint and useful supplied value |
| Missing/inconsistent current retained state | Error identifying the state problem and requesting refitting |
| Valid fit lacks an optional extra needed by a requested feature | Distinct missing-extra error with the correct retention/refit action |
| Explicitly required native engine unavailable | Error describing the selector and available alternative |
| Automatic acceleration unavailable | Documented equivalent portable dispatch, without warning by default |
| Exact distance ties | One inference-qualified advisory per public call |
| Weak localization, ignored extras or infeasible optional candidates | One bounded warning per issue family; preserve applicable rows/reasons |
| Requested estimates/intervals unavailable in status-bearing output | Preserve rows, finite components and statuses; give one bounded summary |
| Negative first stage with a defined ratio | Descriptive note/status, not an automatic identification failure |
| Uncertainty intentionally disabled with nose | Description, not a warning about unavailable inference |
| Unexpected invariant, resource failure or Break | Propagate/stop, preserving the original code and cleanup |

State the problem, the expected/received fact and a useful next action. Name
actual options such as `k()`, `s()` and `x0matrix()`, rather than requiring users
to understand private Mata function names. For example: "Delete-one inference
requires s() < the estimation sample size; got s = 100 and n = 100. Reduce s(),
or specify nose."

Distinguish wrong type, dimensions, missing values, absent fields and semantic
conflicts. A fresh numerical-range failure needs rescaling guidance; it is not
an older-fit problem. Genuine internal errors need a reproducible-example
reporting action. Do not silently repair incomplete current state or add
compatibility readers for older package schemas.

Use complete sentences and consistent punctuation. Show at most five example
locations/names and the omitted count. Explicitly abbreviate individual names
or values longer than 80 characters. Aim for primary console messages under
800 characters, keeping the remedy and full-status destination visible. Do
not dump observations, signals, nuisance arrays or whole fitted states.
Treat dynamic text as literal data, including braces, quotes, percent signs
and control characters; static owned help markup may be separate.

Keep numerical/status contracts distinct. A paired ratio can validly be
unavailable at every target while retaining its component/status information;
an ordinary nonfinite mean need not follow the same return policy. Preserve
existing localization and denominator criteria. Signed TDNN concentration is
not probability-weight ESS. Consolidating messages does not justify new
thresholds, tie tolerances, jitter or statistical guarantees.

## Stata, Mata and plugin signaling

Use standard return codes with explanatory text: for example, 198 for invalid
syntax/options, 301 for absent estimation results, and 498 for package-specific
state/statistical problems. Preserve meaningful file, memory and Break codes;
plugin allocation failure should not masquerade as invalid syntax. Do not
replace `e()` or `r()` with a transient last-error store.

Fatal ado errors use `display as error` and exit with the selected code. Mata
fatal helpers choose the error display level, render dynamic text literally
and terminate; ordinary user errors should not require a raw Mata traceback.
A custom explanation does not need a second generic "invalid syntax" message.
Keep plugin `SF_error()` for failures, not routine nonfatal advisories, and
contain every native exception at the ABI with allocation-safe failure handling.

Nonfatal warnings/notes use text output with an explicit `warning:` or `note:`
prefix and respect `quietly`. Reserve error output for errors: it remains
visible under `quietly`, whereas `capture` may suppress it. Save `_rc` before
cleanup or another command can overwrite it; `capture noisily` preserves
output while allowing handling. Propagate Break (1) and genuine resource errors.

Recover only from specifically understood failures. Missing acceleration can
select the portable implementation; an unexpected native failure cannot simply
become "engine unavailable". A recognized optional-candidate failure may retain
a failed row, but the requested main result must not quietly disappear.
Do not choose recovery by parsing English text. The R companion must also
preserve user warning handlers and warnings promoted by `options(warn = 2)`.

## Reuse grouping information and bound reporting cost

Prepare selected-metric groups once for unchanged geometry, retaining derived
tie information from that work. Pass the existing groups/facts through fitting,
grids, prediction, diagnostics and plot orchestration. Reporting must not
recompute distances, re-sort, reconstruct groups or repeatedly inspect the
same group columns for each consumer. Current fits already retain
`e(model_groups)`; validate complete current state before deriving/reusing
call-local reporting facts. Do not trust an unchecked persisted tie flag.

Deletion-refitted geometry needs its own actual groups. Frozen geometry still
uses the true surviving block counts in numerical calculations. Reuse does
not mean carrying full-sample facts into a different geometry or changing the
distance definition. Preserve exact represented equality and signed group mass.

Emit at most one notice per recurring issue family per public call, including
nested candidates and plot consumers. Independent calls report independently.
Do not add a session-global or fitted-object "already warned" flag. Complete
counts come from the owning result/status batch, not overlapping panel counts.
Preserve the qualification that the theory assumes a continuous distance
distribution and does not establish inference under persistent ties.

Construct text only for output that will be emitted. Reuse existing diagnostic
passes and bound examples; output should scale with issue families rather than
observations or target count. Gate expensive display-only work under quiet
execution. Keep reporting context out of scientific geometry and numerical
inner loops. Never remove safety validation or numerical guards to save
reporting work, or replace an original failure by flushing pending warnings
during cleanup.

## Unit and regression tests

New behavior needs valid, invalid and boundary cases in its owning suites.
Each real regression should fail for the original defect at the most stable
public boundary available. Reuse established fixtures, independent references
and caller-state comparators; do not add copied corruption runners or one-off
frameworks. Register suites once and update family filters and test maps.

Review a small set of complete canonical messages, then test consumer return
codes, meaningful problem/remedy text, affected identities and notice counts.
For R structured conditions, prefer class/code and useful fields over repeated
English snapshots. Wording must not become a recovery protocol.

Cover ordinary, quiet, capture and capture-noisily execution; exact return codes;
cold loading with cleared data; literal strings; long/bounded detail lists; and
promised `e()`/`r()` and caller-data preservation on success and failure. Test
independent calls and nested consumers. Deliberately disabled uncertainty must
not warn. Partial/all-unavailable ratios retain their original rows/statuses.
Test recognized recovery separately from invalid explicit options, corrupted
state, unexpected errors, Break and allocation failure. Inject unrelated
warnings/errors to prove they remain observable rather than broadly suppressed.

Prove grouping reuse with counters or fault injection: after preparation,
reporting causes zero distance evaluations, sorts or group rebuilds. A cache
field existing is insufficient evidence. Exercise scalar/batch calls, repeated
and reordered targets, grids, restored predictions and plots. Include a
deletion-refitted tie pattern differing from the full sample and corruption
rejected before cached facts are used.

Retain exact rank, support, group, target, status and missingness checks; use
justified tolerances for continuous quantities and independent mathematical
or explicit-deletion references. Force Mata/plugin engines and prove actual
execution; test automatic capability fallback separately from real execution
failure. Do not enshrine avoidable overflow as the correct answer to an exactly
finite calculation. Guard tests need independently genuine overflow and finite
controls, while permitting valid numerical improvements.

Measure reporting cost/allocations on ordinary and tied scalar, batch, grid
and reuse routes before/after implementation. Keep structural tests against
repeated grouping, per-observation condition allocation and unbounded output.
Timings are review evidence, not flaky wall-clock CI thresholds or permission
to weaken mathematical checks. Use [CONTRIBUTING.md](../CONTRIBUTING.md#verification)
for final verification appropriate to the changed source.

## Reference basis

See the official [Stata error codes](https://www.stata.com/manuals/perror.pdf),
[display conventions](https://www.stata.com/manuals/pdisplay.pdf),
[capture behavior](https://www.stata.com/manuals/pcapture.pdf),
[Mata errors](https://www.stata.com/manuals/m-5error.pdf) and
[literal Mata display](https://www.stata.com/manuals/m-5display.pdf).
The companion design uses [base R conditions](https://stat.ethz.ch/R-manual/R-patched/library/base/html/conditions.html)
and [R warning protocols](https://stat.ethz.ch/R-manual/R-patched/library/base/html/warning.html).
