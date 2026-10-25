# Error reporting and testing guidance

This note records the reporting and testing contract for implementation and
maintenance. Package-owned custom errors use the existing ado and Mata message
owners throughout validation, fitting and reuse. Native validation and runtime
errors retain their host behavior; the boundaries below define their roles.
Use [AGENTS.md](../AGENTS.md) for editing rules,
[CONTRIBUTING.md](../CONTRIBUTING.md) for commands and checks, and
[the test map](../tests/README.md) for suite selection. Record newly adopted
principles in the maintained note that owns them, with authoritative sources,
the applicable boundary and regression expectations. Keep execution receipts
private; an intended policy alone does not establish implemented behavior. The
[numerical reference](../docs/numerical-stability.md) owns arithmetic and module
responsibilities. This note keeps durable requirements rather than session
plans, execution receipts or historical-fit migration procedures.

## Shared ownership and parity

Keep input validation, saved-state validation, optional-extra access and
scientific statuses at their existing shared owners. A small common message
utility may centralize rendering and termination; do not create an independent
message framework for each estimator. Kernels return numerical facts/statuses;
public orchestration owns user context, warning frequency and summaries.

Mata calculations used by ado commands run in typed functions with local
scratch. Top-level assignments create external globals, so neither estimation
nor plotting may retain training arrays or overwrite a caller binding. Dropping
or renaming a fixed scratch global is not a repair. Fresh-command IO adapters
capture validated ado locals at entry and publish only temporary Stata outputs;
the base adapter module loads after its numerical dependencies, while estimator
extensions keep their own adapters. Test existing and absent caller bindings
and the external-array inventory on success and after genuine public errors.
See Stata's [ado/Mata boundary](https://www.stata.com/manuals/m-1ado.pdf),
[typed local declarations](https://www.stata.com/manuals/m-2declarations.pdf) and
[external globals](https://www.stata.com/manuals/m-5findexternal.pdf).

Shared R/Stata methods should agree on causes, consequences, recovery and
affected target identities. They use native argument syntax and host controls:
R condition objects/handlers and Stata return codes/quiet output are distinct
mechanisms. No new runtime dependency or parallel numeric error-code space is
needed. NP, TDNN, CATE and CLATE lead shared changes; PLQR remains experimental.

Give each returned count one current unit. Plot data use `n_rows`, omissions
use `n_omitted`, and dashboards use `n_panels` plus `panel_id n_rows` metadata.
Do not preserve mixed-unit historical counts beside those owners. Target and
candidate identities describe the selected content directly; no older-output
alias is needed. Private plot builders publish the same actual-row contract,
so dashboard metadata can read counts without copying full plot matrices.

Capture scalar metadata before `return matrix` or `ereturn matrix` transfers
ownership: by default these commands destructively rename the source matrix.
Use `copy` only when the original must remain available, rather than retaining
every matrix as a workaround. See the [returned results manual](https://www.stata.com/manuals/preturn.pdf).

## Package version identity

`rdistnn.pkg`'s single `d Version:` declaration owns the package version. The
README, first release-note heading, overview help and every installed `.ado`
and `.sthlp` source stamp must agree with it, including the full development
suffix. Private ado files containing Mata definitions follow the same rule.
Use exactly one `*! version VERSION  YYYY-MM-DD` comment at the start of each
ado file; help keeps `{smcl}` first and wraps that stamp in a SMCL comment on
the second line. Other descriptive `*!` notes may remain. Do not maintain
independent per-command revision counters; Git records individual histories.

This is a package policy. Stata's [which manual](https://www.stata.com/manuals/rwhich.pdf)
describes `*!` lines as notes reporting file revisions, rather than executable
interpreter controls. Its [version manual](https://www.stata.com/manuals/pversion.pdf)
describes the separate command-language setting. Keep executable `version 17`
statements, saved-model schema numbers, native API versions and `v 3` package-file
syntax in their own roles. Neither the package version nor a comment stamp alone
establishes executable provenance; retain commit/archive and installed identities.

Windows `net install` can convert text payloads from LF to CRLF. Record raw
archive/source and installed hashes; when comparing ado/help content, permit
only that explicit line-ending conversion and verify the complete text. Native
binary payloads still require byte equality. A changed text hash alone does
not establish a different implementation, and a matching version alone does
not establish matching content.

The existing package-hygiene owner supplies both stamp validation and the
explicit `--sync-version-stamps` operation. Read the version from the manifest;
never hardcode a second current version in a checker. Synchronization preserves
executable bytes, line endings, authors and other comments, and rejects ambiguous
headers before writing any file. Its ISO date records synchronization, separately
from the package's publication date. Run it after each development or production
version change; use [RELEASING.md](../RELEASING.md) for the commands.

CI must reject stale, missing, duplicate, misplaced and malformed stamps in public
commands, private helpers and help files. Source archives apply the same shared
validator to their own manifest and payload, without relying on checkout versions.
Regression controls cover production and development versions, read-only checking,
idempotent synchronization, preserved executable/schema declarations and unchanged
historical records/pins. For a comment-only correction, compare executable bodies
and check metadata visibility in Stata; reserve numerical-suite reruns for changes
to executable inputs.

## Current user documentation

During active experimental development, describe current inputs and outputs
and briefly tell users to refit saved models after updates. Release notes own
retired options and schemas; current tutorials do not need migration catalogues.
Preserve published replication pins and current saved-state/reuse requirements.

Use the existing division between a first-fit walkthrough, task guides and full
installed command help. Give long shared explanations one primary home, while
retaining short interpretation limits beside outputs and complete local syntax,
result shapes and command-specific constraints in help. Source deduplication is
not by itself a shorter reading path. This applies the companion's
[documentation-role policy](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/11_examples_vignettes_and_pkgdown_policy.md),
based on [R Packages' workflow guidance](https://r-pkgs.org/vignettes.html)
and [Diataxis' documentation roles](https://diataxis.fr/), to Stata's installed
help and runnable do-files.

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

Use Stata's established return codes for the actual cause. Its
[error manual](https://www.stata.com/manuals/perror.pdf) groups syntax, stored
results, statistical, matrix, file, resource and Mata failures. Do not invent a
package code registry or replace `e()`/`r()` with a transient last-error store.

| Cause | Code and boundary |
| --- | --- |
| Invalid option combination, range, rank contents or finite-value requirement | 198, with a package explanation when native syntax does not express the constraint |
| Native type or existence validation | Let `syntax` and `confirm` report their own code; for example, integer confirmation can return 7 and a missing named matrix returns 111 |
| No active estimates for a requested operation | 301 |
| Incomplete current fitted state, unavailable retained feature or package-specific statistical/numerical limitation | 498 with an explanation and appropriate refit, retention or rescaling action |
| Purely incompatible matrix dimensions | 503; distinguish shape from malformed numlist/option contents |
| Empty usable sample or insufficient observations | 2000 or 2001 respectively when that is the actual cause |
| Break, file, allocation or unexpected runtime failure | Preserve the original code, including 1, file codes, 909 and native Mata runtime codes |

Native [confirm](https://www.stata.com/manuals/pconfirm.pdf) already supplies
messages and codes for ordinary verification. Do not capture and relabel those
as generic 198. A custom guard can use a standard code with more context; a
custom explanation needs no second generic message. A bare native `error 198`,
`error 301` or `error 2000` remains appropriate when its standard message is
sufficient. A combined malformed-request guard is not automatically a pure
matrix-conformability failure.

All custom ado fatal messages go through `_rdistnn_messages error`, which owns
`display as error _asis` and termination. All custom Mata guards go through
`rdistnn_fail(rc, message)`, which selects error output, calls literal
`display(message, 1)` and exits with the chosen code. Do not introduce parallel
`display as error` or `errprintf()` sites. Native Mata errors and invariants
retain their native signaling. The separate ado/Mata owners reflect host
language boundaries, rather than an incomplete migration. See the Mata
[error](https://www.stata.com/manuals/m-5error.pdf),
[display](https://www.stata.com/manuals/m-5display.pdf) and
[errprintf](https://www.stata.com/manuals/m-5errprintf.pdf) manuals.

There are exactly two direct ado rendering roles: the message owner itself
and the core loader's pre-helper plugin-declaration failure. The latter prints
a fixed literal explanation and its original code, then exits without trying
to load another helper under a possible resource failure. No blanket exception
exists for historical or one-off validators. C++ keeps native `SF_error()` for
failures and contains exceptions at the ABI with allocation-safe handling;
routine advisories belong to public orchestration.

Both renderers normalize ASCII control characters and preserve dynamic text
literally. The ado owner uses native Unicode string replacement without a Mata
dependency, so cold validation works before numerical modules load. Bound
individual dynamic values at their existing context owner rather than
truncating an entire remedy. Use compound quotes when passing dynamic text.

Nonfatal warnings/notes use text output with an explicit `warning:` or `note:`
prefix and respect `quietly`. Reserve error output for errors: it remains
visible under `quietly`, whereas `capture` may suppress it. Save `_rc` before
cleanup or another command can overwrite it; `capture noisily` preserves
output while allowing handling. Propagate Break (1) and genuine resource errors.
These are native [display](https://www.stata.com/manuals/pdisplay.pdf) and
[capture](https://www.stata.com/manuals/pcapture.pdf) semantics.

Program and matrix presence probes recognize only native missing-object 111;
Confidence-bound presence assertions recognize only success 0 or false assertion 9. Save and
propagate every other result before another command overwrites `_rc`. Keep the
explicit existing loader for multi-program plot files: ordinary ado autoload
does not make their auxiliary programs available. Do not replace it with a
new general loader or treat every probe failure as an absent helper.

Interpret recovery codes at the operation that generated them. Optional plugin
declaration recognizes missing/nonloadable binaries at that dynamic-load
boundary; it propagates unexpected declaration failures. Its loader mapping
must not be reused for plugin execution or matrix copies. Likewise, 198 denotes
an unsupported capability only at the documented capability probe, not at every
native call. See [plugin loading](https://www.stata.com/manuals/pplugin.pdf).

Check optional returned-matrix presence through `r(matrices)` before copying.
Copying an absent `r(geometry)` can succeed with a missing `1 x 1` placeholder,
so the assignment return code is not a presence test. Once a return is present,
copy failures propagate; do not turn allocation or Break into optional absence.
For ordinary temporary matrix names, `confirm matrix` can check existence;
inspect its saved result before output generation or returned-result commands.
Keep complete current fitted-state validation at the shared report/reuse owners.

The existing source/package-hygiene check enforces custom-renderer ownership
while permitting native validation and standard signals. It is an architecture
guard, not a Stata parser or runtime proof. Runtime regressions exercise cold
commands, literal quoted/control text, actual native codes, quiet/capture
channels and caller-state preservation. Use the existing message and optional
boundary suites, and inject only the host operation when testing unsafe resource
or Break outcomes; do not add production fault switches or phrase-locked tests.

At native bridges, retrieve or decode output only after a successful native
return. On failure, clear owned temporary matrices and propagate the original
code before attempting another large allocation or output copy. Test a failing
output retrieval after an injected native resource error or Break; a successful
copy in a small error fixture does not establish preservation under memory
pressure. These are review requirements for existing and new bridge paths.

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

Validate one unchanged snapshot once at the trusted public boundary, then pass
that call-local verification explicitly through private reporting adapters.
Standalone private adapters validate by default. Public syntax must not offer
a validation bypass, and no verification flag belongs in `e()` or a session
global. A derived post is a different snapshot: validate it before reporting
even when its parent was already checked. Separate public calls validate
independently, including after restored or changed fitted state.

Candidate grids belong to that snapshot. Check their dimensions, target and
computed-order identities, default anchor, aliases, feasibility, availability
and counts before any table is sliced. Preserve legitimate unavailable ratio
rows and optional infeasible orders; a missing grid on a nongrid fit is a
missing-extra request, rather than corrupt current state. Do not recompute
weights or fit results to validate retained metadata. Once full state is
validated, reconstruct weight matrices only for the selected target columns,
retaining original target IDs, and allocate known output shapes once.

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

## Numerical diagnostic recovery

Preserve ordinary representable arithmetic when correcting an exceptional
range failure. Within-arm probability-weight ESS is unchanged by a positive
common scaling. Scale before squaring when the squared arm mass is subnormal
or zero despite positive weight; retain unavailable ESS for genuine zero mass.
Use one owner for both treatment/instrument arms and all fresh/derived consumers.
This diagnostic recovery does not change weights, estimates or covariance.

CLATE's near-zero denominator flag retains its existing relative criterion and
all-original-signal scale, including zero-weight observations. A finite negative
signal can have an unrepresentable positive absolute magnitude in Mata. Use a
safe common binary scale in that exceptional case while preserving ordinary
cutoff decisions, zero-denominator behavior and the meaning of tolerance zero.
The flag remains descriptive. Mata's
[numeric-limit documentation](https://www.stata.com/manuals/m-5mindouble.pdf)
establishes the asymmetric nonmissing range; the shared recovery convention is
the packages' numerical policy, not a new statistical guarantee.

Detect unsafe magnitudes from the original signed values. Mata's
[minimum/maximum functions](https://www.stata.com/manuals/m-5minmax.pdf) ignore
missing elements; testing only the resulting maximum would miss an overflowed
absolute value that was dropped from the reduction.

Tests need independent concentration/relative-magnitude references, normal
controls, positive subnormal and zero-square cases, both arms and common
rescaling. Follow the corrected facts through saved prediction, candidates and
cold restoration. Keep exact discrete statuses separate from continuous numeric
tolerances and preserve the reported cause of genuine numerical failures.

## Unit and regression tests

Current schema 4 retains `e(target_notices)` (distance-tie and weak-localization
columns) and the aggregate `e(n_tied_targets)` / `e(n_weak_targets)` counters.
Tie checks reuse the rank/group validation pass; weak flags use existing
diagnostic concentration and its unchanged `.9*n`, `n>=2` rule. Shared
fresh/derived posting sets facts before reporting. Selected results and saved
prediction/candidate operations return their batch's counters, preserving the
active fit for nonposting calls. Standalone weights return analogous `r()`
counts. These are diagnostic facts, not warning-history flags or a last-error
store. All-target counts operate on the matrix directly rather than expanding
a potentially oversized numlist. Protect notice rendering inside caller-state
cleanup as well as computation; unexpected notice errors and Break must restore
`e()` before propagating.

An omitted target selector means the full numeric target index at shared reuse
owners. Do not expand it through `numlist`: Stata's 2,500-element numlist limit
is smaller than SE/MP result-matrix capacity. Explicit `targets()` retains the
host's documented syntax and ordering rules. Test the shared numeric selector
beyond that limit and exercise full public reporting where the edition's
matrix capacity permits it; large-panel graph rendering is a separate concern.

New behavior needs valid, invalid and boundary cases in its owning suites.
Each real regression should fail for the original defect at the most stable
public boundary available. Reuse established fixtures, independent references
and caller-state comparators; do not add copied corruption runners or one-off
frameworks. Register suites once and update family filters and test maps.

Name suites and section labels for the current observable behavior rather than
an audit, merge milestone or superseded rule. The test map covers every runtime
suite explicitly or through a meaningful filename pattern. Add a short purpose
comment where numerical references or neighboring owners need explanation;
descriptive assertion/section labels carry the individual scenario contract.
Keep independent kernel references separate from public orchestration checks.

Split unrelated responsibilities when that improves ownership and focused
execution; a large cohesive suite does not need splitting merely for its size.
Moved suites create their own deterministic data and load their own helpers,
restore instrumented providers/selectors and work through an exact master-runner
filter. Preserve assertion bodies, independent oracles, tolerances and actual
engine requirements; verify independent runs and a combined ordered run. Update
registration, current crosswalks and execution examples together. Preserve
historical release records and frozen fixture provenance. This applies the
companion's [test-design policy](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/09_testing_and_regression_policy.md)
to Stata's explicit runner and session-state model.

Source instrumentation uses unique temporary files/directories with cleanup
that also runs after an unexpected fixture error. Keep dated session receipt
paths out of executable tests and restore provider definitions and search paths
after replacing them. Structural efficiency checks count actual validated
snapshots and columns supplied to real weight primitives; timings remain
supporting review evidence.

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

## Removing unused implementations and tests

Check production reachability before pruning. Ado autoload, explicitly loaded
multi-program modules, typed Mata calls and function-pointer lookups are roots
that a plain name count can miss. In particular, `rdistnn_loo_prefix()` remains
active through `findexternal("rdistnn_loo_prefix()")`; its accumulation helpers
remain necessary. The official [findexternal manual](https://www.stata.com/manuals/m-5findexternal.pdf)
documents the function-name suffix used by that lookup.

Remove an unused private adapter together with its reload guard and tests of
its obsolete shape. Preserve substantive cases by calling the current shared
localization, public command or native transport. Numerical oracles stay
independent, and invalid current saved-state tests remain necessary. This
applies the companion's [test-maintenance policy](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/09_testing_and_regression_policy.md)
without imposing R's loading conventions on Stata.

Apply the same rule to repository/tooling tests: protect the current runner,
archive and validation behavior rather than maintaining a catalogue of removed
symbols or workflow filenames. A deliberate invalid-schema/policy fixture still
tests the current state boundary; it does not provide historical compatibility.

The current plugin exposes packed deletion transport only. Its numerical
kernels remain directly testable; automatic Mata execution when optional native
capability is absent remains supported current behavior. It is separate from
compatibility with earlier plugins or package objects. Instrument performance
at the actual operation and prove counters with positive controls before
accepting zero work. Retain timings outside CI pass/fail thresholds.

## Developer tooling contracts

Use private current-source copies for Clang, Bear and sanitizer builds. Preserve
the release Makefile, installed plugin, generated binaries and global compiler
configuration. Bound compilation to two jobs and keep runtime thread controls
separate. Standalone C++ tests establish native kernel/protocol behavior;
licensed Stata tests establish ado/Mata dispatch and host integration.

SDK fixtures must reproduce Stata's exact missing-value boundary, including
finite neighbors on both sides of the signed range. Calibrate the callback
model against real `maxdouble()`/missing-value behavior; an abbreviated decimal
cutoff can incorrectly reject valid inputs. The
[numeric-limit manual](https://www.stata.com/manuals/m-5mindouble.pdf)
owns those host limits. Keep public input validity distinct from a subsequent
computed-distance range failure.

Test live interruption separately from propagation of an injected return code
1. Long native work needs safe cooperative cancellation and cleanup using the
host's supported mechanism; a successful short SDK fixture does not establish
responsiveness to Break in the running Stata process. Keep host API calls out
of worker threads and validate cancellation on the actual supported host.

Capture compile databases from actual uncached compilation, including the SDK,
and restrict findings to handwritten translation units and headers. Retain the
captured source and reject stale native/config/compiler identities before tidy.
Cache compiler objects with strict identity and separate profile keys, never
test results. Keep cache/workspace writes inside safe private paths and reject
symbolic-link traversal. Source snapshots include current uncommitted inputs
while excluding private and generated output.

A successful capture binds the retained project native/configuration bytes
and compilation database to its receipt. Check agreement around source copying
and verify retained identity before accepting analysis; a failed, incomplete or
changed capture requires a fresh `compdb`. Never repair the receipt silently.
Clang's [database format](https://clang.llvm.org/docs/JSONCompilationDatabase.html)
supplies commands/paths/flags; these content checks are the package's additional
current-source provenance contract.

ASan/UBSan flags belong at both compile and link steps, with symbolization and
nonzero failures preserved. A negative detector fixture proves instrumentation
is active. Standalone leak checking requires process inspection; do not hide a
sandbox limitation behind a successful receipt. These profiles supplement
parity tests and supported-platform release gates.

Run Python correctness rules on maintained tools/tests with a pinned Ruff
version. Formatting requires explicitly selected files; no automatic mass
conversion or baseline suppression. Exclude generated/vendor/output paths even
when supplied explicitly. PowerShell analysis uses separate Windows PowerShell
5.1 launcher/common and PowerShell 7 plugin-builder compatibility targets.
Check syntax, commands and types against bundled profiles and fail all selected
warning/error diagnostics. Conservative profiles are static evidence, not a
replacement for real Windows execution. Analyzer helpers never launch Stata.

Requested checks fail actionably for missing tools, zero eligible files and
real diagnostics. Infrastructure tests use temporary positive/negative controls
for tool failures, selection/exclusion and source/output safety; avoid tests
that merely repeat command construction. Keep fast analyzers in the existing
static CI job and heavy profiles explicit. Optional developer tools remain
outside installed package dependencies. [CONTRIBUTING.md](../CONTRIBUTING.md#developer-analysis-tools)
owns installation and commands.

Reference basis: [clang-tidy](https://clang.llvm.org/extra/clang-tidy/),
[AddressSanitizer](https://clang.llvm.org/docs/AddressSanitizer.html),
[UndefinedBehaviorSanitizer](https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html),
[ccache 4.9](https://ccache.dev/manual/4.9.html),
[Ruff configuration](https://docs.astral.sh/ruff/configuration/) and
[PSScriptAnalyzer](https://learn.microsoft.com/en-us/powershell/utility-modules/psscriptanalyzer/overview).

## Reference basis

See the official [Stata error codes](https://www.stata.com/manuals/perror.pdf),
[display conventions](https://www.stata.com/manuals/pdisplay.pdf),
[capture behavior](https://www.stata.com/manuals/pcapture.pdf),
[Mata errors](https://www.stata.com/manuals/m-5error.pdf) and
[literal Mata display](https://www.stata.com/manuals/m-5display.pdf).
The official [size limits](https://www.stata.com/products/detailed-size-limits/)
distinguish numlist capacity from matrix dimensions.
Stata's 498 is reserved for situation-specific explanations, so a raw count of
that code does not establish misuse. Use the specific cause and documented
recovery boundary when selecting a code.
The companion design uses [base R conditions](https://stat.ethz.ch/R-manual/R-patched/library/base/html/conditions.html)
and [R warning protocols](https://stat.ethz.ch/R-manual/R-patched/library/base/html/warning.html),
[Advanced R conditions](https://adv-r.hadley.nz/conditions.html#custom-conditions)
and [R Packages' condition tests](https://r-pkgs.org/testing-basics.html#testing-errors).
