# Contributing to RdistNN for Stata

The package provides CATE, single/two-scale nonparametric regression, CLATE,
weights, diagnostics and plots. PLQR moment evaluation remains experimental. Read
[AGENTS.md](AGENTS.md) before editing, and use installed command help and
[the article collection](docs/index.md) and [command guide](docs/command-guide.md) for current contracts. Open an issue
before proposing a new estimator or a public-interface change.

The maintained [implementation guidance](markdown_notes/implementation_guidance.md)
records error/warning severity, actionable copy, grouping reuse, selective
recovery and unit/performance testing requirements. Read it before changing
those surfaces. It records intended design; code, help and passing regressions
establish which behavior has been implemented.

## Development and production workflow

Use local `development` for ongoing work and keep it unpublished until
publication is authorized. GitHub branches in a public repository are public
when pushed. Historical branch names remain part of the Git record.

For an authorized production update, publish a review branch and open a pull
request targeting `master`. Describe the final behavior and verification,
complete the required CI checks and review, and squash-merge the PR. This keeps
intermediate development commits out of the production history users follow.
Direct pushes to `master` require an explicit maintainer override. Keep local
development current with the reviewed production state before the next change.

The first development commit after production `0.1.x` establishes
`0.1.x.9000`. Each subsequent development commit increments the fourth component,
including small changes: `0.1.x.9001`, `0.1.x.9002`, and so on. Uncommitted drafts
do not advance this sequence. Read the current version from `rdistnn.pkg`'s
`Version` and the current `RELEASE_NOTES.md` heading; update both together.
By default, each production merge becomes `0.1.(x+1)`; the maintainer may choose another
version. After that release, append `.9000` to the selected production version
to start development. For example, the cycle after production `0.1.2` starts at
`0.1.2.9000`, then `0.1.2.9001`; its next default production version is `0.1.3`.

## Current development compatibility

Use the current documented interfaces and complete current fit schemas.
During this early development stage, remove compatibility layers for older
package versions: retired aliases and wrappers, saved-state upgraders and
migration shims. Keep one current implementation unless the maintainer
explicitly requests an exception. Older, incomplete or inconsistent fits
require an actionable refit error; do not silently reconstruct or repair them.
Current fits must still survive serialization and support documented reuse.

Document intentional interface changes and update their regression tests.
Preserve historical changelogs, released tags/assets and frozen scientific
evidence. Compatibility and deprecation commitments can be revisited when the
packages stabilize; prior production releases do not change the current policy.

## Local setup

Use a repository checkout, Stata 17 or later, and Python 3.10 or later for static
maintenance/documentation checks and source archives. From Stata, change to the repository
root before adding local commands:

```stata
cd "path/to/RdistNN_stata"
do "dev-setup.do"
do "examples/getting-started.do"
```

`dev-setup.do` resolves `src/` from the current working directory. It is a
session setup step; running it by an absolute filename from another directory
does not change the working directory automatically. Native plugin development
additionally requires a C++17 compiler; see [src/plugin/README.md](src/plugin/README.md).

## Numerical parity contract

Require agreement both between R and Stata for shared methods and between
compiled C++ and R/Mata implementations of each method. Compare identical
numeric inputs and row order, targets, fitted geometry, signals, weight
normalization and covariance definitions. Keep exact versus asymptotic weights
and frozen-rank versus geometry-refitting deletions distinct when they define
different computations. Export numeric fixture data across languages; using
the same seed in R and Stata does not establish identical random draws.

Require exact ranks ordered by represented distance then original row index,
target-specific exact-distance group membership and rank-weight averaging,
support, target identities,
statuses and missing-value masks after documented representation mappings.
Compare continuous estimates, distances, weights and covariance at explicit,
justified tolerances, using exact numeric equality where arithmetic is fully
specified. A small numeric difference cannot excuse a changed neighbor,
support or status. Investigate discrepancies before widening tolerances or
regenerating fixtures.

The [known asymptotic-tail limitation](docs/numerical-stability.md#asymptotic-tail-limits)
remains unresolved. Fixed-table engine agreement does not establish correct
weight construction, and bounded Beta fixtures do not certify every tail scale.
Document the discrepancy rather than broadening tolerances or treating changed
support as continuous rounding.

Keep decision-sensitive operation order and precision aligned. For distances,
subtract before applying frozen inverse scales, accumulate coordinates in
double precision and round L2 squared contributions before addition. Preserve
the shared support-aware centering and deletion/covariance contracts. Reuse
helpers across estimators. Compiler contraction, wider accumulation and
parallel reductions require regression coverage when they can affect decisions.
Align engines toward stable independent references: equivalent formulas using
`ln1p`/`expm1`, direct tail probabilities and scaled sums of squares can avoid
cancellation and intermediate range failures. Use compensated accumulation
selectively and benchmark constant runtime/allocation costs before replacing
vectorized operations. Apply a stability correction across affected engines;
agreement between two inaccurate implementations is insufficient.
Platform math libraries and linear algebra can still cause small numeric
differences; universal bit-for-bit identity is not a package-wide promise.

Force and establish each engine independently, and test automatic fallback
separately. Native acceptance requires the real plugin API; a portable-only
run must report native checks as unavailable. Retain independent analytic or
explicit-loop references. For an affected primitive, cover relevant ties,
offsets, tiny/constant coordinates, scale clamping, signed weights, zero
variance, missingness and overflow. Follow it through batching, target
reordering/duplicates, candidates, covariance and restored prediction where
applicable. Run relevant native regressions on supported build architectures
and record compiler/math-library/thread settings when investigating differences.
The [test-suite map](tests/README.md) and
[R test crosswalk](tests/R_TEST_CROSSWALK.md) identify existing coverage.

## Verification

Choose checks for the changed surface; use the full suite for shared Mata or
command infrastructure. Reuse successful results while source is unchanged.

| Surface | Check from the repository root |
| --- | --- |
| Package structure, public documentation and policy | Commands below |
| Estimator or numerical behavior | Focused suites from the [test map](tests/README.md), then full runtime for shared changes |
| Native code | [Plugin build/tests](src/plugin/README.md#build-from-source) and relevant licensed native/Mata regressions |
| Examples/help | Execute changed examples; the documentation helper runs all eight and renders the overview |
| Installation/release | Isolated-install check and [release gates](RELEASING.md) |

```sh
bash tools/ci/pkg_check.sh
python3 -m unittest discover -s tests/ci -v
python3 tools/ci/documentation_check.py
```

`pkg_check.sh` includes manifest/source hygiene: public and private ado/help
inputs, version declarations, the `rdistnn_` Mata prefix and unique typed
function ownership. With sibling R available, also run
`python3 tools/ci/package_hygiene.py --r-root ../RdistNN` for byte-identical
shared headers. Optional `--r-root ../RdistNN` on the documentation checker
compares the [documentation crosswalk](docs/r-documentation-crosswalk.md)
with R's live article/reference inventory.

With licensed Windows Stata:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_runtime_check.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_install_check.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_documentation_check.ps1
```

The helpers accept `-StataExe` and `-RepoRoot`; resolve native paths on the
machine running them. Discovery supports `STATA_EXE` and Stata editions,
including the local `C:\Program Files\StataNow19\StataSE-64.exe`. Paths with
spaces and Windows/UNC checkout paths are supported. Each helper requires
both exit code zero and a fresh completion marker. Runtime `-KeepLog` archives
only the preceding log as `tests.previous.log`; current execution must create
its own `tests.log`.

The isolated-install helper uses temporary PLUS/PERSONAL directories and runs
all eight examples using only the installed commands. The documentation helper
writes ignored `build/documentation/` output and `documentation_examples.log`.
Add `-UpdateFigures` to regenerate the four gallery PNGs in `docs/figures/`,
then inspect them. The default leaves committed figures unchanged.
After the full runtime suite, `python3 tools/ci/plot_svg_check.py` verifies four
exported filled bands and three profiles without uncertainty against numeric
CSV data; receipts stay under ignored `build/premerge/`.

The [test map](tests/README.md) owns suite selection, fixture generation and
engine controls; the [numerical reference](docs/numerical-stability.md) owns
module responsibilities and arithmetic. Keep fixed R values unless their
source contract changes. User guides should teach tasks, installed help owns
syntax/results, and the paper owns theory. Retain practical interpretation
limits beside the outputs. Articles/examples are source-archive inputs;
checkout-only tests/tools need production-branch URLs from those articles.

## Reporting and review

Report issues at [GitHub](https://github.com/JakobJuergens/RdistNN_stata/issues).
Include a reproducible do-file, Stata version/platform, expected and actual
results, and whether the discrepancy persists with the Mata path forced.
Describe changes and verification in pull requests. Update
[RELEASE_NOTES.md](RELEASE_NOTES.md) for user-visible changes. Generated logs,
archives, compiled development binaries, and workflow debug artifacts stay
untracked. Release preparation and publication follow
[RELEASING.md](RELEASING.md).
Earlier plans and audits are indexed in [historical records](markdown_notes/history.md).
