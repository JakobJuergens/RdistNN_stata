# Contributing to RdistNN for Stata

The package provides CATE, single/two-scale nonparametric regression, CLATE,
weights, diagnostics and plots. PLQR moment evaluation remains experimental. Read
[AGENTS.md](AGENTS.md) before editing, and use installed command help and
[the article collection](docs/index.md) and [command guide](docs/command-guide.md) for current contracts. Open an issue
before proposing a new estimator or a public-interface change.

## Development and production workflow

Use the local `development` branch, renamed from `feature/clate`, for ongoing
work. It has no upstream until publication is authorized. GitHub cannot make
an individual branch private in a public repository, so unpublished work stays
local. The already published historical `feature/clate` branch is separate
from this local branch and is not automatically renamed or deleted remotely.

For an authorized production update, publish a review branch and open a pull
request targeting `master`. Describe the final behavior and verification,
complete the required CI checks and review, and squash-merge the PR. This keeps
intermediate development commits out of the production history users follow.
Direct pushes to `master` require an explicit maintainer override. Keep local
development current with the reviewed production state before the next change.

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

Require exact ranks, stable original-row tie order, support, target identities,
statuses and missing-value masks after documented representation mappings.
Compare continuous estimates, distances, weights and covariance at explicit,
justified tolerances, using exact numeric equality where arithmetic is fully
specified. A small numeric difference cannot excuse a changed neighbor,
support or status. Investigate discrepancies before widening tolerances or
regenerating fixtures.

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

Static checks from the repository root:

```bash
bash tools/ci/pkg_check.sh
python3 -m unittest discover -s tests/ci -v
python3 tools/ci/documentation_check.py
git diff --check
```

The [test-suite map](tests/README.md) explains focused runs and fixture
ownership. On Windows with licensed Stata installed, run the full suite and
a clean installation:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_runtime_check.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_install_check.ps1
```

The isolated-install helper also executes all eight public examples using only
the temporary PLUS command payload and verifies the installed overview.
Both helpers accept `-StataExe` and `-RepoRoot`. The expected local Windows
installation is `C:\Program Files\StataNow19\StataSE-64.exe`; discovery also
supports other Stata editions and `STATA_EXE`. Resolve the executable on the
machine actually running the check. PowerShell uses native filesystem paths
and quotes do-file names, including paths containing spaces.
The Stata setup and test runner normalize Windows backslashes for ado discovery
from UNC shares, including a WSL checkout.
All three Stata verification helpers require a zero process exit code and a
fresh completion marker. Runtime `-KeepLog` copies the previous run to
`tests.previous.log`; the current run must create its own `tests.log`.
After the full licensed suite, verify the actual exported confidence-band
and `nose` geometry against its numerical CSV data with
`python3 tools/ci/plot_svg_check.py`. It checks four filled bands and three
profiles without uncertainty, and writes an ignored JSON receipt under
`build/premerge/`.

Use `global rdistnn_distance_engine mata` or `plugin` for explicit distance
parity checks. Drop it for auto dispatch, where `rdistnn_use_plugin 0` skips
native calls and `1` requires them. Reset test globals after use. The native
build helper executes core, distance and real-entry-point protocol tests and
inspects the host binary:

```powershell
pwsh -NoProfile -File tools/ci/plugin_build_check.ps1
```

Run numerical or validation tests affected by a behavior change, then the full
suite for changes to shared Mata or command infrastructure. Documentation-only
changes should execute changed examples and verify help and links. The
article and installed-help checks require no R installation; for the optional
live R article/reference inventory, pass `--r-root` to
`tools/ci/documentation_check.py`. Keep [the documentation crosswalk](docs/r-documentation-crosswalk.md)
in sync when either package's public topics change. From a licensed Windows
checkout, execute all eight seeded examples and render the package overview with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_documentation_check.ps1
```

Add `-UpdateFigures` to regenerate the four committed gallery PNGs in
`docs/figures/` from those seeded examples; inspect them before committing.
The default check leaves the published figures untouched.

The helper accepts `-StataExe`/`-RepoRoot`, stops on any failed example and
checks a dedicated success marker. It writes ignored output under
`build/documentation/` and `documentation_examples.log`. The ordinary Stata
process exit code alone is not proof that a do-file completed. The examples
can also be run individually from the repository root after `dev-setup.do`.
Documentation articles and examples are public source-archive inputs; developer
tests/tools remain checkout-only, so article links to those use matching
production-branch URLs rather than nonexistent archive paths. Preserve
fixed R golden values; regenerate them only when their source contract changes,
using `Rscript tools/generate_r_priority_fixtures.R` with the sibling R checkout.
Reuse successful verification while the content snapshot is unchanged.

## Reporting and review

Report issues at [GitHub](https://github.com/JakobJuergens/RdistNN_stata/issues).
Include a reproducible do-file, Stata version/platform, expected and actual
results, and whether the discrepancy persists with the Mata path forced.
Describe changes and verification in pull requests. Update
[RELEASE_NOTES.md](RELEASE_NOTES.md) for user-visible changes. Generated logs,
archives, compiled development binaries, and workflow debug artifacts stay
untracked. Release preparation and publication follow
[RELEASING.md](RELEASING.md).
