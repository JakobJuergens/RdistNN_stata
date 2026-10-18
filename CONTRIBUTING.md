# Contributing to RdistNN for Stata

The initial public release focuses on CATE, nonparametric regression, weights,
diagnostics, and plots. PLQR moment evaluation is experimental. Read
[AGENTS.md](AGENTS.md) before editing, and use installed command help and
[the article collection](docs/index.md) and [command guide](docs/command-guide.md) for current contracts. Open an issue
before proposing a new estimator or a public-interface change.

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

The isolated-install helper also executes all six public examples using only
the temporary PLUS command payload and verifies the installed overview.
Both helpers accept `-StataExe` and `-RepoRoot`. The expected local Windows
installation is `C:\Program Files\StataNow19\StataSE-64.exe`; discovery also
supports other Stata editions and `STATA_EXE`. Resolve the executable on the
machine actually running the check. PowerShell uses native filesystem paths
and quotes do-file names, including paths containing spaces.
The Stata setup and test runner normalize Windows backslashes for ado discovery
from UNC shares, including a WSL checkout.

Force `global rdistnn_use_plugin 0` for Mata reference checks or `1` to require
the L2 plugin. Clear the global afterward with `macro drop rdistnn_use_plugin`.
The native build helper executes deterministic core tests and inspects the
host binary:

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
checkout, execute all six seeded examples and render the package overview with:

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
feature-branch URLs rather than nonexistent archive paths. Preserve
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
