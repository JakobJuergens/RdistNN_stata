# Releasing RdistNN for Stata

The initial public release focuses on CATE and nonparametric regression,
with DNN weights, tuning, diagnostics, and plots. Experimental PLQR evaluates
localized score moments at supplied coefficients; publication does not turn
it into a coefficient estimator or add coefficient inference.

## Establish the release snapshot

Use a clean repository checkout with the tools in
[CONTRIBUTING.md](CONTRIBUTING.md). Review the version/tag and release notes.
The command-specific `*! version` headers have their own revision histories;
do not rewrite all of them merely to match a package release tag.

1. Run static package integrity and maintenance checks:

   ```bash
   bash tools/ci/pkg_check.sh
   python3 -m unittest discover -s tests/ci -v
   git diff --check
   ```

2. Run the full licensed Stata suite, onboarding example, and isolated install
   check on a Windows machine:

   ```powershell
   powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_runtime_check.ps1
   powershell -NoProfile -ExecutionPolicy Bypass -File tools/ci/stata_install_check.ps1
   ```

   In Stata, after `do dev-setup.do`, execute
   `do examples/getting-started.do`. Use `graph display cate_example` and
   `graph display np_example` to view its stored graphs. Inspect plots and stored
   results. The install helper uses temporary `PLUS` and `PERSONAL`
   directories and verifies installed command resolution and plugin parity,
   rather than relying on the development adopath. Record runtime results
   separately from static checks; static checks do not establish Stata execution.

3. Run the host plugin build/core check:

   ```powershell
   pwsh -NoProfile -File tools/ci/plugin_build_check.ps1
   ```

   Require the existing hosted Windows x86-64, Linux x86-64, macOS arm64,
   and macOS x86-64 build matrix on the release commit. A host check covers
   only that host. Review failures and verify the bundled Windows plugin
   separately through the Stata parity/installation tests.

4. Build the source ZIP with the same helper used by CI:

   ```bash
   bash tools/ci/build_source_archive.sh _dist/source.zip
   sha256sum _dist/source.zip > _dist/SHA256SUMS-source.txt
   ```

   The helper includes package metadata, README, release notes, public command
   guides/examples, and canonical ado/Mata/help/native sources. It validates
   the archive with `tools/ci/archive_check.sh`. Ignored build products,
   private notes, tests, CI configuration, and logs are excluded. Confirm that
   all package-manifest files and guide/example links resolve in the extracted
   ZIP. `net install rdistnn, from("path/to/extracted/archive") replace`
   installs its command payload; developer checks require a repository checkout.

Establish the final content before expensive checks. Staging, committing,
and pushing do not invalidate a verified content snapshot. Review the final
diff and hosted checks before tagging. Archives and logs stay under `_dist/`
or temporary directories and are not committed.

## Publish the release

Create the intended tag and GitHub release from the verified commit. Adapt
the public-release entry in `RELEASE_NOTES.md` into the announcement, with
explicit PLQR status and a link to the getting-started example.

The existing package-snapshot workflow attaches the validated source ZIP as
`rdistnn-source.zip`. The native-plugin workflow attaches four named plugin
assets and `SHA256SUMS.txt`. Confirm both workflows finish successfully and
the assets are available before announcing them. Source ZIP and plugin checksums
serve different payloads; preserve both. Retain the default-branch installation
URL in the README and provide the tagged source ZIP for reproducible installation.

GitHub's automatically generated repository archive is a checkout snapshot;
`rdistnn-source.zip` is the curated installable distribution. The ZIP bundles
a Windows plugin; Linux/macOS users install the matching release plugin or
use the portable Mata path. Do not substitute a Linux/macOS binary for the
ordinary `net install` Windows asset.

## Website introduction

[WEBSITE.md](WEBSITE.md) contains software-entry and announcement copy for
jakobjuergens.com. Link to the public repository and release after the release
assets are live. This preparation does not create a tag, publish a release,
or edit the live website.
