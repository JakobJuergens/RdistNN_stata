# Releasing RdistNN for Stata

The package provides CATE and nonparametric regression, with DNN weights,
tuning, diagnostics, and plots. Supported CLATE adds binary-IV conditional
effects with paired frozen-score covariance. Experimental PLQR evaluates
localized score moments at supplied coefficients; publication does not turn
it into a coefficient estimator or add coefficient inference.

## Establish the release snapshot

Use a clean repository checkout with the tools in
[CONTRIBUTING.md](CONTRIBUTING.md). Review the version/tag and release notes.
The next planned release is `0.1.2`, currently an unpublished candidate on
`development`. It adds duplicate-covariate
block averaging, deletion-aware survivor weights and saved-model schema 2,
with a theory warning; see [RELEASE_NOTES.md](RELEASE_NOTES.md). The next
production PR will carry these changes under `0.1.2`; assigning that version
locally does not publish a release.

Review the known [asymptotic-tail limitation](docs/numerical-stability.md#asymptotic-tail-limits)
when assessing readiness. Matching native/Mata deletions for fixed tables does
not establish correct construction or cross-language support parity at
subnormal scales.

The `0.1.1` baseline combines numerical fixes, covariate-based exact-distance tie
ordering, shared modules and a current-only interface, aligned with the R companion.
The previous version `0.1.0.9004` entered production `master` through the reviewed squash
merge of [PR #10](https://github.com/JakobJuergens/RdistNN_stata/pull/10) at
[`3e02115124e2a8f71e73739aa95fa1e34a754f13`](https://github.com/JakobJuergens/RdistNN_stata/commit/3e02115124e2a8f71e73739aa95fa1e34a754f13),
aligned with the R companion. Version `0.1.0.9003` identifies the previous
production snapshot. Exact
tables, compensated reductions, whitening rank policy, shared covariance and
the native Beta nonconvergence guard first shipped in the merged
0.1.0.9003 snapshot. Earlier TDNN/distance work remains under 0.1.0.9002. Retain the
intended release version for candidates; create tagged versions only when
authorized.
Updating `master` does not itself publish a tag or replace assets from the
existing `v0.1.0` release.
Keep `rdistnn.pkg`'s `Version` and `Distribution-Date` metadata aligned with
the release notes and actual publication date. Its `v 3` line identifies the
Stata package-file format.
The command-specific `*! version` headers have their own revision histories;
do not rewrite all of them merely to match a package release tag.

Keep the R and Stata package version families roughly aligned with shared
functional milestones. Stata may trail the R reference as corresponding
features are implemented and verified; patch and development revisions can
advance independently for implementation-specific fixes. A matching version
does not promise identical feature coverage. Review the
[R-to-Stata crosswalk](docs/r-documentation-crosswalk.md) when assigning a
version and retain its explicit R-only and experimental boundaries.

The next planned version-pinned baseline is `v0.1.2` in both packages. After
the production PRs are merged and the required checks pass on their final
commits, create each repository's `v0.1.2` tag and publish its verified assets.
Keep the `0.1.2` changelog entry marked unreleased until publication, then add
the publication date and update `Distribution-Date`. For installation after
that release is published:

```stata
net install rdistnn, from("https://raw.githubusercontent.com/JakobJuergens/RdistNN_stata/v0.1.2/") replace
```

This planned tag is not an installation source yet. Keep the ordinary README
installation on production `master` and retain historical tags, assets and
scientific dependency pins. Record the paired final commits and asset
checksums in the release-set manifest below.

Commit candidates locally. Publishing review branches/PRs, tags, assets or
websites requires separate authorization; do not push to `master`. Install the
previous immutable `0.1.0.9004` source as a historical reference with:

```stata
net install rdistnn, from("https://raw.githubusercontent.com/JakobJuergens/RdistNN_stata/3e02115124e2a8f71e73739aa95fa1e34a754f13/") replace
```

This pins the reviewed source snapshot. Unpinned `master` follows production
source; `v0.1.0` selects the existing older public tag and release assets.
Different source snapshots can share a development version, so retain the
full commit or checked archive identity. Unpublished future candidates require
a local checkout at their full commit or the checked curated ZIP; their GitHub
URLs do not exist yet.

After both candidate commits are final, bind them in one external release-set
manifest under ignored outputs. Include tested trees, artifact filenames/sizes/SHA-256, runtime and
compiler information, commands, hosted run IDs/tested SHAs and explicit
PASS/FAIL/BLOCKED/NOT_RUN gates. Keep historical scientific dependency identity
separate, and avoid circular committed HEAD references.

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

   Run the entire registered runtime suite without a test-name filter on the
   final source tree. Isolate stale logs first, require process exit zero and
   the framework's actual zero-failure summary, and bind the fresh log hash and
   invocation timestamps to that tree. A failed full run plus focused passes
   does not meet this gate. Do not assume a fixed assertion count. If licensed
   Stata cannot run, mark the gate BLOCKED and preserve that distinction from
   native builds and static checks. Also run the licensed documentation helper
   and inspect the curated ZIP through the isolated-install helper using
   temporary PLUS/PERSONAL directories and the archive payload as its source.

   In Stata, after `do dev-setup.do`, execute
   `do examples/getting-started.do`. Use `graph display cate_example` and
   `graph display np_example` to view its stored graphs. Inspect plots and stored
   results. The install helper uses temporary `PLUS` and `PERSONAL`
   directories and verifies installed command resolution and plugin parity,
   rather than relying on the development adopath. Record runtime results
   separately from static checks; static checks do not establish Stata execution.

3. Run the host plugin build, core/distance and protocol checks:

   ```powershell
   pwsh -NoProfile -File tools/ci/plugin_build_check.ps1
   ```

   Require the existing hosted Windows x86-64, Linux x86-64, macOS arm64,
   and macOS x86-64 build matrix on the release commit. A host check covers
   only that host. Review failures and verify the bundled Windows plugin
   separately through the Stata parity/installation tests.

   Before merging a feature branch into production, run the hosted checks on
   its final commit using manual dispatch (`gh workflow run WORKFLOW --ref
   BRANCH`). Run `pkg-check.yml`, `plugin-build.yml`, `package-snapshot.yml`,
   and `repository-policy.yml`, and require their successful completion.

4. Build the source ZIP with the same helper used by CI:

   ```bash
   bash tools/ci/build_release_assets.sh _dist
   (cd _dist && sha256sum --check SHA256SUMS-source.txt)
   ```

   The helper includes package metadata, README, release notes, public command
   guides/examples, and canonical ado/Mata/help/native sources. It validates
   the archive with `tools/ci/archive_check.sh`. Ignored build products,
   private notes, tests, CI configuration, and logs are excluded. Confirm that
   all package-manifest files and guide/example links resolve in the extracted
   ZIP, including `LICENSE.md` and `LINKING_EXCEPTION.md`.
   Keep the package-specific installed copies `src/rdistnn_license.md` and
   `src/rdistnn_linking_exception.md` synchronized with these notices, adjusting
   their relative links to the installed filenames. The maintenance tests check
   that they match. Uppercase `F` manifest entries make Stata install these
   notices into PLUS rather than treating them as optional ancillary files.
   `net install rdistnn, from("path/to/extracted/archive") replace`
   installs its command payload; developer checks require a repository checkout.

Establish the final content before expensive checks. Staging, committing,
and pushing do not invalidate a verified content snapshot. Review the final
diff and hosted checks before tagging. Archives and logs stay under `_dist/`
or temporary directories and are not committed.

## Documentation snapshot

The public article collection is indexed by [docs/index.md](docs/index.md)
and the installed overview is `help rdistnn`. Before release, run
`python3 tools/ci/documentation_check.py` and the licensed
`tools/ci/stata_documentation_check.ps1` helper, then verify the public docs,
help and examples from the curated ZIP. Command help must be in `rdistnn.pkg`.
Article/example links should resolve inside the ZIP; developer-only controls,
tests and tools use explicit repository URLs. An optional `--r-root` check
compares the crosswalk with the live R article/reference inventory.

The current articles describe the production interface in this development
snapshot. Keep ordinary installation links on `master` and use immutable tags
for pinned installs.
CLATE and the shared reporting additions belong to version 0.1.0.9000; the unchanged
v0.1.0 tag predates them. Preserve CLATE's documented pointwise inference boundaries
and PLQR's supplied-coefficient, moment-only scope in release copy.

## Publish the release

When a tagged release is requested, create the matching version tag and
GitHub release from the verified commit. Adapt its entry in
`RELEASE_NOTES.md` into the announcement, with explicit CLATE/PLQR status
and a link to the getting-started example. Keep existing version tags and
their assets unchanged.

The existing package-snapshot workflow attaches the validated source ZIP as
`rdistnn-source.zip` and its `SHA256SUMS-source.txt`. The native-plugin workflow
attaches four named plugin assets and `SHA256SUMS.txt`. Confirm both workflows finish successfully and
the assets are available before announcing them. Source ZIP and plugin checksums
serve different payloads; preserve both. Retain the default-branch installation
URL in the README and provide the tagged source ZIP for reproducible installation.
Before announcing the release, download its assets without credentials, verify
both checksum files from the download directory, and run the isolated installation
check using the extracted source ZIP. Verify the tagged `net install` URL too.

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
