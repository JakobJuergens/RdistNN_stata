# Stata companion public-release preparation

Status: prepared and verified locally, 2026-10-14. Publication is a separate step.

The user's objective is preparation analogous to the R package's release
polish and test-maintenance passes. This is an existing-package documentation,
packaging, and verification pass on `master`, preserving estimator behavior.

## Scope and plan

- [x] Refocus README and package metadata on CATE and nonparametric regression;
  preserve the pointwise inference, nuisance, tuning, and plugin contracts.
  Route advanced command details to a public guide and verify deterministic
  onboarding examples. Use the author's public research page for the paper.
- [x] Add contributor instructions, a release checklist, and website copy;
  retain explicit experimental PLQR supplied-coefficient moment semantics.
- [x] Keep the already-focused regression files and their numerical assertions.
  Add a test map, filtered master runs, reset counters for each master run,
  and a check that every regression file is registered exactly once.
- [x] Correct the verification helpers' Windows UNC path handling and quote
  do-file arguments. The baseline runtime attempt failed because Stata received
  `Microsoft.PowerShell.Core\\FileSystem::` rather than a native filename.
- [x] Provide one source-archive builder shared by local release preparation
  and CI, including public guides/examples and excluding ignored build output.
- [x] Run static checks, the complete Windows Stata suite, focused runner and
  onboarding checks, an isolated installation, native plugin/core checks,
  archive-content checks, and an independent final review.

## Guardrails

Keep estimation/Mata/native source and the bundled plugin byte-identical. Preserve
all fixed R fixtures, numerical tolerances, command defaults, aliases,
returned results, and optional plugin controls. Do not introduce Stata
features merely because the broader R package has them. Keep historical
review notes intact. Local archives/logs stay ignored. Preparing this work
does not create a release, deploy a website, or publish repository changes.

## Baseline

The checkout was clean on `master`. Copies and hashes of all 86 tracked files
are in `/tmp/rdistnn-stata-release-baseline/`. Static package checks and four
Python repository-policy tests passed. Windows StataNow was found locally;
the first runtime attempt exposed the UNC path bug before executing tests.

The literal README check additionally exposed CATE-then-NP plot autoloading:
NP attempted to reload CATE's public ado to access locally scoped routines.
A focused regression fails against the original code. Shared plot-builder
bodies move byte for byte into `src/_rdistnn_plot_helpers.ado`; only loading
changes in the two public plot commands. This exception is necessary for
working onboarding and is separate from estimation or numerical changes.

## Final verification

- Licensed Windows StataNow 19.5, invoked from this WSL checkout: **683
  assertions passed** in the complete suite. The original 678 assertions
  remain; five plotting-order assertions were added. The focused regression
  failed two assertions before the shared-helper move and passes afterward.
- Every original regression file and all fixed R golden values remain
  byte-identical. Estimation commands, diagnostics commands, Mata, native
  plugin source, and the bundled Windows binary remain byte-identical.
  All shared plot-builder bodies also match the baseline byte for byte.
- Literal README/command-guide examples execute successfully. Oracle inputs,
  pointwise uncertainty, and the narrower experimental PLQR contract remain
  explicit. Repeated graph examples use `replace`; the corrected overlap
  example follows a nuisance-based CATE fit. Updated PLQR SMCL help renders
  successfully with Stata's translator.
- Fresh focused runs and repeated-run counter resets pass. An unknown pattern
  returns error 198. Native UNC paths and a do-file path containing spaces
  execute successfully. Both source-checkout installation and installation
  from the **extracted final source ZIP** pass, including onboarding before
  other plots and installed shared-helper resolution.
- Static package integrity and six Python repository/maintenance tests pass.
  The maintenance tests are wired into the existing package-check workflow.
  Workflow lint, Markdown path/anchor checks, and whitespace checks pass.
- The Linux x86-64 native plugin builds and its deterministic core tests pass.
  Windows plugin parser and forced plugin/Mata parity run in the Stata suite
  and isolated-install check. Other native platforms remain hosted release
  gates; this local preparation does not claim a current four-platform result.
- `_dist/source.zip` contains 33 public source/documentation files, including
  the canonical plugin sources, Windows install binary, public guide, and
  example. Every ZIP entry matches its current source bytes, every local
  Markdown link resolves inside the ZIP, and required/forbidden/duplicate
  checks pass. SHA256:
  `06b4326f644874843745d1f9be406f025a52002ba0e1f7643479479b02d75474`.
- Independent review confirmed numerical-body preservation and found no further
  helper-loading issue. Its two
  maintenance findings were fixed: Python guards now run in CI, and the
  command guide's stale inference reference routes to README. The R/Stata
  crosswalk also now accurately describes R's supplied-vector PLQR inference
  without attributing local root solving to either companion.

Receipts: `/tmp/rdistnn-stata-final-runtime.log`,
`/tmp/rdistnn-stata-final-install.log`,
`/tmp/rdistnn-stata-archive-install.log`,
`/tmp/rdistnn-stata-final-examples.log`,
`/tmp/rdistnn-stata-final-help.log`,
`/tmp/rdistnn-stata-final-static.log`,
`/tmp/rdistnn-stata-final-python.log`,
`/tmp/rdistnn-stata-final-native.log`, and
`/tmp/rdistnn-stata-final-archive.log`.
The source ZIP, checksum, rendered help, and native build output remain ignored
under `_dist/`. No tag, release, website deployment, commit, or push was made
as part of this preparation pass.
