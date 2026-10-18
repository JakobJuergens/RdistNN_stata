# Repository Instructions for Codex

This file applies to the `RdistNN_stata` repository. Read the enclosing
workspace's `AGENTS.md` when present: `/home/jakob/research/AGENTS.md` in WSL,
or the applicable `U_Statistics/AGENTS.md` in a Windows checkout.

## Development and production branches

Develop on the local `development` branch (renamed from `feature/clate`). Keep
it unpublished until publication is authorized; branches in a public repository
are public when pushed. Production is `master`. Submit production changes
through a pull request and squash-merge after review and required checks.
Do not push development commits directly to `master` without an explicit
maintainer instruction. See `CONTRIBUTING.md` for the workflow.

## Project Shape

- `RdistNN_stata` is the Stata companion package for the `RdistNN` R package
  and the `Cond_Z_Estimation` paper.
- Estimation/weight commands are `np_dnn`, `cate_dnn`, experimental
  `plqr_dnn`, experimental `clate_dnn`, and `weights_dnn`; NP, CATE and CLATE
  also have diagnostics/plot commands.
- Read `README.md` before changing command behavior, package layout, or
  installation instructions. User documentation starts at `docs/index.md`
  and installed `help rdistnn`; see CONTRIBUTING for documentation checks.
- In WSL, use `../RdistNN` and `../Cond_Z_Estimation` for theory/software
  alignment. In the Windows `U_Statistics/02_Code/RdistNN_stata` layout,
  the corresponding paths are `..\RdistNN` and `..\..\Cond_Z_Estimation`.
  Resolve the current layout before using either set.

## Local CLI Tooling

Resolve tools in the environment actually running the command. In WSL, use
`command -v TOOL`; in Windows PowerShell, use `Get-Command TOOL`. A tool
installed on one machine is not evidence that it is on another machine's PATH.

- Use `rg` for source search; `jq` or Python for structured metadata.
- Use the package helpers under `tools/ci/` for static, runtime, installation,
  and native-plugin checks; their platform requirements differ.
- Use `Rscript` for R work, `actionlint` for workflows, and `wolframscript`
  when a mathematical check calls for it. Follow the task-specific checks below.
- Poppler (`pdfinfo`, `pdftotext`, `pdftoppm`) and `qpdf` support document
  inspection. Windows Ghostscript is commonly `gswin64c`; use `magick` for
  ImageMagick. Check availability before promising a build or check.
- For Literature documents, use the shared catalog/extraction workflow below;
  do not replace its preserved evidence with ad hoc or in-place OCR.

## Working Rules

- Keep ado/Mata/plugin behavior aligned with the documented command contracts
  in `README.md`.
- Use `bash tools/ci/pkg_check.sh` for static package integrity checks.
- Generate source archives under a temporary or ignored build directory, then
  run `bash tools/ci/archive_check.sh <path-to-source.zip>`. Validated release
  archives are published as workflow artifacts and GitHub release assets, not
  committed under `archives/`.
- GitHub Actions failure logs are finite-retention workflow artifacts. Do not
  commit generated `GH_logs/` content or downloaded debug copies.
- Use `powershell -NoProfile -ExecutionPolicy Bypass -File
  tools/ci/stata_runtime_check.ps1` for the full local Stata runtime suite.
- Use `powershell -NoProfile -ExecutionPolicy Bypass -File
  tools/ci/stata_install_check.ps1` to verify a clean `net install` into
  temporary `PLUS`/`PERSONAL` directories.
- For `weights_dnn` plugin parity work, use `global rdistnn_use_plugin 0` to
  force the Mata reference path and `global rdistnn_use_plugin 1` to require
  the L2 plugin path.  The default is to use the plugin when installed.
- Do not treat generated logs such as `tests.log` as source of truth unless the
  task is explicitly about test-output archaeology.
- If Stata itself is needed for verification, first check whether a CLI command
  is available locally; otherwise record that Stata is not available on `PATH`.
  On the Windows development installation, the expected StataNow 19 executable is
  `C:\Program Files\StataNow19\StataSE-64.exe`.

## Literature and source lookup

For local literature work, use `literature-vault`; use `obsidian-paper-notes`
for Obsidian paper/book notes. If a skill is not exposed, read its maintained
`SKILL.md` under the vault's `_ai/skills/`. Read the vault's `AGENTS.md` and
`_ai/AGENTS.md` before source or citation work.

The WSL vault is `/home/jakob/research/Literature`; Windows uses
`D:\Jakob_Clouds\Dropbox\Literature` on the desktop and
`C:\Users\jakob\Dropbox\Literature` on the laptop. Preserve both Windows
layouts and resolve the current machine's paths. Project-local `01_Literature`
collections and extraction caches are separate from this full vault.

Use the vault's `_bibliography/scripts/registry.py` by absolute path: `search`
for passages, `search --find` for discovery, and `lookup` for citation identity.
Preserve existing project keys and unresolved identities. Check `status`'s
`catalog_schema_version` and `variants DOCUMENT_SHA`: baseline-only catalogs
use original sidecars plus `supplemental` MinerU reads; registered v5 catalogs
support preferred `read` and `--variant all --compare`. Follow the vault's
`_extractions/README.md` and `_bibliography/QUALITY.md` for exact commands and
review scope. Preferred text is not verification; inspect relied-on mathematics,
tables, figures, and uncertain exact wording unless a current scoped visual
review covers it. Keep extraction identities and physical/printed pages distinct.
The personal vault is a research resource, not a repository build dependency.
