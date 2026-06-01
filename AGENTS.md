# Repository Instructions for Codex

This file applies to the `RdistNN_stata` repository. The parent
`U_Statistics/AGENTS.md` also applies when working from the multi-repo
workspace.

## Project Shape

- `RdistNN_stata` is the Stata companion package for the `RdistNN` R package
  and the `Cond_Z_Estimation` paper.
- User-facing commands include `np_dnn`, `cate_dnn`, and `weights_dnn`.
- Read `README.md` before changing command behavior, package layout, or
  installation instructions.
- Prefer keeping theory-facing behavior aligned with `..\RdistNN` and
  `..\..\Cond_Z_Estimation`.

## Local CLI Tooling

The parent Windows workspace has these tools on `PATH` as of 2026-05-06:

- Search and discovery: `rg`, `rga` for PDFs/Office files/archives, and `fd`.
- Structured data: `jq`, `yq`, and `duckdb` for JSON, YAML, CSV, Parquet, and
  tabular scratch queries.
- Code and CI quality: `actionlint` for GitHub workflows, `shellcheck` and
  `shfmt` for shell scripts and workflow `run:` blocks, `pre-commit` for
  multi-tool hooks, and `ruff` for Python helper scripts.
- PDF/OCR/image work: `pdfinfo`, `pdftotext`, `pdftoppm`, `pdffonts`, `qpdf`,
  Ghostscript as `gswin64c`, ImageMagick as `magick`, Tesseract as
  `tesseract`, `ocrmypdf`, and ExifTool as `exiftool`.
- Package/build helpers: `Rscript`, `pandoc`, `gh`, `gcc`, `g++`, `make`, and
  `cmake`.
- Python tooling: `pip`, `pip3`, `pre-commit`, `ruff`, and `exiftool` are
  exposed through small shims in `C:\Users\jakob\.local\bin`, which is already
  on the user `PATH`.
- Use `magick`, not Windows `convert.exe`. Tesseract has `eng` and `osd`
  data; use `ocrmypdf --skip-text --deskew --rotate-pages` for scanned PDFs
  and `pdftotext -layout` for searchable PDFs.

## Working Rules

- Keep ado/Mata/plugin behavior aligned with the documented command contracts
  in `README.md`.
- Use `bash tools/ci/pkg_check.sh` for static package integrity checks.
- Use `bash tools/ci/archive_check.sh archives/source.zip` to validate the
  generated source archive when checking packaging output.
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
  On this Windows machine, the expected StataNow 19 executable is
  `C:\Program Files\StataNow19\StataSE-64.exe`.
