# Stata documentation parity plan

> Use executing-plans for parent work and dispatching-parallel-agents for the independent article/help domains. Track evidence in build/documentation/.

**Goal:** Bring all Stata user documentation to the coverage and usability of the R package, commit and push feature/clate.
**Architecture:** Task-oriented articles plus installed command references and deterministic examples; one R-to-Stata coverage map exposes supported and R-only surfaces.
**Tech Stack:** Markdown, SMCL help, Stata17 do files, Python/shell, licensed WindowsStataNow19.
**Spec:** ../specs/2026-10-17-documentation-parity-design.md

## Constraints and review focus

Documentation/examples only; preserve numerical code/APIs. Cover every public
Stata command and R article/reference family. Explain unavailable diagnostics,
missing targets/SEs, confidence levels and inference boundaries accurately.
Validate reader workflows in both checkout and clean installation; no sibling
paths in public user instructions. Work on the existing feature/clate branch.

- [x] Core articles and CATE/NP help/examples: docs/getting-started.md, docs/cate-and-regression.md, examples/cate-and-regression.do, core CATE/NP helps.
- [x] Reporting/plots articles, six postestimation helps and examples: docs/diagnostics-and-reporting.md, docs/visualising-fits.md, examples/diagnostics-and-plots.do.
- [x] Experimental articles/core helps/examples: docs/conditional-late.md, docs/plqr-moments.md, examples/clate.do and examples/plqr-moments.do.
- [x] Parent: article index, notation, weights/tuning, full R documentation coverage map, installed rdistnn overview, weights help/example, README/command guide/developer routing, manifest and doc checker.
- [x] Audit all topic/option/result boundaries; execute all examples, check links/help/install/archive, independent review, verification receipt.
- [x] Commit and push feature/clate.

Acceptance evidence: [documentation verification](../../documentation-verification.md).
