# Pre-merge testing audit

Audit `feature/clate` (starting 5744599) against production `origin/master`.
Scope: CLATE, changes shared by CATE/NP, postestimation, summaries, numerical
R parity, documentation/examples, installation/source archive, native paths
and test-harness reliability. R remains the numerical reference. Production
integration itself is outside this audit; changes stay on the feature branch.

- [x] Refresh production reference and inventory changed interfaces/tests.
- [x] Independently inspect paired estimation/covariance and broaden fixed/R-generated edge fixtures.
- [x] Independently inspect postestimation state, selectors/returns/rendering and broaden behavioral regressions.
- [x] Independently inspect test oracles/runner/packaging/CI and test false-pass/failure paths.
- [x] Parent: compare actual production/feature fits, restored production estimates, regenerate R fixtures reproducibly, run full licensed suite, isolated installation and source-archive/doc workflows.
- [x] Build/test native plugin/core paths on available platforms, report untested platform/version limits precisely.
- [x] Reproduce defects before fixes, verify corrected paths, run final affected/full gates on final snapshot.
- [x] Record findings, coverage, commands/results, remaining risks and evidence-based merge readiness in docs/premerge-testing-audit.md.

Three independent workers own separate new audit tests/fixture generators;
the parent owns production fixes, shared runner/manifests and final evidence.
Do not adjust R goldens to match a Stata error. Existing uncommitted changes
are user-owned (none at audit start). Keep raw logs/build outputs ignored.
