# Documentation parity across the Stata package

The user requests documentation comparable to the live R package for all
Stata functionality, on feature/clate, followed by commit and push. Baseline
ca690eb. This is a documentation/examples task; numerical interfaces and
estimators remain unchanged. The R sibling supplies the topic structure and
scientific reference; Stata source supplies exact implemented syntax/results.

Deliver a navigable article collection covering onboarding, notation, CATE
and regression, weights and tuning, diagnostics/reporting, visualization,
CLATE and PLQR moments. Add a package-level installed help entry point, improve
all applicable command helps, and supply deterministic runnable examples.
Map every R article/reference family to its Stata coverage or explicitly state
that the feature is R-only. Preserve the distinctions between estimate-only
and inference paths, pointwise and joint/simultaneous inference, ratio status
and weak-IV claims, NP outcome distributions and confidence intervals, and
PLQR moments and coefficient inference. Privilege preferred Stata options and
show compatibility aliases as such. Explain sample/target/matrix conventions,
returned results, selection/saving/export, installation and native fallback.

Public documentation must work from a standalone Stata checkout/install,
without sibling filesystem paths. These development docs target feature/clate;
release0.1.0 remains distinguished. Preserve prior verification receipts and
historical plans as dated records. No website deployment or release is required.

Acceptance: source-grounded topic/command coverage audit, all runnable examples
executed with licensed Stata (including installed overview help), local links,
help anchors/syntax and manifest checked, static/Python/source archive and
fresh isolated installation passed. No full numerical rerun for docs-only
changes unless behavior changes or a relevant failure requires one. Obtain
an independent documentation review, then commit and push feature/clate.
