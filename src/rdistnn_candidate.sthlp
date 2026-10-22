{smcl}
{* *! version 0.1.0  2026-10-18}{...}
{title:rdistnn_candidate — Extract a reusable candidate fit}
{p 8 12 2}{cmd:rdistnn_candidate, candidate(}{it:#}{cmd:)}
[{cmd:level(}{it:#}{cmd:) post store(}{it:name}{cmd:)}]{p_end}
{pstd}
Select an available original row of e(tuning_candidates) from NP, CATE or CLATE.
The command reuses saved ranks, distances, frozen signals and the fitted sample
at the candidate's resolved s. It does not enumerate subsamples or rebuild
geometry. Duplicate resolved orders give identical fits while retaining the
requested candidate ID. The chosen fit has no nested tuning grid.
{pstd}
Weights are reconstructed from saved raw covariates as well as ranks. Exactly
identical raw covariate rows receive their block's average weight, and each
deletion redistributes the n-1 table's block mass over the surviving members.
Commands warn that this modification is outside the paper's theory. Refit
older saved results without e(model_version)=2 before extracting candidates;
see {help weights_dnn} for the shared block rule.
{pstd}
The parent's covariance policy is retained, including nose. post replaces the
active fit; store(name) retains a named standalone fit. By default e() and the
active data are preserved. Errors preserve the active fit. Old results without
a numeric model snapshot require a refit. PLQR has no integrated Stata candidate
grid and therefore cannot supply a candidate view.
{pstd}
Current snapshots are validated before reconstruction, including dimensions,
controls, identities, fitted geometry and cached ordering. Refit an incomplete
or inconsistent snapshot. Validation errors preserve the active data and e().
{title:Returned results}
{pstd}
Returns r(candidate), r(s), and the reporting/prediction matrices documented
in {help rdistnn_results}. r(stored) names an explicitly stored fit. level()
controls normal pointwise limits and uses percentages. A stored or posted fit
can be replayed, diagnosed, plotted, saved and predicted independently of its
parent; its e(candidate_id) identifies the requested original grid row.
{title:See also}
{pstd}{help rdistnn}, {help rdistnn_results}, {help rdistnn_predict}{p_end}
