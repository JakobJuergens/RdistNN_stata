{smcl}
{* *! version 0.1.1  2026-07-12}{...}
{viewerjumpto "Syntax" "cate_dnn_diagnostics##syntax"}{...}
{viewerjumpto "Description" "cate_dnn_diagnostics##description"}{...}
{viewerjumpto "Options" "cate_dnn_diagnostics##options"}{...}
{viewerjumpto "Stored results" "cate_dnn_diagnostics##results"}{...}
{viewerjumpto "Examples" "cate_dnn_diagnostics##examples"}{...}

{title:Title}

{phang}
{bf:cate_dnn_diagnostics} {hline 2} Inspect CATE weight, distance, overlap,
and tuning diagnostics

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:cate_dnn_diagnostics}
[{cmd:,} {opt candidates(numlist)} {opt allcandidates} {opt detail}
{opt saving(filename)} {opt replace}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:cate_dnn_diagnostics} is a postestimation command for {help cate_dnn}.
By default it reports the R-aligned per-target diagnostics stored by the
active fit: effective sample size, weight concentration, localization radii,
weighted distance, treatment-overlap summaries, the ordinary-Euclidean radius
containing 90 percent of DNN weight, normalized within-arm treated/control
ESS, and top-five absolute score-contribution share.  The physical radius is
not the rank-s radius; arm ESS is not overall ESS; score-contribution share is
not top-five raw-weight share.  {opt candidates()} or {opt allcandidates}
instead reports tuning-grid diagnostics.
The default sensitivity grid is lambda = (.3, .4, .5, .75), with .4 marked
as the finite-sample default candidate.

{pstd}
The command does not modify the active estimation results or dataset.  The
full table is returned in {cmd:r(diagnostics)} even when the compact display
is used.

{pstd}
Arm ESS is available only when {cmd:cate_dnn} received nuisance-based treatment
inputs.  It is missing when an arm has zero local weight mass.  Physical radius
and score-contribution concentration are available for both nuisance-based and
precomputed-score fits.

{marker options}{...}
{title:Options}

{phang}
{opt candidates(numlist)} selects tuning-candidate row numbers.  Active
estimates must come from {cmd:cate_dnn, tuningmode(grid)}.

{phang}
{opt allcandidates} reports every stored tuning candidate, including
infeasible rows.

{phang}
{opt detail} displays the full table rather than the compact diagnostic
columns.

{phang}
{opt saving(filename)} saves the full diagnostic table as a Stata dataset.
{opt replace} permits overwriting an existing file.

{marker results}{...}
{title:Stored results}

{synoptset 22 tabbed}{...}
{synopt:{cmd:r(diagnostics)}}full selected diagnostics matrix{p_end}
{synopt:{cmd:r(compact)}}compact displayed matrix{p_end}
{synopt:{cmd:r(n_rows)}}number of reported rows{p_end}
{synopt:{cmd:r(scope)}}{cmd:targets} or {cmd:candidates}{p_end}
{synopt:{cmd:r(saved)}}saved filename, if supplied{p_end}

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. cate_dnn_diagnostics}{p_end}
{phang2}{cmd:. cate_dnn_diagnostics, detail}{p_end}
{phang2}{cmd:. cate_dnn_diagnostics, allcandidates}{p_end}
{phang2}{cmd:. cate_dnn_diagnostics, saving(cate_diagnostics.dta) replace}{p_end}
