{smcl}
{* *! version 0.1.0  2026-07-12}{...}
{title:Title}

{phang}{bf:np_dnn_diagnostics} {hline 2} Inspect NP localization and outcome diagnostics

{title:Syntax}

{p 8 17 2}{cmd:np_dnn_diagnostics}
[{cmd:,} {opt candidates(numlist)} {opt allcandidates} {opt detail}
{opt saving(filename)} {opt replace}]

{title:Description}

{pstd}
{cmd:np_dnn_diagnostics} is a postestimation command for {help np_dnn}.
By default it reports one row per target with target coordinates, weight
concentration, distance scales, and localized outcome summaries.  Specify
{opt candidates()} or {opt allcandidates} after a grid fit to report
candidate-level tuning diagnostics.
The built-in sensitivity grid is lambda = (.3, .4, .5, .75), with .4 marked
as the finite-sample default candidate.

{title:Options}

{phang}{opt candidates(numlist)} selects candidate row numbers from
{cmd:e(tuning_candidates)}.{p_end}
{phang}{opt allcandidates} selects every candidate row.{p_end}
{phang}{opt detail} displays the full table rather than a compact selection.{p_end}
{phang}{opt saving(filename)} saves the full selected table as a Stata dataset.{p_end}
{phang}{opt replace} permits overwriting the saved dataset.{p_end}

{title:Stored results}

{synoptset 20}{...}
{synopt:{cmd:r(diagnostics)}}full selected diagnostics matrix{p_end}
{synopt:{cmd:r(compact)}}compact displayed matrix{p_end}
{synopt:{cmd:r(n_rows)}}number of returned rows{p_end}
{synopt:{cmd:r(scope)}}{cmd:targets} or {cmd:candidates}{p_end}
{synopt:{cmd:r(saved)}}saved filename, if any{p_end}

{title:Examples}

{phang2}{cmd:. np_dnn_diagnostics}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, detail}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, allcandidates}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, saving(np_diagnostics.dta) replace}{p_end}
