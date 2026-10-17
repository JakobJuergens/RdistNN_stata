{smcl}
{* *! version 0.1.2  2026-10-17}{...}
{title:Title}

{phang}{bf:np_dnn_diagnostics} {hline 2} Inspect NP localization and outcome diagnostics

{title:Syntax}

{p 8 17 2}{cmd:np_dnn_diagnostics}
[{cmd:,} {opt targets(numlist)} {opt candidates(numlist)} {opt allcandidates} {opt detail}
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

{pstd}
The command preserves the active estimation results and dataset. The full
selected table is returned even when the compact display is used.

{title:Options}

{phang}{opt targets(numlist)} selects original target row numbers in the
supplied order. It cannot be combined with {opt candidates()} or
{opt allcandidates}.{p_end}
{phang}{opt candidates(numlist)} selects candidate row numbers from
{cmd:e(tuning_candidates)} after {cmd:np_dnn, tuningmode(grid)}.{p_end}
{phang}{opt allcandidates} selects every candidate row, including infeasible
rows. It cannot be combined with {opt candidates()}.{p_end}
{phang}{opt detail} displays the full table rather than a compact selection.{p_end}
{phang}{opt saving(filename)} saves the full selected table as a Stata dataset.{p_end}
{phang}{opt replace} permits overwriting the saved dataset.{p_end}

{title:Stored results}

{synoptset 20}{...}
{synopt:{cmd:r(diagnostics)}}full selected diagnostics matrix{p_end}
{synopt:{cmd:r(localization)}}common localization matrix (target scope only){p_end}
{synopt:{cmd:r(compact)}}compact displayed matrix{p_end}
{synopt:{cmd:r(n_rows)}}number of returned rows{p_end}
{synopt:{cmd:r(scope)}}{cmd:targets} or {cmd:candidates}{p_end}
{synopt:{cmd:r(selected)}}selected original row numbers in returned order{p_end}
{synopt:{cmd:r(saved)}}saved filename, if any{p_end}

{pstd}
{cmd:r(localization)} contains {cmd:target_id}, target coordinates
{cmd:x0_1}, ..., {cmd:x0_p}, then {cmd:ess}, {cmd:max_weight},
{cmd:top1_weight_share}, {cmd:top5_weight_share}, {cmd:top10_weight_share},
{cmd:n_nonzero_weight}, {cmd:distance_available}, {cmd:radius_rank_k},
{cmd:radius_rank_s}, {cmd:weighted_mean_distance}, and
{cmd:weighted_median_distance}. This schema is shared with CATE and CLATE.
IDs and row labels retain their original target identities after selection.

{title:Examples}

{phang2}{cmd:. np_dnn_diagnostics}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, detail}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, targets(3 1)}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, allcandidates}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, saving(np_diagnostics.dta) replace}{p_end}
