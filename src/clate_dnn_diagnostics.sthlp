{smcl}
{* *! version 0.1.0  2026-10-17}{...}
{viewerjumpto "Syntax" "clate_dnn_diagnostics##syntax"}{...}
{viewerjumpto "Description" "clate_dnn_diagnostics##description"}{...}
{viewerjumpto "Options" "clate_dnn_diagnostics##options"}{...}
{viewerjumpto "Stored results" "clate_dnn_diagnostics##results"}{...}
{viewerjumpto "Examples" "clate_dnn_diagnostics##examples"}{...}

{title:Title}

{phang}
{bf:clate_dnn_diagnostics} {hline 2} Inspect CLATE weight, distance, overlap,
and tuning diagnostics

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:clate_dnn_diagnostics}
[{cmd:,} {opt targets(numlist)} {opt candidates(numlist)} {opt allcandidates} {opt detail}
{opt saving(filename)} {opt replace}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:clate_dnn_diagnostics} is a postestimation command for {help clate_dnn}.
By default it reports the R-aligned per-target diagnostics stored by the
active fit: effective sample size, weight concentration, localization radii,
weighted distance, instrument-overlap summaries, the ordinary-Euclidean radius
containing 90 percent of DNN weight, normalized within-arm instrument-1/instrument-0
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
Arm ESS is available only when {cmd:clate_dnn} received nuisance-based instrument
inputs.  It is missing when an arm has zero local weight mass.  Physical radius
and score-contribution concentration are available for both nuisance-based and
precomputed-score fits.

{pstd}
The 25-column {cmd:e(diagnostics)} schema is shared with CATE: its internal
{cmd:treated_*}/{cmd:control_*} arm slots refer here to A=1/A=0. This command
renames the displayed and saved arm fields {cmd:instrument_1_*} and
{cmd:instrument_0_*}. Arm summaries are unavailable for paired-signal fits.
The score-contribution concentration uses the paired score
{it:phi_Y-theta*phi_D}.

{pstd}
For target diagnostics the command also displays and returns reduced form,
first stage, their standard errors, covariance and correlation, the signed
first-stage estimate divided by its standard error, and the ratio estimate.
Zero standard errors give missing correlation/statistic. The five status
columns preserve {cmd:e(status)}: ratio_status (0 ok, 1 zero denominator,
2 numerical failure), negative_denominator, near_zero_denominator,
covariance_status (0 ok, 1 not computed, 2 zero score variance,
3 numerical failure), and wald_available. A near-zero flag does not change
the ratio. These are pointwise Wald calculations without weak-IV-uniform
guarantees.

{marker options}{...}
{title:Options}

{phang}
{opt targets(numlist)} selects original target row numbers. It cannot be
combined with candidate diagnostics.

{phang}
{opt candidates(numlist)} selects tuning-candidate row numbers.  Active
estimates must come from {cmd:clate_dnn, tuningmode(grid)}.

{phang}
{opt allcandidates} reports every stored tuning candidate, including
infeasible rows.

{phang}
{opt detail} displays the full table rather than the compact diagnostic
columns.

{phang}
{opt saving(filename)} saves the full diagnostic table as a Stata dataset,
appending component and status columns for target scope.
{opt replace} permits overwriting an existing file.

{marker results}{...}
{title:Stored results}

{synoptset 22 tabbed}{...}
{synopt:{cmd:r(diagnostics)}}full selected diagnostics matrix{p_end}
{synopt:{cmd:r(components)}}target_id, reduced_form, first_stage, numerator_se,
denominator_se, component_covariance, component_correlation,
denominator_stat, estimate (target scope only){p_end}
{synopt:{cmd:r(status)}}five status columns for selected targets{p_end}
{synopt:{cmd:r(compact)}}compact displayed matrix{p_end}
{synopt:{cmd:r(n_rows)}}number of reported rows{p_end}
{synopt:{cmd:r(scope)}}{cmd:targets} or {cmd:candidates}{p_end}
{synopt:{cmd:r(saved)}}saved filename, if supplied{p_end}

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. clate_dnn_diagnostics}{p_end}
{phang2}{cmd:. clate_dnn_diagnostics, detail}{p_end}
{phang2}{cmd:. clate_dnn_diagnostics, allcandidates}{p_end}
{phang2}{cmd:. clate_dnn_diagnostics, saving(cate_diagnostics.dta) replace}{p_end}
