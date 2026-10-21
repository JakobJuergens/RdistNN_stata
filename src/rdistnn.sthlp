{smcl}
{* *! version 0.1.1  2026-10-21}{...}
{viewerjumpto "Commands" "rdistnn##commands"}{...}
{viewerjumpto "Conventions" "rdistnn##conventions"}{...}
{viewerjumpto "Workflow" "rdistnn##workflow"}{...}
{viewerjumpto "Inference" "rdistnn##inference"}{...}
{viewerjumpto "Articles" "rdistnn##articles"}{...}

{title:Title}

{phang}
{bf:RdistNN for Stata} {hline 2} Distributional nearest-neighbor estimation,
weights, diagnostics and native plots

{pstd}
Stata 17 or later. Package {bf:0.1.1} provides conditional treatment
estimation, nonparametric and two-scale regression, and shared reporting.
RdistNN is the numerical reference for shared functionality; the Stata
companion has its own documented public surface.

{marker commands}{...}
{title:Commands}

{synoptset 28}{...}
{synopt:{help cate_dnn}}CATE using AIPW or supplied scores{p_end}
{synopt:{help np_dnn}}nonparametric conditional means{p_end}
{synopt:{help weights_dnn}}single/batch exact/asymptotic normalized weights{p_end}
{synopt:{help np_tdnn}}two-scale regression with signed weights{p_end}
{synopt:{help weights_tdnn}}single/batch signed two-scale rank weights{p_end}
{synopt:{help np_tdnn_diagnostics}}signed-weight mass and concentration tables{p_end}
{synopt:{help np_tdnn_plot}}two-scale estimates, signed weights and diagnostics{p_end}
{synopt:{help clate_dnn}}binary-IV or paired-signal ratios{p_end}
{synopt:{help plqr_dnn}}experimental moments at supplied PLQR coefficients{p_end}
{synopt:{help cate_dnn_diagnostics}}CATE target/candidate tables and saving{p_end}
{synopt:{help np_dnn_diagnostics}}regression localization/outcome tables{p_end}
{synopt:{help clate_dnn_diagnostics}}instrument/component/status tables{p_end}
{synopt:{help cate_dnn_plot}}CATE estimates, diagnostics, overlap, influence and weights{p_end}
{synopt:{help np_dnn_plot}}regression panels including outcome distributions{p_end}
{synopt:{help clate_dnn_plot}}CLATE panels with instrument and paired-score meaning{p_end}

{synopt:{help rdistnn_results}}graph-free numeric reports and exports{p_end}
{synopt:{help rdistnn_candidate}}reusable candidate fit from saved ranks{p_end}
{synopt:{help rdistnn_predict}}evaluate saved frozen signals at new targets{p_end}

{marker conventions}{...}
{title:Conventions}

{pstd}
The active dataset supplies numeric observations and nuisance prediction columns.
CATE, CLATE and NP use {cmd:xcols()} plus exactly one of {cmd:x0()} or
{cmd:x0matrix()}. The target matrix has one row per target and one column per
covariate, in xcols order. Standalone weights take covariates positionally and
one x0 or a target matrix; PLQR takes x(), d(), g(), md() and one x0.
None fits nuisance learners.

{pstd}
The default nearest-neighbor count is k=1. Kernel order s is the subset size,
not a request to keep only the first s full-sample observations. Exact support
ends at rank n-s+k. Larger s tightens localization; ESS=1/sum(w_i^2) measures
weight concentration. Exact/asymptotic methods and L1/L2, scaled L1/L2 and Mahalanobis metrics are
available subject to each command's sample/dimension checks.

{pstd}
{help np_tdnn} combines two orders s1<s2 for conditional mean regression.
Its weights can be negative. Use its signed-weight diagnostics rather than
probability-weighted outcome quantiles or cumulative probability shares.
Both orders, normalization and fitted geometry remain fixed in deletions and
saved predictions. CLATE and CATE retain single-scale DNN.

{pstd}
CATE/CLATE/NP support batch targets, lambda-based defaults and integrated
sensitivity grids. A grid retains a chosen default rather than optimizing it.
Standalone weights and PLQR have explicit s controls without integrated grids.
Covariance needs s<n; CATE/CLATE/NP {cmd:nose} and standalone weights permit s=n.

{marker workflow}{...}
{title:Workflow}

{pstd}On a dataset containing Y, x1 and x2:{p_end}
{phang2}{cmd:. np_dnn Y, xcols(x1 x2) x0(0 0)}{p_end}
{phang2}{cmd:. matrix estimates = e(estimates)}{p_end}
{phang2}{cmd:. np_dnn, level(90)}{p_end}
{phang2}{cmd:. matrix summary = r(summary)}{p_end}
{phang2}{cmd:. np_dnn_diagnostics, detail}{p_end}
{phang2}{cmd:. np_dnn_plot, name(local_mean) replace}{p_end}

{pstd}
Fit/replay of CATE/CLATE/NP returns full-target {cmd:r(table)} (9 by m) and
{cmd:r(summary)} (m by 6). Canonical {cmd:e(estimates)}, {cmd:e(se_targets)}
and {cmd:e(estimate_vcov)} preserve original target order, including undefined
CLATE ratios. Copy r() results before another command overwrites them; use
{cmd:estimates store}/{cmd:estimates restore} for multiple fits.

{pstd}
Diagnostic selectors preserve original target/candidate IDs and can save .dta
reports. Single-view plot commands return selected {cmd:r(plot_data)} and
{cmd:r(n_rows)}; dashboards return panel metadata. Effective-kernel and influence
views use the current numeric snapshot. Older development fits require refitting.
Stored-result views can be used after data changes. See each command's help
for allowed selectors, view-specific options and result schemas.
Use {help rdistnn_results} for numeric extraction without graphs,
{help rdistnn_candidate} for standalone grid fits, and {help rdistnn_predict}
for saved-fit prediction. estimates save/use retains the numeric model snapshot.
PLQR returns moments with supplied coefficients reported separately.

{marker inference}{...}
{title:Inference and experimental scope}

{pstd}
Reported normal intervals are pointwise. Joint target covariance does not
supply simultaneous coverage. Frozen-score jackknife calculations do not
refit nuisances, remove localization bias or verify the design, nuisance-rate,
delete-one stability, variance or undersmoothing conditions. Supplied predictions
may use a fuller adjustment state than the variables defining localization.
For identification, nuisance-learning and localization conditions, see the
{browse "https://jakobjuergens.com/research/#orthogonal-inference-conditional-z":paper}
and the estimator-specific help.

{pstd}
CLATE additionally needs the documented binary-IV causal assumptions for a
conditional complier interpretation. Negative or near-zero first-stage signals
are retained with statuses; pointwise Wald inference is not weak-IV-uniform.
PLQR evaluates score moments at thetaeval and supplies moment covariance;
it does not solve coefficients or provide coefficient inference. Weighted NP
outcome quantiles describe the outcome distribution, not a confidence interval.

{marker articles}{...}
{title:Articles, examples and installation}

{pstd}
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/index.md":Article index}
links getting started, notation, CATE/regression, weights/tuning, diagnostics,
visualization, CLATE and PLQR guides with deterministic standalone do-files.
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/r-documentation-crosswalk.md":R documentation crosswalk}
maps every R article/reference family and identifies R-only GMM,
custom registries, confidence-set inversion and replication workflows.

{pstd}
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/README.md":README}
provides installation from {cmd:master};
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/RELEASING.md":RELEASING.md}
provides an immutable installation pin. The portable ado/Mata baseline needs
no compiler. Optional native plugins accelerate distance and frozen-score
jackknife calculations. See the
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/src/plugin/README.md":plugin guide}
for platform binaries and engine controls. This overview and command help
install with the package; articles and examples are available from the
repository or source ZIP.

{title:Numerical implementation}

{pstd}
Compensated sums, numerical centering and scaled covariance calculations
protect represented results from rounding and intermediate range failures.
The {browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/numerical-stability.md":numerical implementation reference}
describes their behavior and limits. Exact constants retain zero variance;
unavailable ratio uncertainty retains its numerical status. These safeguards
do not verify the statistical assumptions for inference.

{title:Contributors}

{pstd}
Quality assurance, testing and validation:
{browse "https://github.com/Zakareya-Chehime":Zakareya Chehime}.
