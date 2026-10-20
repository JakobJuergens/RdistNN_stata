{smcl}
{* *! version 0.1.0.9004  2026-10-20}{...}
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
Stata 17 or later. This overview describes package {bf:0.1.0.9004} in the development snapshot,
including supported supplied-prediction conditional LATE, two-scale regression
and shared reporting. These are development-branch additions.
RdistNN is the numerical reference for shared
functionality; the Stata companion has its own documented public surface.

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
views use the numeric snapshot in current fits. Legacy fits require the original
fitted values, sample and observation order.
Stored-result views can be used after data changes. See each command's help
for allowed selectors, view-specific options and result schemas.
Use {help rdistnn_results} for numeric extraction without graphs,
{help rdistnn_candidate} for standalone grid fits, and {help rdistnn_predict}
for saved-fit prediction. estimates save/use retains the numeric model snapshot.
PLQR returns moments with supplied coefficients reported separately.

{marker inference}{...}
{title:Inference and experimental scope}

{pstd}
For numerical stability, all estimator jackknife engines, including naive
deletions, subtract a common reference per target and signal before summing.
Only observations with nonzero weight in at least one deletion determine the
reference. A rank-r weight in the n-1 table can use original ranks r and r+1;
the latter can gain weight when an earlier observation is deleted. Exact DNN
deletion support therefore ends at n-s+k (n-s+1 for k=1). Unused observations
are excluded before subtraction, preventing zero-weight extremes from
affecting the reference or causing overflow. All n deletions and their
original covariance scaling remain in the calculation.

{pstd}
Weighted signal means and rank-weight normalization use compensated sums.
Native and Mata deletion engines retain correction terms until the final
prefix/suffix combination. The fallback uses bounded merge buffers, and
candidate fits and saved prediction reuse the same arithmetic helpers.
These double-precision safeguards preserve the estimator's statistical definition.

{pstd}
The reference is the smallest active value if all are positive, the largest
if all are negative, and zero otherwise. This common numerical shift cancels
in jackknife covariance; point estimates use the original signals. It is
separate from the usual centering of deletion estimates around their mean.
Asymptotic weights retain every represented nonzero entry, without the exact
support cutoff or additional renormalization. See the
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/weights-and-tuning.md#numerical-centering":numerical centering guide}
for the shared R/Stata variance-9 example.

{pstd}
Reported normal intervals are pointwise. Joint target covariance does not
supply simultaneous coverage. Frozen-score jackknife calculations do not
refit nuisances, remove localization bias or verify the design, nuisance-rate,
delete-one stability, variance or undersmoothing conditions. Supplied predictions
may use a fuller adjustment state than the variables defining localization.
Consult the estimator-specific inference contract before interpreting results.

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
provides installation from {cmd:master} and version-pinned alternatives. The portable
ado/Mata baseline needs no compiler. Optional native plugins accelerate
L1/L2 distances, including fitted scales and batches, through the common
weight, estimation and saved-prediction path. Scale fitting, stable ranks and
weight tables and Mahalanobis factorization remain in Mata. Projected
Mahalanobis norms use the selected L2 engine. Legacy saved inverses retain
their original Mata calculation.
Use {cmd:global rdistnn_distance_engine plugin} or {cmd:mata} for explicit
distance routes; drop the global to restore automatic dispatch, which keeps
small jobs in Mata and attempts native distances for substantial batches. Plugins also
accelerate frozen-score jackknife deletions for all estimators.
Both C++ and Mata use prefix/suffix sums; naive deletions remain a final
fallback and can be selected with {cmd:global rdistnn_jackknife_engine naive}
for parity checks. This overview and command help install with the package;
the article/example collection is available from the repository or source ZIP.

{title:Numeric covariance range}

{pstd}
Shared covariance helpers normalize exceptional columns by binary powers
before centering or squaring and preserve full paired/target dependence.
Representable extreme variances are recovered after intermediate range
failures. A nonzero variance outside the finite numeric range gives an error
requesting rescaling or point-only estimation; CLATE retains its per-target
numerical-failure status. Exact constants still have zero variance. Stata's
finite nonmissing range is smaller than R's IEEE range. These guards do not
change the documented inference conditions.
