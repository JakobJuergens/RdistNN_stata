{smcl}
{* *! version 0.1.1  17oct2026}{...}
{viewerjumpto "Syntax" "clate_dnn##syntax"}{...}
{viewerjumpto "Signals and targets" "clate_dnn##description"}{...}
{viewerjumpto "Options" "clate_dnn##options"}{...}
{viewerjumpto "Paired inference" "clate_dnn##inference"}{...}
{viewerjumpto "Numerical status" "clate_dnn##status"}{...}
{viewerjumpto "Stored results" "clate_dnn##results"}{...}
{viewerjumpto "Examples" "clate_dnn##examples"}{...}

{title:Title}

{phang}
{cmd:clate_dnn} {hline 2} Conditional LATE from supplied binary-IV predictions

{pstd}
Stata 17 or later is required. Supplied-prediction CLATE and ordinary paired-signal
ratios are supported workflows in the development version {bf:0.1.0.9004} under the documented
pointwise inference conditions. See {help rdistnn} for package navigation.

{marker syntax}{...}
{title:Syntax}

{pstd}Supplied instrument-arm nuisance predictions:

{p 8 12 2}
{cmd:clate_dnn} {it:Y} {ifin}, {cmd:xcols(}{it:varlist}{cmd:)}
{it:target_specification}
{cmd:dcol(}{it:D}{cmd:)} {cmd:acol(}{it:A}{cmd:)} {cmd:ecol(}{it:e}{cmd:)}
{cmd:my0col(}{it:my0}{cmd:)} {cmd:my1col(}{it:my1}{cmd:)}
{cmd:md1col(}{it:md1}{cmd:)}
[{cmd:md0col(}{it:md0}{cmd:)} {cmd:onesided} {it:options}]

{pstd}Two supplied signals:

{p 8 12 2}
{cmd:clate_dnn} {it:phi_y} {ifin}, {cmd:pseudooutcome}
{cmd:denominator(}{it:phi_d}{cmd:)} {cmd:xcols(}{it:varlist}{cmd:)}
{it:target_specification} [{it:options}]

{pstd}
{it:target_specification} is exactly one of {cmd:x0(}{it:numlist}{cmd:)} and
{cmd:x0matrix(}{it:matrix_name}{cmd:)}. Targets occupy matrix rows;
columns follow {cmd:xcols()} order. {cmd:md0col()} is required unless
{cmd:onesided} is specified. {cmd:x()} is a compatibility alias for
{cmd:xcols()}; do not combine them.

{pstd}Replay, with a confidence level expressed as a percentage:

{p 8 12 2}{cmd:clate_dnn} [{cmd:, level(}{it:#}{cmd:)}]

{synoptset 26 tabbed}{...}
{synopthdr:options}
{synoptline}
{syntab:Signals}
{synopt:{opt clip(#)}}instrument-propensity clipping; default 1e-8{p_end}
{synopt:{opt denominatortol(#)}}numerical first-stage cancellation flag; default about 1.49e-8{p_end}
{synopt:{opt training(string)}}supplied, oracle, full_sample, oob or cross_fit; default supplied{p_end}
{syntab:Localization and reporting}
{synopt:{opt k(#)}}neighbors averaged within each size-s subsample; default 1{p_end}
{synopt:{opt s(#)}}deterministic integer subsample order{p_end}
{synopt:{opt lambda(#)}}paper-rule multiplier; default .4{p_end}
{synopt:{opt tuningmode(string)}}default or grid{p_end}
{synopt:{opt lambdagrid(numlist)}}replace grid multipliers; default .3 .4 .5 .75{p_end}
{synopt:{opt sgrid(numlist)}}direct integer candidate orders, replacing multiplier grid{p_end}
{synopt:{opt sdefault(#)}}default order for sgrid(), not a candidate row number{p_end}
{synopt:{opt distance(string)}}l1, l2 (default), l1_scaled, l2_scaled, or mahalanobis{p_end}
{synopt:{opt scalerule(string)}}sd (default), mad, iqr, range, or custom{p_end}
{synopt:{opt scalevec(numlist)}}custom scales, one per covariate{p_end}
{synopt:{opt scaleeps(#)}}positive scale floor; default 1e-12{p_end}
{synopt:{opt method(string)}}exact (default) or asymp{p_end}
{synopt:{opt nose}}estimates without jackknife covariance{p_end}
{synopt:{opt level(#)}}pointwise normal confidence level; defaults to {cmd:c(level)}, normally 95{p_end}
{synoptline}

{marker description}{...}
{title:Signals and targets}

{pstd}
A is a binary instrument and D a binary treatment. Supplied predictions
condition on all pre-instrument adjustment variables V=(X,Z), while weights
use only the reporting covariates in {cmd:xcols()}. The predictions are
e(V)=Pr(A=1|V), m_ya(V)=E[Y|A=a,V] and m_da(V)=E[D|A=a,V]. Arm indices
refer to A. Outcome predictions do not additionally condition on realized D.
The command fits no nuisance learners.

{pstd}
For each response R=Y,D it constructs the augmented signal
phi_R=m_R1-m_R0+A*(R-m_R1)/e_used-(1-A)*(R-m_R0)/(1-e_used), with
e_used=min(max(e,clip),1-clip). Both signals receive the same DNN weights.
The estimate is Delta_Y(x)/Delta_D(x), where Delta_R=sum_i w_i(x)*phi_Ri.
It solves the affine equation sum_i w_i(x)*(theta*phi_Di-phi_Yi)=0.
It never averages observation-level signal ratios.

{pstd}
The target averages both instrument-arm contrasts over the same conditional
law of Z given X=x. Under the IV assumptions it is a compliance-weighted
conditional effect, generally different from population CATE, an unweighted
average of stratum effects, or residualized IV with instrument-variance weights.
Signed first-stage signals and finite negative aggregate denominators remain
unchanged. Augmented first-stage signals are not probabilities.

{marker options}{...}
{title:Options and sample}

{phang}
{cmd:xcols()} lists numeric reporting covariates. Supply one finite target
coordinate per covariate in {cmd:x0()}, or a nonmissing m-by-p matrix in
{cmd:x0matrix()}. All targets use the same complete estimation sample.

{phang}
On the nuisance route, raw e, md0 and md1 must be finite in [0,1] throughout
the requested if/in sample, before exclusions for other missing variables.
Missing outcome, covariate, instrument, treatment or outcome-prediction values
are excluded. D and A must be 0/1 on the resulting complete sample. External
predictions must already be correctly aligned with the observations; check
{cmd:e(N)} and {cmd:e(sample)}. Overflow in constructed signals stops estimation.

{phang}
{cmd:onesided} explicitly selects the structural restriction D(0)=0. It
requires observed D<=A, permits omitted {cmd:md0col()} interpreted as zero,
and requires a supplied md0 column to be zero. The observed inequality alone
does not establish the counterfactual restriction. Without {cmd:onesided},
all four instrument-arm conditional means are required.

{phang}
{cmd:pseudooutcome denominator()} accepts two aligned signals, requiring
complete signal/covariate rows. It cannot be combined with nuisance inputs or
{cmd:onesided}. This route also computes generic paired local ratios.
Instrument overlap information is unavailable, and no propensity clipping is
performed on supplied signals. There is no separate Stata affine-ratio command.

{phang}
{cmd:clip()} must be finite and strictly between zero and .5. Only e is
clipped; md0 and md1 remain unchanged. Clipping does not establish overlap.
{cmd:training()} records a declaration only. Prediction arrays do not verify
oracle, full-sample, out-of-bag or cross-fitted provenance.

{phang}
{cmd:k()} defaults to 1. The default subsample order is
floor(.4*n^((p+1)/(p+2))), with complete-sample n and reporting dimension p.
{cmd:lambda()} changes the .4 multiplier; {cmd:s()} overrides the order.
Exact top-k positive support has size n-s+k, or n-s+1 under k(1).
Larger s tightens localization; exact k=s gives uniform weights.
{cmd:method(asymp)} uses normalized Beta-rank weights.
{cmd:distance(mahalanobis)} holds full-estimation-sample covariance geometry
fixed across deletions, using generalized-inverse behavior for singular covariance.
Input observation position breaks distance ties.

{phang}
{cmd:tuningmode(grid)} compares multipliers .3 .4 .5 .75 while retaining the
ordinary .4 fit. {cmd:lambdagrid()} replaces the multiplier list and
{cmd:lambda()} identifies the default. {cmd:sgrid()} instead supplies integer
orders; its first order is the default unless {cmd:sdefault()} supplies an
order, inserted if absent. Grid options require grid mode. Do not combine
{cmd:s()} with grid mode, or {cmd:lambdagrid()} with {cmd:sgrid()}.
Infeasible nondefault candidates remain in metadata; an infeasible default
stops estimation. Rounded duplicate orders retain their candidate identities.
Grid comparisons are sensitivity analyses, without automatic order selection.
See {help cate_dnn} and {help weights_dnn} for the shared weight controls.

{phang}
Covariance requires 1<=k<=s<n. {cmd:nose} permits s=n and leaves covariance
and SEs missing. {cmd:level()} affects pointwise normal bounds on fit/replay,
without changing the fitted arrays. It does not create simultaneous bands.

{marker inference}{...}
{title:Inference contract}

{pstd}
The documented sampling framework uses independent observations; cluster
and survey-design inference are not supplied. Input row order breaks distance
ties, so permuting tied observations can change estimates.


{pstd}
Covariance uses paired frozen-score rank-shift deletions with the same s,
k, predictions, target points and metric. It retains both components and
all cross-covariances. Each ratio is linearized at its full-sample value.
Neither deleted ratios nor refitted nuisance learners are used. For one
target, V_theta=(V_YY+theta^2*V_DD-2*theta*V_YD)/Delta_D^2.
Dividing the numerator SE by the first stage omits denominator uncertainty.
{cmd:e(theta_vcov)} and {cmd:e(component_vcov)} are already on the estimator
scale: do not divide them by ESS. ESS is inverse weight concentration, not a
number of compliers. Cross-target covariance does not imply simultaneous coverage.

{pstd}
Causal interpretation requires consistency, exclusion, conditional instrument
independence given V, monotonicity, instrument overlap and positive local
compliance. Without those assumptions the calculation is a standardized
observed-data Wald ratio. Predictions alone cannot establish causal validity.

{pstd}
Orthogonality controls first-order population nuisance drift. Full-sample
inference additionally requires the affine theorem's design, oracle,
nuisance-accuracy, quantitative deletion-stability, envelope and exceptional-event
controls, positive score variance, fixed relevance and undersmoothing conditions.
The theorem's deleted nuisance states differ from the software's frozen-array
variance deletions. A training declaration alone does not establish these
conditions. The exact k=1 theorem does not automatically extend to every
supported k or weight method. Jackknife uncertainty does not remove target bias.

{pstd}
Intervals are pointwise Wald intervals, without weak-instrument-uniform
guarantees. There is no weak-stage trimming, first-stage cutoff or denominator
regularization. The R companion's separate {cmd:ratio_confidence_set()}
inversion is not implemented in Stata and its algebra alone does not
establish uniform weak-IV coverage.

{marker status}{...}
{title:Numerical status}

{pstd}
{cmd:e(status)} is m-by-5, in the following column order:

{synoptset 27 tabbed}{...}
{synopt:{cmd:ratio_status}}0=ok; 1=zero_denominator; 2=numerical_failure{p_end}
{synopt:{cmd:negative_denominator}}1 for a finite negative local first stage{p_end}
{synopt:{cmd:near_zero_denominator}}1 for the numerical cancellation condition below{p_end}
{synopt:{cmd:covariance_status}}0=ok; 1=not_computed; 2=zero_score_variance; 3=numerical_failure{p_end}
{synopt:{cmd:wald_available}}1 for available positive-variance ratio inference; 0 otherwise{p_end}

{pstd}
An exact zero aggregate denominator gives a missing ratio and no Wald bounds,
while components and their available covariance remain. A finite negative
first stage gives a signed ratio without replacement. Ratio and covariance
status must be read together: component covariance can be computed when no
ratio exists. Numerical failures remain distinct from exact zeros.
Finite component covariance can also remain available when ratio inference
is numerically unavailable.

{pstd}
{cmd:denominatortol()} must be finite and nonnegative. For nonzero first
stages, the flag checks abs(Delta_D)/max_i abs(phi_Di)<=tolerance. Its default
is sqrt(machine epsilon), about 1.49e-8; zero disables this flag. It changes
no ratio, covariance or Wald mask. It is not a statistical weak-IV threshold.

{pstd}
Exact zero paired-score variance is a computed zero SE, covariance status 2,
and no Wald statistic, p-value or interval. {cmd:nose} instead gives missing
uncertainty and covariance status 1. No zero coefficient substitutes for an
undefined ratio. All requested target positions remain in canonical tables.

{marker results}{...}

{title:Distance scaling}

{pstd}
{cmd:distance(l1)} uses Manhattan distance. {cmd:distance(l2)} uses Euclidean
distance (default). {cmd:l1_scaled} and {cmd:l2_scaled} divide each coordinate
by a scale fitted on the complete estimation sample before computing the
corresponding distance. {cmd:scalerule(sd)} is the default; alternatives are
{cmd:mad} (median absolute deviation, constant 1.4826), {cmd:iqr} (R type-7
quartiles), {cmd:range}, and {cmd:custom}. Sample SD uses the n-1 divisor and
requires two observations. Supply {cmd:scalevec(numlist)} only with
{cmd:scalerule(custom)}, with one finite scale per covariate. Custom values
are absolute scales, not inverse scales. Each scale is floored at the positive
finite {cmd:scaleeps()} (default 1e-12), including constant covariates and zero
custom scales. Scaling controls do not change unscaled or Mahalanobis metrics.
The fitted scales or Mahalanobis whitening factor remain fixed across
jackknife deletions and saved-model predictions. Ties retain original row order.

{title:Stored results}

{pstd}
With m targets, p reporting covariates and G tuning candidates:

{synoptset 27 tabbed}{...}
{p2col 5 27 31 2:Full-target matrices}{p_end}
{synopt:{cmd:e(targets)}, {cmd:e(x0_matrix)}}m-by-p targets, in input row order{p_end}
{synopt:{cmd:e(estimates)}}1-by-m ratios, including missing entries{p_end}
{synopt:{cmd:e(se_targets)}}1-by-m ratio SEs{p_end}
{synopt:{cmd:e(components)}}m-by-2, reduced_form first_stage{p_end}
{synopt:{cmd:e(theta_vcov)}}m-by-m joint ratio covariance{p_end}
{synopt:{cmd:e(estimate_vcov)}}copy of joint ratio covariance{p_end}
{synopt:{cmd:e(component_vcov)}}2m-by-2m joint covariance, interleaved y1,d1,y2,d2,...{p_end}
{synopt:{cmd:e(status)}}m-by-5 statuses described above{p_end}
{synopt:{cmd:e(neff_targets)}, {cmd:e(ess_targets)}}1-by-m effective sample sizes{p_end}
{synopt:{cmd:e(diagnostics)}}m-by-25 localization, overlap and score diagnostics{p_end}
{synopt:{cmd:e(target_diagnostics)}}copy of e(diagnostics){p_end}

{pstd}
Stata does not permit missing coefficients in {cmd:ereturn post}.
{cmd:e(b)} therefore contains only defined ratios (theta for one target;
theta1,theta2,... for a batch, preserving gaps in target numbers).
{cmd:e(V)} is posted only if covariance is available for every posted ratio.
All-undefined fits retain sample, components and full tables without
coefficients. Use canonical full matrices for batch export.

{p2col 5 27 31 2:Grid matrices}{p_end}
{synopt:{cmd:e(tuning_candidates)}}G-by-8 metadata: candidate_id computed_id is_default feasible reason status lambda s{p_end}
{synopt:{cmd:e(tuning_estimates)}, {cmd:e(tuning_se)}}G-by-m ratio and SE tables{p_end}
{synopt:{cmd:e(tuning_reduced_form)}}G-by-m reduced forms{p_end}
{synopt:{cmd:e(tuning_first_stage)}}G-by-m first stages{p_end}
{synopt:{cmd:e(tuning_target_status)}}(G*m)-by-7; candidate_id target_id then five status columns{p_end}
{synopt:{cmd:e(tuning_diagnostics)}}G-by-5 shared weight-concentration diagnostics{p_end}

{pstd}
Grid status rows are ordered by candidate then target; skipped candidates
have missing component/status entries. {cmd:e(tuning_default_candidate)}
indexes the retained fit. Scalars {cmd:e(tuning_n_candidates)},
{cmd:e(tuning_n_feasible)}, {cmd:e(tuning_n_skipped)},
{cmd:e(tuning_n_failed)} and {cmd:e(tuning_n_unique)} summarize the grid.
Candidate status is 0=infeasible or 1=ok; reasons are 0=ok, 1=k_gt_s,
2=jackknife_requires_s_lt_n, 3=s_gt_n.

{p2col 5 27 31 2:Scalars and metadata}{p_end}
{synopt:{cmd:e(N)}, {cmd:e(m)}, {cmd:e(kx)}, {cmd:e(p)}}sample, target and reporting dimensions{p_end}
{synopt:{cmd:e(k)}, {cmd:e(s)}, {cmd:e(lambda)}}resolved localization parameters{p_end}
{synopt:{cmd:e(neff_min)}, {cmd:e(neff_max)}}target ESS range; ess_min/ess_max are aliases{p_end}
{synopt:{cmd:e(neff)}, {cmd:e(ess)}}single-target ESS only{p_end}
{synopt:{cmd:e(clip)}, {cmd:e(denominator_tol)}, {cmd:e(one_sided)}}signal/status settings{p_end}
{synopt:{cmd:e(score_method)}}augmented_iv or paired_signals{p_end}
{synopt:{cmd:e(se_method)}}paired_rank_shift or none{p_end}
{synopt:{cmd:e(covariance)}}paired_frozen_score or none{p_end}
{synopt:{cmd:e(training)}, {cmd:e(overlap_role)}}declared provenance; instrument{p_end}
{synopt:{cmd:e(distance)}, {cmd:e(method)}}metric and rank-weight method{p_end}
{synopt:{cmd:e(tuning_mode)}, {cmd:e(tuning_source)}}default/grid and paper/s/lambda_grid/sgrid{p_end}

{pstd}
Input names are recorded in {cmd:e(xvars)}, {cmd:e(x_cols)},
{cmd:e(outcome)}, {cmd:e(signal)}, {cmd:e(denominator)}, {cmd:e(d_col)},
{cmd:e(a_col)}, {cmd:e(e_col)}, {cmd:e(my0_col)}, {cmd:e(my1_col)},
{cmd:e(md0_col)} and {cmd:e(md1_col)}. {cmd:e(treat)} refers to the
instrument in shared postestimation. {cmd:e(propensity_version)} records
clipped or not_used. {cmd:e(signature_vars)} and {cmd:e(ordered_signature)}
protect observation-level plots against changed values or row order.

{pstd}
Fit/replay return a 9-by-m {cmd:r(table)} with rows
{cmd:b se z pvalue ll ul df crit eform}, an m-by-6 {cmd:r(summary)} with
{cmd:estimate std_error z pvalue ci_lower ci_upper}, and {cmd:r(level)}.
Copy these temporary r-class results before another r-class command.
Undefined ratios retain target positions; unavailable bounds remain missing.

{pstd}
{help clate_dnn_diagnostics} supplies instrument-labelled target diagnostics,
first-stage SE/statistic, paired components and statuses. Shared stored
{cmd:e(diagnostics)} retains CATE slot names for compatibility; treated/control
slots mean instrument 1/0. Influence concentration uses
abs(w*(phi_Y-theta*phi_D)), not individual effects or deletion pseudo-values.
{help clate_dnn_plot} covers estimates, overlap, weights, influence and grids.
Current kernel/influence views use saved numeric observations and signals.
Legacy fits without snapshots require unchanged original sample, values and
row order. {help rdistnn_predict} relocalizes the saved frozen signals at new
targets; {help rdistnn_candidate} extracts reusable grid fits. Both preserve
active estimates by default. {help rdistnn_results} supplies graph-free tables.

{marker examples}{...}
{title:Examples}

{pstd}A self-contained oracle-prediction example; this illustrates mechanics:

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 80}{p_end}
{phang2}{cmd:. set seed 4711}{p_end}
{phang2}{cmd:. gen double x = rnormal()}{p_end}
{phang2}{cmd:. gen double e = invlogit(.3*x)}{p_end}
{phang2}{cmd:. gen byte A = runiform()<e}{p_end}
{phang2}{cmd:. gen double md1 = invlogit(.4+.2*x)}{p_end}
{phang2}{cmd:. gen byte D = A*(runiform()<md1)}{p_end}
{phang2}{cmd:. gen double my0 = .5*x}{p_end}
{phang2}{cmd:. gen double my1 = my0+md1}{p_end}
{phang2}{cmd:. gen double Y = my0+D+.3*rnormal()}{p_end}
{phang2}{cmd:. clate_dnn Y, xcols(x) x0(0) dcol(D) acol(A) ecol(e) my0col(my0) my1col(my1) md1col(md1) onesided training(oracle)}{p_end}
{phang2}{cmd:. matrix list e(components)}{p_end}
{phang2}{cmd:. matrix list e(status)}{p_end}
{phang2}{cmd:. clate_dnn_diagnostics, detail}{p_end}
{phang2}{cmd:. clate_dnn_plot, type(estimates) replace}{p_end}

{pstd}
The source archive includes {cmd:examples/clate.do}, with a
four-observation paired covariance calculation, an undefined-ratio example,
and equality checks between nuisance and supplied-signal routes.

{title:Also see}

{psee}{help rdistnn}, {help clate_dnn_diagnostics}, {help clate_dnn_plot}, {help cate_dnn}, {help weights_dnn}
