{smcl}
{* *! version 0.1.0  17oct2026}{...}
{title:Title}

{p 4 4 2}
{cmd:clate_dnn} {hline 2} Experimental conditional LATE with DNN localization

{title:Syntax}

{p 8 12 2}
{cmd:clate_dnn} {it:Y} {ifin}, {cmd:xcols(}{it:varlist}{cmd:)}
{cmd:x0(}{it:numlist}{cmd:)} | {cmd:x0matrix(}{it:matrix}{cmd:)}
{cmd:dcol(}{it:D}{cmd:)} {cmd:acol(}{it:A}{cmd:)} {cmd:ecol(}{it:e}{cmd:)}
{cmd:my0col(}{it:my0}{cmd:)} {cmd:my1col(}{it:my1}{cmd:)}
{cmd:md0col(}{it:md0}{cmd:)} {cmd:md1col(}{it:md1}{cmd:)}
[{cmd:onesided} {cmd:clip(}{it:real}{cmd:)} {it:options}]

{p 8 12 2}
{cmd:clate_dnn} {it:phi_y} {ifin}, {cmd:pseudooutcome}
{cmd:denominator(}{it:phi_d}{cmd:)} {cmd:xcols(}{it:varlist}{cmd:)}
{cmd:x0(}{it:numlist}{cmd:)} | {cmd:x0matrix(}{it:matrix}{cmd:)} [{it:options}]

{p 8 12 2}
{cmd:clate_dnn} [{cmd:, level(}{it:#}{cmd:)}]

{title:Description}

{pstd}
Constructs augmented reduced-form and first-stage signals from externally
supplied predictions. A is a binary instrument and D is a binary treatment.
Nuisance predictions condition on all adjustment covariates V=(X,Z); weights
use only the reporting covariates in {cmd:xcols()}. No learners are fitted.
The target averages both instrument-arm contrasts over the same conditional
distribution of Z given X, then takes their ratio.

{pstd}
For R=Y,D the signal is
phi_R=m_R1-m_R0+A*(R-m_R1)/e-(1-A)*(R-m_R0)/(1-e).
Arm indices refer to A, not D. The same clipped instrument propensity is used
for both signals. The estimator divides their common-weight local averages;
it never averages observation-level ratios. Signed first-stage signals and
finite negative denominators remain unchanged.

{title:Options}

{phang}
{cmd:xcols()} lists numeric reporting covariates; {cmd:x()} is an alias.
Supply exactly one of {cmd:x0()} (one coordinate per covariate) and
{cmd:x0matrix()} (targets in rows). Missing targets are rejected. Complete
observations in the requested if/in sample are used; missing outcome,
covariate, instrument, treatment or non-probability prediction values are
excluded. Probability predictions must be finite on the requested sample.

{phang}
{cmd:onesided} imposes D<=A. It permits omitted {cmd:md0col()}, interpreted as
zero, and requires a supplied md0 column to be zero. Without it all four
instrument-arm conditional means are required. D and A must be coded 0/1;
md0, md1 and raw e must lie in [0,1]. Treatment predictions are not clipped.

{phang}
{cmd:pseudooutcome denominator()} accepts two aligned finite signals instead
of nuisance predictions. It cannot be combined with nuisance options or
{cmd:onesided}. Instrument-overlap information is unavailable on this route.

{phang}
{cmd:clip()} defaults to 1e-8 and must lie strictly between zero and .5.
{cmd:denominatortol()} defaults to sqrt(machine epsilon). It flags a nonzero
local first stage whose absolute value divided by the maximum absolute input
first-stage signal is at most the tolerance. This flag changes no estimates
or covariance and supplies no statistical weak-instrument cutoff.

{phang}
{cmd:k()} defaults to 1; {cmd:s()} is the size of each subsample. The exact
top-k positive support has size n - s + k, or n - s + 1 under k(1).
Larger s tightens localization: it is opposite to the usual bandwidth convention.
The default s is floor(.4*n^((p+1)/(p+2))). {cmd:method(exact|asymp)} and
{cmd:distance(l2|mahalanobis)} share the CATE implementation. Mahalanobis
geometry is fixed at the full estimation sample throughout deletions.

{phang}
{cmd:tuningmode(grid)} uses the CATE paper grid (.3,.4,.5,.75), retaining the
ordinary .4 fit. {cmd:lambda()}, {cmd:lambdagrid()}, {cmd:sgrid()} and
{cmd:sdefault()} have the same controls, conflicts and feasibility metadata
as {help cate_dnn}. Grid comparisons do not perform automatic selection.
{cmd:nose} disables covariance and permits s=n; otherwise s<n is required.
{cmd:level()} controls pointwise normal intervals, including replay.

{phang}
{cmd:training(supplied|oracle|full_sample|oob|cross_fit)} records a declaration
only. Its default is supplied; prediction arrays do not verify that declaration.

{title:Inference contract}

{pstd}
Causal CLATE interpretation requires consistency, exclusion, conditional
instrument independence given V, monotonicity, instrument overlap and positive
local compliance. Without those assumptions the computation is a standardized
observed-data Wald ratio. These assumptions cannot be checked from predictions.

{pstd}
Orthogonality controls first-order population nuisance drift. Full-sample
inference additionally requires the affine theorem's design, oracle,
nuisance-accuracy, deletion-stability and undersmoothing conditions.
Supplying fitted predictions does not certify those conditions. The k(1)
exact-DNN theorem is not automatically extended to every supported k/method.

{pstd}
Covariance uses paired frozen-score rank-shift deletions with the same s,
predictions, targets and metric. It retains both components and their
cross-covariance, linearized at the full-sample ratio. Neither deleted
ratios nor refitted nuisance learners are used. {cmd:e(theta_vcov)} and
{cmd:e(V)} are already on the estimator scale: do not divide them by ESS.
Cross-target covariance does not make the displayed intervals simultaneous.
Wald intervals do not have weak-instrument-uniform validity, and jackknife
uncertainty does not remove localization bias.

{title:Stored results}

{pstd}
{cmd:e(estimates)} and {cmd:e(se_targets)} are full 1-by-m target tables,
including missing entries for undefined ratios or unavailable uncertainty.
{cmd:e(components)} is m-by-2 with columns reduced_form and first_stage.
{cmd:e(theta_vcov)} is m-by-m. {cmd:e(component_vcov)} is 2m-by-2m, ordered
y1,d1,y2,d2,..., retaining cross-target component covariance.

{pstd}
Stata does not permit missing coefficients in {cmd:ereturn post}.
{cmd:e(b)} therefore contains only defined ratios (theta for one target,
theta1, theta2,... for a batch); omitted target numbers remain identifiable.
{cmd:e(V)} is posted when uncertainty is available for every posted ratio.
If all ratios are undefined, the sample, components and canonical tables
remain available without coefficient results. No zero coefficient is
substituted for an undefined ratio.

{pstd}
{cmd:e(status)} has columns ratio_status (0=ok, 1=zero denominator,
2=numerical failure), negative_denominator, near_zero_denominator,
covariance_status (0=ok, 1=not computed, 2=zero score variance,
3=numerical failure), and wald_available. Zero score variance remains a
computed zero SE with no Wald test. Flags are computational diagnostics.

{pstd}
Targets, ESS, diagnostics, tuning candidates, estimates and SE tables use the
CATE conventions. Additional grid tables {cmd:e(tuning_reduced_form)},
{cmd:e(tuning_first_stage)} and {cmd:e(tuning_target_status)} retain paired
components and target-specific statuses. In the shared 25-column
{cmd:e(diagnostics)} schema, treated/control slots refer to instrument 1/0;
CLATE postestimation tables and plots label them as instrument arms.
Score concentration uses absolute w*(phi_Y-theta*phi_D), not individual effects.

{title:Example}

{pstd}See {cmd:examples/clate.do} in the source archive for a runnable oracle-nuisance example.

{title:Also see}

{psee}{help clate_dnn_diagnostics}, {help clate_dnn_plot}, {help cate_dnn}, {help weights_dnn}
