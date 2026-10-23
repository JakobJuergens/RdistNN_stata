{smcl}
{* *! version 0.1.1  2026-10-17}{...}
{viewerjumpto "Syntax" "weights_dnn##syntax"}{...}
{viewerjumpto "Description" "weights_dnn##description"}{...}
{viewerjumpto "Options" "weights_dnn##options"}{...}
{viewerjumpto "Numerical limits" "weights_dnn##numerical"}{...}
{viewerjumpto "Stored results" "weights_dnn##results"}{...}
{viewerjumpto "Examples" "weights_dnn##examples"}{...}

{title:Title}

{phang}
{bf:weights_dnn} {hline 2} Construct normalized DNN rank weights at one or more targets

{pstd}
Stata 17 or later. For the package overview and article collection, see
{help rdistnn}. This help describes package {bf:0.1.3}.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:weights_dnn} {it:xvars} {ifin}{cmd:,} {opt x0(numlist)} | {opt x0matrix(name)}
[{opt k(#)} {opt s(#)} {opt generate(newvar)}
{opt distance(string)} {opt scalerule(string)} {opt scalevec(numlist)}
{opt scaleeps(#)} {opt method(exact|asymp)}]

{marker description}{...}
{title:Description}

{pstd}
The command computes normalized rank-weight vectors around {cmd:x0()} or
the rows of {cmd:x0matrix()}.
Use numeric covariates in the order of the target coordinates. The usable
sample consists of observations selected by {cmd:if}/{cmd:in} with complete
covariates. {cmd:generate()} creates a double weight variable in original
observation order, missing outside that sample. Numeric weights, ranks and ordinary distances are always returned as
n-by-m matrices, with one usable observation per row and one target per column.
The command is r-class and does not replace the active estimation result.

{pstd}
The exact DNN construction averages the k nearest observations within every
size-s subset. It uses finite-population combination probabilities without
enumerating subsets. For global rank r, the coefficient is proportional to
sum over l=1,...,k of C(r-1,l-1) C(n-r,s-l), divided by k C(n,s), with
impossible choices zero. Positive support ends at rank n - s + k, or n - s + 1
under the default k=1. Thus k is not the number of positive full-sample weights,
and s is not a truncation to the nearest s observations.

{pstd}
The support formula describes the unaveraged rank table. Averaging assigned
weights across exact-distance groups can extend nonzero observation weights to
the end of a group crossing that edge. Use actual returned weights for realized
support and concentration; rank-only tables do not contain target-specific groups.

{pstd}
The asymptotic rule uses averaged Beta(l,s-l+1) CDF increments over adjacent
rank fractions and renormalizes. Lower or direct upper tails avoid subtracting
CDF values rounded to one. It is generally a full-rank approximation;
floating-point tail weights can be zero. Both methods normalize to sum one.
Ranks use represented selected-metric distance then original row index.
At each target, observations with exactly equal distance receive the mean of
their assigned rank weights, including distinct profiles and zero-weight ranks
beyond nominal support. No tolerance or rounding groups nearby distances.
Signed and unnormalized group mass is preserved; individual weights, support
and ESS can change. The paper assumes a continuous distance distribution, under
which old and new rules agree almost surely. Inference with persistent distance
ties and a uniform MSE improvement are not established.

{pstd}
The order s has direction opposite to the usual bandwidth convention.
Larger s concentrates exact support nearer the target for fixed n and k;
smaller s spreads weight and increases effective local sample size. ESS is
1/sum(w_i^2), a concentration measure rather than a support count or coverage
check. Fixed k>1 is implemented; the companion theory writes out k=1.
For exact DNN, averaging corresponds to uniform orderings within distance ties
and uniform subsets of size s. A tied boundary group shares only the remaining
top-k slots. With k=2, one closer observation and two tied boundary observations
in a subset receive weights (.5,.25,.25); including all three with equal weight
would define a different rule.

{marker numerical}{...}
{title:Numerical limits}

{pstd}
At subnormal scales, rounded host Beta endpoints can erase a tail mass or
produce different zero versus nonzero support in R and Stata. Averaging order
statistics and signed TDNN blending can also lose tiny masses. This remains
unresolved and can affect estimates and jackknife covariance for extreme
signals. Normalization or native deletion engines cannot restore probability
already lost during table construction. Exact weights use a different rule;
changing methods changes the calculation. See the
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/numerical-stability.md#asymptotic-tail-limits":numerical reference}.

{marker options}{...}
{title:Options}

{phang}
Specify exactly one of {opt x0(numlist)} or {opt x0matrix(name)}. Supply
one coordinate per positional covariate, in exactly that order; a target matrix
has one target per row. {cmd:generate()} requires a single target.

{phang}
{opt k(#)} is the positive integer nearest-neighbor count within each subset;
default 1. It must satisfy k<=s.

{phang}
{opt s(#)} is the integer subset/kernel order. Default 0 resolves
floor(.4*n^((p+1)/(p+2))), where p is the number of covariates and n is the
usable sample size. Supply a feasible explicit order if the default is below k.
The full endpoint s=n is valid because this command computes no jackknife.
There are no {cmd:lambda()}, tuning-grid or {cmd:nose} options here.

{phang}
{opt generate(newvar)} creates a new double variable containing observation
weights. Use a name not already in the dataset. Weights sum to one over the
usable sample; excluded observations have missing weights.

{phang}
{opt distance(string)} chooses l1, l2 (default), l1_scaled, l2_scaled, or the
Mahalanobis
full-sample covariance metric. Mahalanobis requires at least two covariates
and two complete observations for this command; singular covariance uses a
centered-design SVD whitening factor. Singular values greater than
max(n,p)*2^-52 times the largest singular value are retained; null directions
contribute zero. The covariance matrix and its inverse are not formed.
The estimator wrappers have their own dimension checks.
Euclidean neighborhoods depend on covariate units.

{pstd}
SVD fitting and projection rounding can change represented Mahalanobis distance
equality or ordering across packages or platforms, even for well-conditioned
full-rank designs far from the cutoff. Ranks, exact-distance groups, weights
and estimates can then differ. Group averaging uses each fit's computed
distances, with no tolerance or additional rounding, and cannot reconcile
these differences. See the
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/numerical-stability.md#fitted-mahalanobis-factorization":numerical geometry reference}.

{phang}
{opt method(exact|asymp)} chooses the finite-population exact rule (default)
or Beta-increment approximation. Changing it changes the weight vector.

{title:Native acceleration}

{pstd}
L1, L2, scaled L1 and scaled L2 distances use the installed platform-specific
plugin for substantial batches when its distance API is available, with Mata
fallback by default. Automatic dispatch keeps small jobs in Mata to avoid
matrix-transfer overhead; explicit plugin requirements bypass that policy.
The ordinary package bundles a Windows x86-64 binary;
Linux/macOS users can install matching release assets or build from source.
Scale fitting, stable ranks and weight tables use Mata; all five estimators can use
the compiled prefix/suffix jackknife, including under frozen Mahalanobis ranks.
Set global
{cmd:rdistnn_use_plugin} to {cmd:0} to force Mata or {cmd:1} to require a working
plugin for eligible calls in auto mode. Drop the global to restore default dispatch.
{cmd:global rdistnn_distance_engine plugin} requires compiled distances;
{cmd:mata} forces their reference path, and {cmd:auto} restores automatic dispatch.
This selector takes precedence over {cmd:rdistnn_use_plugin}.
{cmd:rdistnn_last_distance_engine} records the last actual distance route.
Mahalanobis factorization uses Mata; its projected L2 norms use the selected
engine.
For jackknife parity checks, {cmd:global rdistnn_jackknife_engine naive}
selects explicit deletions; {cmd:mata} selects fast Mata and {cmd:plugin}
requires compiled prefix/suffix sums. Drop that global for automatic compiled,
fast Mata, then naive fallback. Native distance and jackknife matrices created
through Mata may exceed {cmd:c(matsize)}; dispatch preserves that setting.

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
jackknife deletions and saved-model predictions. Ranks use represented selected-metric distance then original row index.
At each target, observations with exactly equal distance share the mean of their
assigned rank weights, including different covariate profiles and zero-weight
ranks beyond nominal support. No tolerance or rounding groups nearby distances.
Each deletion removes the actual observation before averaging the n-1 table
over each remaining group using its survivor count. Signed and unnormalized
mass is preserved; individual weights, support and ESS can change. Aligned row
permutations preserve estimates and covariance up to numerical precision.
The paper assumes a continuous distance distribution, under which old and new
rules agree almost surely. Inference with persistent distance ties and a uniform
MSE improvement are not established; outcomes and nuisances do not define groups.

{marker results}{...}
{title:Stored results}

{synoptset 18}{...}
{synopt:{cmd:r(neff)}}ESS, 1/sum(w_i^2){p_end}
{synopt:{cmd:r(n)}}usable sample size{p_end}
{synopt:{cmd:r(k)}}resolved nearest-neighbor count{p_end}
{synopt:{cmd:r(s)}}resolved order{p_end}
{synopt:{cmd:r(kx)}}number of covariates{p_end}

{pstd}
Copy r-class results before running another command that overwrites them.
{cmd:r(weights)}, {cmd:r(ranks)}, {cmd:r(distances)} and {cmd:r(groups)} retain
n-by-m numeric matrices; groups are integer IDs for exact represented distance
equality within each target. {cmd:r(row_ids)} maps rows to original observations.
{cmd:r(targets)} and {cmd:r(neff_targets)} retain target coordinates and ESS.
{cmd:r(geometry)} retains fitted scales or packed Mahalanobis whitening
geometry when used.
{cmd:r(m)} is the target count; scalar {cmd:r(neff)} refers to the first target. To reconstruct an estimator's weighted signal,
match its complete-case sample using a copy of {cmd:e(sample)}, as well as
its target, k, s, method and distance. Covariate-only weights may otherwise
use rows dropped for a missing outcome/nuisance and produce different ranks.

{marker examples}{...}
{title:Examples}

{pstd}These examples assume numeric x1 and x2 in the current dataset.{p_end}
{phang2}{cmd:. weights_dnn x1 x2, x0(0 0) s(8) generate(w)}{p_end}
{phang2}{cmd:. scalar weight_ess = r(neff)}{p_end}
{phang2}{cmd:. summarize w}{p_end}
{phang2}{cmd:. weights_dnn x1 x2, x0(0 0) s(8) method(asymp) generate(wa)}{p_end}
{phang2}{cmd:. weights_dnn x1 x2, x0(0 0) s(8) distance(mahalanobis) generate(wm)}{p_end}

{pstd}
A matched regression reconstruction starts with a fit and a saved sample:
{p_end}
{phang2}{cmd:. np_dnn Y, xcols(x1 x2) x0(0 0) s(8)}{p_end}
{phang2}{cmd:. generate byte used = e(sample)}{p_end}
{phang2}{cmd:. weights_dnn x1 x2 if used, x0(0 0) s(8) generate(wfit)}{p_end}
{phang2}{cmd:. generate double wy = wfit*Y}{p_end}
{phang2}{cmd:. summarize wy, meanonly}{p_end}
{pstd}The sum r(sum) equals the fitted NP mean up to numerical precision.{p_end}

{title:Further reading}

{pstd}
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/weights-and-tuning.md":Weights and tuning article}
and its standalone do-file cover exact/asymptotic comparison, matched samples
and candidate grids. Custom distance registries are available in the R interface.

{pstd}Also see {help cate_dnn}, {help clate_dnn}, {help np_dnn}, {help plqr_dnn}.{p_end}
