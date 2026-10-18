{smcl}
{* *! version 0.1.1  2026-10-17}{...}
{viewerjumpto "Syntax" "weights_dnn##syntax"}{...}
{viewerjumpto "Description" "weights_dnn##description"}{...}
{viewerjumpto "Options" "weights_dnn##options"}{...}
{viewerjumpto "Stored results" "weights_dnn##results"}{...}
{viewerjumpto "Examples" "weights_dnn##examples"}{...}

{title:Title}

{phang}
{bf:weights_dnn} {hline 2} Construct normalized DNN rank weights at one target

{pstd}
Stata 17 or later. For the package overview and article collection, see
{help rdistnn}. This help describes the development branch.

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:weights_dnn} {it:xvars} {ifin}{cmd:,} {opt x0(numlist)}
[{opt k(#)} {opt s(#)} {opt generate(newvar)}
{opt distance(l2|mahalanobis)} {opt method(exact|asymp)}]

{marker description}{...}
{title:Description}

{pstd}
The command computes one normalized rank-weight vector around {cmd:x0()}.
Use numeric covariates in the order of the target coordinates. The usable
sample consists of observations selected by {cmd:if}/{cmd:in} with complete
covariates. {cmd:generate()} creates a double weight variable in original
observation order, missing outside that sample. Without it, weights are
computed for the reported effective sample size but no vector is returned.
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
The asymptotic rule uses averaged Beta(l,s-l+1) CDF increments over adjacent
rank fractions and renormalizes. It is generally a full-rank approximation;
floating-point tail weights can be zero. Both methods normalize to sum one.
Distance ties are broken by observation order. Match row order when comparing
with R or reconstructing a fit.

{pstd}
The order s has direction opposite to the usual bandwidth convention.
Larger s concentrates exact support nearer the target for fixed n and k;
smaller s spreads weight and increases effective local sample size. ESS is
1/sum(w_i^2), a concentration measure rather than a support count or coverage
check. Fixed k>1 is implemented; the companion theory writes out k=1.

{marker options}{...}
{title:Options}

{phang}
{opt x0(numlist)} is required. Supply one coordinate per positional covariate,
in exactly that order. This command has no {cmd:x0matrix()} interface.

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
{opt distance(l2|mahalanobis)} chooses Euclidean distance (default) or the
full-sample covariance metric. Mahalanobis requires at least two covariates
and two complete observations for this command; singular covariance uses a
generalized inverse. The estimator wrappers have their own dimension checks.
Euclidean neighborhoods depend on covariate units.

{phang}
{opt method(exact|asymp)} chooses the finite-population exact rule (default)
or Beta-increment approximation. Changing it changes the weight vector.

{title:Native acceleration}

{pstd}
Eligible L2 calls try the installed platform-specific plugin, with Mata fallback
on load/execution failure. The ordinary package bundles a Windows x86-64 binary;
Linux/macOS users can install matching release assets or build from source.
The estimator wrappers and Mahalanobis path use Mata. Set global
{cmd:rdistnn_use_plugin} to {cmd:0} to force Mata or {cmd:1} to require a working
plugin for eligible L2 calls. Drop the global to restore default dispatch.
Requiring the plugin does not add Mahalanobis plugin support.

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
No W, rank, order or distance matrix is returned; retain {cmd:generate()}
output for downstream work. To reconstruct an estimator's weighted signal,
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
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/feature/clate/docs/weights-and-tuning.md":Weights and tuning article}
and its standalone do-file cover exact/asymptotic comparison, matched samples
and candidate grids. R's batch weight objects, custom registries and TDNN
weights have no Stata public equivalent.

{pstd}Also see {help cate_dnn}, {help clate_dnn}, {help np_dnn}, {help plqr_dnn}.{p_end}
