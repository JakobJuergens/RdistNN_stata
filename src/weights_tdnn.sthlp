{smcl}
{* *! version 0.1.0 18oct2026}{...}
{title:weights_tdnn — Signed two-scale rank weights}

{p 4 4 2}
{cmd:weights_tdnn} {it:varlist} {ifin},
{cmd:x0(}{it:numlist}{cmd:)} [{cmd:s1(}{it:#}{cmd:)} {cmd:s2(}{it:#}{cmd:)}
{cmd:dim(}{it:#}{cmd:)} {cmd:method(exact|asymp)} {cmd:nonormalize}
{cmd:distance(}{it:metric}{cmd:)} {cmd:scalerule(}{it:rule}{cmd:)}
{cmd:scalevec(}{it:numlist}{cmd:)} {cmd:scaleeps(}{it:#}{cmd:)}
{cmd:generate(}{it:newvar}{cmd:)}]

{p 4 4 2}
Use {cmd:x0matrix(}{it:M}{cmd:)} for batch targets. One target is required
for {cmd:generate()}. Defaults: s1=1, s2=2, dim=number of covariates.
Require integer 1 <= s1 < s2 <= n and positive integer dim.
Distances and scaling options match {help np_tdnn}.

{p 4 4 2}
The command matches R {cmd:rank_weight_tdnn_exact()} and
{cmd:rank_weight_tdnn_asymp()}. Exact support ends at n-s1+1.
Asymptotic tables use Beta increments with normalization enabled by default.
Weights may be negative and must not be interpreted as probabilities.

{pstd}
The support formula describes the unaveraged rank table. Averaging assigned
weights across exact-distance groups can extend nonzero observation weights to
the end of a group crossing that edge. Use actual returned weights for realized
support and concentration; rank-only tables do not contain target-specific groups.

{p 4 4 2}
At each target, observations with exactly equal represented selected-metric
distance receive the mean of their assigned signed rank weights, without
absolute values or renormalization. No tolerance groups nearby distances.
Group mass is preserved; individual weights, support and concentration can
change. Rank tables alone cannot reconstruct the final weights without the
target-specific groups. Inference with persistent distance ties is not established.

{marker numerical}{...}
{title:Numerical limits}

{pstd}
Asymptotic host Beta endpoints, order-statistic averaging and signed TDNN
blending can lose tiny masses or produce R/Stata support differences at
subnormal scales. This remains unresolved. Extreme signals can amplify the
effect in estimates and jackknife covariance. Normalization and native deletion
engines do not restore probability lost in table construction. See
{help weights_dnn##numerical:weights_dnn} and the
{browse "https://github.com/JakobJuergens/RdistNN_stata/blob/master/docs/numerical-stability.md#asymptotic-tail-limits":numerical reference}.

{title:Saved results}

{pstd}
r(n_tied_targets) and r(n_weak_targets) report target counts even under quietly.
The weak flag uses signed inverse squared-weight concentration, not probability
ESS. These descriptive notices do not change weights or the active e() results.

{p 4 4 2}
{cmd:r(weights)}, {cmd:r(ranks)}, {cmd:r(distances)} and {cmd:r(groups)} are
n-by-m matrices in estimation-sample observation order. Groups are integer IDs
for exact represented distance equality within each target. {cmd:r(row_ids)} maps their rows to
the original observation numbers, including after {cmd:if}, {cmd:in}, and
missing-covariate exclusions. {cmd:r(rank_weights)} is the n-by-1
rank table. {cmd:r(targets)} retains target rows; {cmd:r(geometry)} retains
fitted geometry when needed. {cmd:r(diagnostics)} reports signed sum,
positive mass, negative mass magnitude, absolute mass, inverse squared-weight
concentration, largest absolute weight, nonzero support count, and mean
distance weighted by absolute mass. Concentration is a numerical diagnostic,
not a probability-based effective sample size. Scalars {cmd:r(n)},
{cmd:r(m)}, {cmd:r(s1)}, {cmd:r(s2)}, {cmd:r(dim)}, {cmd:r(normalize)}
record controls.
