{smcl}
{* *! version 0.1.0  2026-10-18}{...}
{title:rdistnn_predict — Evaluate a saved fit at new targets}
{p 8 12 2}{cmd:rdistnn_predict,} {cmd:x0(}{it:numlist}{cmd:)} | {cmd:x0matrix(}{it:name}{cmd:)}
[{cmd:sefit interval(confidence) level(}{it:#}{cmd:) post store(}{it:name}{cmd:)}]{p_end}
{pstd}
Specify exactly one target option. x0matrix() has one target per row and one
column per fitted covariate, in e(xvars) order. The saved training sample,
frozen signals, k, s, weight method and metric are reused. New geometry is
computed using that fixed training sample, including its Mahalanobis metric;
no nuisance model is fitted. No original variables are needed in active data.
{pstd}
At every new target, exactly identical raw training-covariate rows share their
block's average weight. With sefit, each deletion removes the actual member and
divides the n-1 table's remaining block mass among the survivors. The command
warns that this permutation-invariant modification is outside the paper's
theory. Distinct vectors at equal distance remain separate; see {help weights_dnn}.
{pstd}
By default uncertainty is not computed. sefit or interval(confidence) requests
the same frozen-score delete-one covariance as estimation; this requires s<n.
level() controls pointwise normal limits, in percent. This remains localized
target inference and does not remove localization bias or validate nuisances.
{pstd}
The command preserves active data and e() by default. post explicitly replaces
the active fit with a reusable derived fit. store(name) saves a derived fit
without replacing the active fit unless post is also given. Derived fits retain
the numeric model but omit the original tuning grid. They can be saved using
estimates save and reused after estimates use, even with an empty active dataset.
Their e(sample) identifies original rows only when those inputs and their order
can still be verified; otherwise it is empty. Use saved observation IDs for reports.
For np_tdnn, posted coefficients retain target1, target2, and so on (target1
for a single target), preserving coefficient names used by lincom.
{title:Returned results}
{pstd}
Returns r(fit), r(se_fit), r(vcov), r(table), r(summary), r(index), r(data),
r(level), r(s), r(quantity), and r(stored) when requested; see {help rdistnn_results}
for shapes. CLATE covariance is ratio-scale; components and status remain in a
posted or stored derived fit. Undefined ratios are preserved.
PLQR returns evaluated moments, with r(theta_eval) separately retaining supplied
coefficients. It performs no root solving, coefficient inference or density fitting.
{pstd}
Old saved fits without e(model_version)=2 must be refitted once. Current snapshots
are validated for consistent controls, dimensions, row identities, fitted
geometry and cached ordering before prediction. Invalid snapshots request a
refit while preserving the active data and e(). Model snapshots
store fitted observations and signals in .ster files; storage grows with sample
size and the number of targets. These files contain the fitted numeric inputs.
{title:See also}
{pstd}{help rdistnn}, {help rdistnn_results}, {help rdistnn_candidate}{p_end}
