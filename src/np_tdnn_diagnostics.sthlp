{smcl}
{* *! version 0.1.0 18oct2026}{...}
{title:np_tdnn_diagnostics — Signed-weight diagnostics}

{p 4 4 2}
{cmd:np_tdnn_diagnostics} [{cmd:, targets(}{it:numlist}{cmd:)} {cmd:detail}
{cmd:saving(}{it:filename}{cmd:)} {cmd:replace}]

{p 4 4 2}
Requires active {help np_tdnn} estimates. Targets retain requested order.
Reports signed sum, positive mass, magnitude of negative mass, absolute mass,
inverse sum of squared weights ({cmd:inverse_l2_concentration}), largest absolute
weight, nonzero support count, and mean distance weighted by absolute mass.
Signed weights are not probabilities. No weighted outcome quantiles or
signed-weight outcome variances are reported. These diagnostics do not
establish bias removal or valid inference.

{p 4 4 2}
Returns {cmd:r(diagnostics)}, {cmd:r(compact)}, {cmd:r(n_rows)}, and
{cmd:r(selected)}. {cmd:saving()} writes the full diagnostic matrix as a
Stata dataset through the common diagnostic export helper.
