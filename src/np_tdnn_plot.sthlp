{smcl}
{* *! version 0.1.0 18oct2026}{...}
{title:np_tdnn_plot — Two-scale regression plots}

{p 4 4 2}
{cmd:np_tdnn_plot} [{cmd:, type(estimates|weights|diagnostics)}
{cmd:targets(}{it:numlist}{cmd:)} {cmd:style(profile|forest|pointrange)}
{cmd:xindex(}{it:#}{cmd:)} {cmd:level(}{it:#}{cmd:)}
{cmd:name(}{it:name}{cmd:)} {cmd:replace} {cmd:nodraw}
{cmd:saving(}{it:filename}{cmd:)}]

{p 4 4 2}
Requires active {help np_tdnn} estimates. Default {cmd:type(estimates)} uses
the common estimate plot adapter and pointwise confidence bounds.
{cmd:style()} and {cmd:xindex()} apply to estimates only.
{cmd:type(weights)} displays final signed weights around zero, with selected
{cmd:targets()} in separate panels. Duplicate-block averaging can make profiles
differ by target; {cmd:r(plot_data)} includes rank, weight and target_id.
{cmd:type(diagnostics)} plots total absolute weight by selected target.
Signed weights do not define probability distributions; no outcome-quantile
or cumulative-probability panels are provided.

{p 4 4 2}
Returns graph source data in {cmd:r(plot_data)}, graph name in {cmd:r(name)},
and {cmd:r(n_rows)}. Export uses the common graph helper and supports .gph,
.png, .pdf, .svg, .tif, .jpg and .eps. Numeric model state allows these plots
after restoring estimates without the original dataset.
