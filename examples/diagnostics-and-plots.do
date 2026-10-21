// Inspect local support, plot fits and copy numeric reports.
// Install rdistnn first, or run dev-setup.do from the repository root.
// This self-contained example clears the active dataset.
version 17
set more off
clear
// All nuisance values are oracle quantities; intervals are pointwise.
// Remove nodraw to display graphs while stepping through the example.
set seed 20261017
set obs 72
gen double x1 = 2*runiform()-1
gen double x2 = rnormal()
gen double e_oracle = invlogit(.3*x1-.2*x2)
gen byte D = runiform()<e_oracle
gen double m0_oracle = .4*x1-.2*x2
gen double m1_oracle = m0_oracle+.8+.3*sin(x1)
gen double Y = m0_oracle+D*(m1_oracle-m0_oracle)+.3*rnormal()
matrix X0 = (-.8,0 \ -.4,0 \ 0,0 \ .4,0 \ .8,0)

// Candidate 2 is retained as the default; the grid is a sensitivity comparison.
cate_dnn Y, xcols(x1 x2) x0matrix(X0) tcol(D) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) ///
    tuningmode(grid) sgrid(4 8 12) sdefault(8)
estimates store docs_cate
cate_dnn, level(90)
matrix coefficient_report = r(summary)
matrix list coefficient_report

// Select two profiles and save their diagnostics to a temporary dataset.
tempfile cate_report
cate_dnn_diagnostics, targets(5 1) saving("`cate_report'") replace
matrix localization = r(localization)
preserve
use "`cate_report'", clear
list target_id ess treated_mass control_mass top5_score_share, noobs
restore

// A profile, local overlap, kernel geometry and score contributions answer
// different questions. Copy r(plot_data) immediately when making custom reports.
cate_dnn_plot, style(profile) name(docs_cate_profile) replace nodraw
cate_dnn_plot, targets(2 4) level(90) name(docs_cate_est) replace nodraw
matrix plotted_estimates = r(plot_data)
matrix list plotted_estimates
cate_dnn_plot, type(overlap) view(mass) name(docs_cate_mass) replace nodraw
cate_dnn_plot, type(weights) view(kernel) targets(3) ///
    name(docs_cate_kernel) replace nodraw
cate_dnn_plot, type(influence) targets(3) nlabel(0) ///
    name(docs_cate_influence) replace nodraw
cate_dnn_plot, type(tuning) allcandidates name(docs_cate_tune) replace nodraw
// Use the standalone default candidate for a readable four-panel overview;
// inspect the grid in the separate tuning graph above.
rdistnn_candidate, candidate(2) store(docs_dashboard_fit)
estimates restore docs_dashboard_fit
cate_dnn_plot, type(dashboard) name(docs_cate_dashboard) replace nodraw
estimates restore docs_cate

// NP outcome quantiles describe local outcome dispersion, not uncertainty
// of the conditional mean. Compare the two panels.
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid) sgrid(4 8 12) sdefault(8)
estimates store docs_np
np_dnn_plot, name(docs_np_est) replace nodraw
np_dnn_plot, type(outcome) name(docs_np_outcome) replace nodraw

// Generic paired signals carry no instrument data for overlap panels.
// These simulated signals illustrate ratio reporting, not causal identification.
gen double phi_d = .7+.1*x1+.04*rnormal()
gen double phi_y = .9*phi_d+.1*x1+.2*rnormal()
clate_dnn phi_y, denominator(phi_d) pseudooutcome ///
    xcols(x1 x2) x0matrix(X0) tuningmode(grid) sgrid(4 8 12) sdefault(8)
estimates store docs_clate
clate_dnn_diagnostics, targets(1 3 5)
clate_dnn_plot, style(forest) name(docs_clate_est) replace nodraw
clate_dnn_plot, type(dashboard) name(docs_clate_dashboard) replace nodraw

// Optional persistent exports, relative to the current directory.
local export_figures 0
if `export_figures' {
    capture mkdir "build"
    capture mkdir "build/documentation"
    estimates restore docs_cate
    cate_dnn_plot, style(profile) name(docs_export) replace nodraw ///
        saving("build/documentation/cate-profile.svg")
    estimates restore docs_np
    np_dnn_plot, type(outcome) name(docs_export) replace nodraw ///
        saving("build/documentation/np-outcome.pdf")
}
