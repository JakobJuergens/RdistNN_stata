// Reporting, diagnostics, graph data and original-data guards.
// Install the package from master, or run dev-setup.do from the repository root.
// All inputs are generated here; this example clears the active dataset.
// Oracle CATE nuisances illustrate mechanics, not fitted-learner performance.
// Paired CLATE signals illustrate ratio reporting, not causal identification.
version 17
set more off
clear
set seed 20261017
set obs 72
gen long row_id = _n
gen double x1 = 2*runiform()-1
gen double x2 = rnormal()
gen double e_oracle = invlogit(.3*x1-.2*x2)
gen byte D = runiform()<e_oracle
gen double m0_oracle = .4*x1-.2*x2
gen double m1_oracle = m0_oracle+.8+.3*sin(x1)
gen double Y = m0_oracle+D*(m1_oracle-m0_oracle)+.3*rnormal()
matrix X0 = (-.8,0 \ -.4,0 \ 0,0 \ .4,0 \ .8,0)
local original_N = _N

// An explicit small grid makes candidate IDs easy to follow: candidate 2 is fitted.
quietly cate_dnn Y, xcols(x1 x2) x0matrix(X0) tcol(D) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) ///
    tuningmode(grid) sgrid(4 8 12) sdefault(8)
estimates store docs_cate
cate_dnn, level(90)
matrix CT = r(table)
matrix CS = r(summary)
assert rowsof(CT)==9 & colsof(CT)==5
assert rowsof(CS)==5 & colsof(CS)==6
assert abs(CS[1,5]-(CS[1,1]-invnormal(.95)*CS[1,2]))<1e-12

// Diagnostic selection preserves the requested target order and scientific fields.
tempfile cate_report
cate_dnn_diagnostics, targets(5 1) saving("`cate_report'") replace
assert r(n_rows)==2
assert "`r(selected)'"=="5 1"
matrix CD = r(diagnostics)
matrix CL = r(localization)
assert CD[1,1]==5 & CD[2,1]==1
assert colsof(CL)==14
assert CD[1,colnumb(CD,"overlap_available")]==1
assert CD[1,colnumb(CD,"arm_ess_available")]==1
assert CD[1,colnumb(CD,"ess")]>0
assert CD[1,colnumb(CD,"physical_radius_90")]>=0
preserve
use "`cate_report'", clear
assert _N==2
assert target_id[1]==5 & target_id[2]==1
assert treated_mass+control_mass>1-1e-12
confirm variable top5_score_share
restore
assert _N==`original_N'
cate_dnn_diagnostics, candidates(3 1) detail
matrix CG = r(diagnostics)
assert CG[1,1]==3 & CG[2,1]==1
assert colsof(CG)==13

// Every CATE plot type; nodraw supports unattended execution.
cate_dnn_plot, targets(2 4) level(90) name(docs_cate_est) replace nodraw
matrix CP = r(plot_data)
assert r(n_rows)==2
assert CP[1,colnumb(CP,"target_id")]==2
assert abs(CP[1,colnumb(CP,"ci_lower")] - ///
    (CP[1,colnumb(CP,"estimate")]-invnormal(.95)*CP[1,colnumb(CP,"se")]))<1e-12
cate_dnn_plot, style(profile) name(docs_cate_profile) replace nodraw
cate_dnn_plot, style(forest) targets(1 3 5) candidates(1 3) ///
    name(docs_cate_compare) replace nodraw
assert r(n_rows)==6
cate_dnn_plot, style(pointrange) targets(1 5) name(docs_cate_points) replace nodraw
cate_dnn_plot, type(tuning) targets(2 4) candidates(1 3) ///
    name(docs_cate_tune) replace nodraw
assert r(n_rows)==4 & r(n_plotted)==2
cate_dnn_plot, type(diagnostics) targets(1 5) name(docs_cate_diag) replace nodraw
cate_dnn_plot, type(diagnostics) allcandidates name(docs_cate_griddiag) replace nodraw
cate_dnn_plot, type(weights) view(profile) allcandidates logweights ///
    name(docs_cate_rank) replace nodraw
matrix CW = r(plot_data)
mata: assert(all(st_matrix("CW")[,2]:>0))
cate_dnn_plot, type(weights) view(cumulative) name(docs_cate_cum) replace nodraw
cate_dnn_plot, type(weights) view(kernel) targets(2 4) ///
    name(docs_cate_kernel) replace nodraw
assert r(n_rows)==2*`original_N' & r(n_plotted)==2
cate_dnn_plot, type(overlap) view(mass) sharelimits(.1 .9) ///
    name(docs_cate_mass) replace nodraw
cate_dnn_plot, type(overlap) view(armess) targets(1 5) ///
    name(docs_cate_armess) replace nodraw
cate_dnn_plot, type(overlap) view(propensity) name(docs_cate_prop) replace nodraw
cate_dnn_plot, type(influence) targets(3) nlabel(0) ///
    name(docs_cate_influence) replace nodraw
assert r(n_rows)==`original_N' & r(n_plotted)==1
cate_dnn_plot, type(dashboard) name(docs_cate_dashboard) replace nodraw
assert missing(r(n_rows)) & r(n_panels)==5
matrix CPI = r(panel_info)
assert rowsof(CPI)==5 & colsof(CPI)==3

// Invalid view options and changed values/order must fail with useful contracts.
capture noisily cate_dnn_plot, type(weights) targets(1) replace nodraw
assert _rc==198
scalar docs_x1_first = x1[1]
quietly replace x1 = x1+.01 in 1
capture noisily cate_dnn_plot, type(weights) view(kernel) replace nodraw
assert _rc==498
// Stored-result panels remain available after values change.
cate_dnn_plot, targets(1) name(docs_stored_est) replace nodraw
quietly replace x1 = scalar(docs_x1_first) in 1
gsort -row_id
capture noisily cate_dnn_plot, type(influence) targets(1) replace nodraw
assert _rc==498
sort row_id
estimates restore docs_cate
cate_dnn_plot, type(weights) view(kernel) targets(1) ///
    name(docs_restored_kernel) replace nodraw
assert r(n_rows)==`original_N'

// NP distinguishes outcome variation from uncertainty of the mean.
quietly np_dnn Y, xcols(x1 x2) x0matrix(X0) ///
    tuningmode(grid) sgrid(4 8 12) sdefault(8)
estimates store docs_np
matrix NB = e(estimates)
np_dnn_diagnostics, targets(1 3 5) detail
matrix ND = r(diagnostics)
assert abs(ND[1,colnumb(ND,"weighted_mean_y")]-NB[1,1])<1e-12
assert ND[1,colnumb(ND,"weighted_y_q10")]<=ND[1,colnumb(ND,"weighted_y_q50")]
assert ND[1,colnumb(ND,"weighted_y_q50")]<=ND[1,colnumb(ND,"weighted_y_q90")]
np_dnn_plot, level(90) name(docs_np_est) replace nodraw
np_dnn_plot, type(outcome) name(docs_np_outcome) replace nodraw
matrix NO = r(plot_data)
assert colsof(NO)==8 & r(n_rows)==5
assert NO[1,colnumb(NO,"weighted_sd_y")]>=0
np_dnn_plot, type(diagnostics) name(docs_np_diag) replace nodraw
np_dnn_plot, type(tuning) allcandidates name(docs_np_tune) replace nodraw
np_dnn_plot, type(weights) view(profile) name(docs_np_weights) replace nodraw
np_dnn_plot, type(weights) view(kernel) targets(3) name(docs_np_kernel) replace nodraw
np_dnn_plot, type(dashboard) name(docs_np_dashboard) replace nodraw
assert r(n_panels)==5
quietly np_dnn Y, xcols(x1 x2) x0(0 0) s(8) nose
np_dnn_plot, name(docs_np_nose) replace nodraw
matrix NN = r(plot_data)
assert missing(NN[1,colnumb(NN,"se")])
assert missing(NN[1,colnumb(NN,"ci_lower")])
np_dnn_plot, type(outcome) name(docs_np_nose_outcome) replace nodraw

// Paired signals need no instrument columns and therefore have no overlap panel.
gen double phi_d = .7+.1*x1+.04*rnormal()
gen double phi_y = .9*phi_d+.1*x1+.2*rnormal()
quietly clate_dnn phi_y, denominator(phi_d) pseudooutcome ///
    xcols(x1 x2) x0matrix(X0) tuningmode(grid) sgrid(4 8 12) sdefault(8)
estimates store docs_clate
tempfile clate_report
clate_dnn_diagnostics, targets(5 1) saving("`clate_report'") replace
matrix LD = r(diagnostics)
matrix LC = r(components)
matrix LS = r(status)
assert LC[1,1]==5 & LC[2,1]==1
assert abs(LC[1,9]-LC[1,2]/LC[1,3])<1e-12
assert colsof(LS)==5
assert LD[1,colnumb(LD,"overlap_available")]==0
preserve
use "`clate_report'", clear
assert target_id[1]==5 & target_id[2]==1
confirm variable component_covariance
confirm variable wald_available
confirm variable instrument_1_mass
restore
clate_dnn_plot, style(forest) name(docs_clate_est) replace nodraw
clate_dnn_plot, type(tuning) allcandidates name(docs_clate_tune) replace nodraw
clate_dnn_plot, type(diagnostics) name(docs_clate_diag) replace nodraw
clate_dnn_plot, type(weights) view(kernel) targets(3) name(docs_clate_kernel) replace nodraw
clate_dnn_plot, type(influence) targets(3) nlabel(0) name(docs_clate_influence) replace nodraw
matrix LI = r(plot_data)
assert colsof(LI)==9 & r(n_rows)==`original_N'
assert abs(LI[1,7]-LI[1,5]*LI[1,6])<1e-12
capture noisily clate_dnn_plot, type(overlap) replace nodraw
assert _rc==498
clate_dnn_plot, type(dashboard) name(docs_clate_dashboard) replace nodraw
assert r(n_panels)==4
local docs_panels "`r(panels)'"
local docs_panels : list retokenize docs_panels
assert "`docs_panels'"=="estimates diagnostics weights_cumulative tuning"

// Optional exports: run from the repository root and set this local to 1.
// The only persistent outputs are below the ignored build/ directory.
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
    estimates restore docs_clate
    clate_dnn_plot, type(dashboard) name(docs_export) replace nodraw ///
        saving("build/documentation/clate-dashboard.gph")
}
assert _N==`original_N'
display as text "Reporting and graph-data example completed."
