// RdistNN CATE and regression: both score routes, batches and reporting.
// Install the package from master, or run dev-setup.do from the repository root.
// The script clears the active dataset; it writes no files.
// Nuisances below are oracle values, not predictions fitted on these data.
version 17
set more off
clear
set obs 180
set seed 20261017
gen double x1 = rnormal()
gen double x2 = .4*x1 + rnormal()
gen double e_oracle = invlogit(.3*x1 - .2*x2)
gen byte D = runiform() < e_oracle
gen double tau_oracle = .75 + .3*x1
gen double m0_oracle = 1 + x1 - .5*x2
gen double m1_oracle = m0_oracle + tau_oracle
gen double m_oracle = m0_oracle + e_oracle*tau_oracle
gen double Y = m0_oracle + D*tau_oracle + .3*rnormal()
matrix X0 = (0,0 \ .5,0 \ -.5,0)

// AIPW: the main fit is anchored to the default candidate in grid mode.
// Set this report to 95 explicitly so the later 90-percent comparison also
// works when the caller has changed Stata's configured confidence level.
cate_dnn Y, xcols(x1 x2) x0matrix(X0) tcol(D) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) tuningmode(grid) level(95)
matrix AIPW_summary = r(summary)
matrix AIPW_table = r(table)
matrix AIPW_b = e(b)
matrix AIPW_V = e(V)
scalar AIPW_s = e(s)
assert e(N) == 180 & e(m) == 3
assert rowsof(AIPW_summary) == 3 & colsof(AIPW_summary) == 6
assert rowsof(AIPW_table) == 9 & colsof(AIPW_table) == 3
forvalues j = 1/3 {
    assert AIPW_summary[`j',2] >= 0 & AIPW_summary[`j',2] < .
    assert AIPW_summary[`j',5] <= AIPW_summary[`j',1]
    assert AIPW_summary[`j',1] <= AIPW_summary[`j',6]
}
estimates store cate_aipw
matrix list AIPW_summary
cate_dnn_diagnostics, targets(1 3) detail
matrix CATE_diagnostics = r(diagnostics)
assert rowsof(CATE_diagnostics) == 2

// A supplied AIPW signal represents the same fit on the same sample.
gen double eclip = min(max(e_oracle,1e-8),1-1e-8)
gen double rho_aipw = m1_oracle-m0_oracle + D*(Y-m1_oracle)/eclip ///
    - (1-D)*(Y-m0_oracle)/(1-eclip)
cate_dnn rho_aipw, xcols(x1 x2) x0matrix(X0) pseudooutcome
matrix supplied_b = e(b)
matrix supplied_V = e(V)
mata: st_numscalar("route_difference", max(abs(st_matrix("AIPW_b")-st_matrix("supplied_b"))))
assert route_difference < 1e-10
mata: st_numscalar("route_v_difference", max(abs(st_matrix("AIPW_V")-st_matrix("supplied_V"))))
assert route_v_difference < 1e-10
assert e(s) == AIPW_s
assert "`e(score_method)'" == "pseudo_outcome"


// Replay changes the pointwise interval level without changing estimates.
estimates restore cate_aipw
cate_dnn, level(90)
matrix AIPW_summary90 = r(summary)
assert r(level) == 90
assert AIPW_summary90[1,1] == AIPW_b[1,1]
assert AIPW_summary90[1,6]-AIPW_summary90[1,5] < ///
    AIPW_summary[1,6]-AIPW_summary[1,5]

// The batch covariance supports a prespecified first-minus-second contrast.
matrix L = (1,-1,0)
matrix CATE_difference = AIPW_b*L'
matrix CATE_difference_variance = L*AIPW_V*L'
assert CATE_difference_variance[1,1] >= 0
matrix list CATE_difference
matrix list CATE_difference_variance

// Estimate-only inference is explicitly missing, rather than a zero SE.
cate_dnn rho_aipw, xcols(x1 x2) x0(0 0) pseudooutcome nose
matrix estimate_only = r(summary)
assert "`e(se_method)'" == "none"
assert estimate_only[1,1] < . & missing(estimate_only[1,2])
local stored_matrices : e(matrices)
assert strpos(" `stored_matrices' "," V ") == 0

// NP regression estimates E[Y|X], rather than a treatment effect.
// One artificial missing outcome illustrates complete-case handling after if.
gen double Y_np = Y
replace Y_np = . in 180
np_dnn Y_np if _n > 5, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix NP_summary = r(summary)
matrix NP_b = e(b)
matrix NP_V = e(V)
assert e(N) == 174 & e(m) == 3
assert e(sample) == (_n > 5 & !missing(Y_np,x1,x2))
forvalues j = 1/3 {
    assert NP_summary[`j',1] < . & NP_summary[`j',2] >= 0 & NP_summary[`j',2] < .
}
estimates store np_fit
np_dnn_diagnostics, detail
matrix NP_diagnostics = r(diagnostics)
np_dnn_diagnostics, allcandidates
matrix NP_candidates = r(diagnostics)
assert rowsof(NP_candidates) == 4

// Restore the fit matching the postestimation command you want to use.
estimates restore cate_aipw
cate_dnn
matrix list e(targets)
matrix list e(se_targets)
estimates restore np_fit
np_dnn
matrix list NP_summary
display as result "CATE routes, batch reports and NP sample checks passed."
