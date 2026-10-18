// Standalone exact/asymptotic weights and integrated-order sensitivity.
// Clears the active dataset. Install from master or run dev-setup.do first.
version 17
clear
set more off
set obs 60
set seed 8731
gen double x1 = rnormal()
gen double x2 = .5*x1 + rnormal()
gen double Y = sin(x1)+.4*x2+rnormal()/5
replace Y = . in 3

// Reconstruct the fitted NP mean using precisely its complete-case sample.
np_dnn Y, xcols(x1 x2) x0(0 0) s(8)
scalar np_mean = el(e(estimates),1,1)
scalar np_ess = e(ess)
gen byte used = e(sample)
weights_dnn x1 x2 if used, x0(0 0) s(8) generate(w_exact)
scalar weight_n = r(n)
scalar weight_ess = r(neff)
assert abs(weight_ess-np_ess)<1e-10
quietly summarize w_exact, meanonly
assert abs(r(sum)-1)<1e-12
assert missing(w_exact) if !used
gen double weighted_y = w_exact*Y
quietly summarize weighted_y, meanonly
assert abs(r(sum)-np_mean)<1e-12

weights_dnn x1 x2 if used, x0(0 0) s(8) method(asymp) generate(w_asymp)
quietly summarize w_asymp, meanonly
assert abs(r(sum)-1)<1e-12
weights_dnn x1 x2 if used, x0(0 0) s(8) distance(mahalanobis) generate(w_mahal)
quietly summarize w_mahal, meanonly
assert abs(r(sum)-1)<1e-12

// Full order is a valid weight construction: exact k=1 selects one row.
weights_dnn x1 x2 if used, x0(0 0) s(`=scalar(weight_n)') generate(w_full)
quietly count if w_full>0 & !missing(w_full)
assert r(N)==1

matrix X0 = (-.5,0 \ 0,0 \ .5,0)
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix list e(tuning_candidates)
np_dnn_diagnostics, allcandidates
np_dnn_plot, type(tuning) allcandidates name(weights_paper_grid) replace nodraw

// Keep an infeasible nondefault row visible without selecting it for a plot.
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid) ///
    sgrid(4 8 59) sdefault(8)
matrix candidate_meta = e(tuning_candidates)
assert candidate_meta[2,3]==1
assert candidate_meta[3,4]==0 & candidate_meta[3,5]==2
np_dnn_diagnostics, allcandidates
np_dnn_plot, type(tuning) candidates(1 2) name(weights_direct_grid) replace nodraw

// Explicit estimate-only order n is valid; its SE must remain unavailable.
np_dnn Y, xcols(x1 x2) x0(0 0) s(59) nose
matrix no_se = e(se_targets)
assert missing(no_se[1,1])
di as result "WEIGHTS_AND_TUNING_EXAMPLE_PASSED"
