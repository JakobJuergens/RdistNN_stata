// Reconstruct a local mean and compare localization orders.
// Install rdistnn first, or run dev-setup.do from the repository root.
// This self-contained example clears the active dataset.
version 17
set more off
clear
set obs 60
set seed 8731
gen double x1 = rnormal()
gen double x2 = .5*x1 + rnormal()
gen double Y = sin(x1)+.4*x2+rnormal()/5
replace Y = . in 3

// Match the estimator's sample and order when reusing its weights.
np_dnn Y, xcols(x1 x2) x0(0 0) s(8)
scalar fitted_mean = el(e(estimates),1,1)
gen byte used = e(sample)
weights_dnn x1 x2 if used, x0(0 0) s(8) generate(w_exact)
display "Effective sample size: " r(neff)
gen double weighted_y = w_exact*Y
quietly summarize weighted_y, meanonly
display "Fitted mean: " fitted_mean "  Reconstructed mean: " r(sum)
assert abs(r(sum)-fitted_mean)<1e-12

// The asymptotic rule changes weights; it is not an engine switch.
weights_dnn x1 x2 if used, x0(0 0) s(8) method(asymp) generate(w_asymp)
summarize w_exact w_asymp if used

// A grid compares estimates and retains a specified default, without selection.
matrix X0 = (-.5,0 \ 0,0 \ .5,0)
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid)
matrix list e(tuning_candidates)
np_dnn_plot, type(tuning) allcandidates name(weights_paper_grid) replace nodraw

// Use explicit candidate orders when they match the intended sensitivity check.
np_dnn Y, xcols(x1 x2) x0matrix(X0) tuningmode(grid) ///
    sgrid(4 8 16) sdefault(8)
np_dnn_diagnostics, allcandidates
np_dnn_plot, type(tuning) allcandidates name(weights_direct_grid) replace nodraw
