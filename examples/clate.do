// Experimental CLATE with oracle instrument-arm predictions.
// Run after installing the feature/clate development package, or with src/
// on the adopath. This example illustrates mechanics, not fitted learners.
version 17
clear
set more off
set obs 160
set seed 4711
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double e = invlogit(.3*x1-.2*x2)
gen byte A = runiform()<e
gen double my0 = .5*x1-.2*x2
gen double md1 = invlogit(.4+.2*x1)
gen double md0 = 0
gen double tau = 1+.3*sin(x1)
gen byte D = A*(runiform()<md1)
gen double Y = my0 + tau*D + rnormal()*.3
gen double my1 = my0 + md1*tau
matrix X0 = (-.5,0 \ 0,0 \ .5,0)
clate_dnn Y, xcols(x1 x2) x0matrix(X0) dcol(D) acol(A) ecol(e) ///
    my0col(my0) my1col(my1) md1col(md1) onesided ///
    training(oracle) tuningmode(grid)
matrix list e(components)
matrix list e(status)
clate_dnn_diagnostics, detail
clate_dnn_plot, type(estimates) name(clate_example_estimates) replace nodraw
clate_dnn_plot, type(overlap) view(mass) name(clate_example_overlap) replace nodraw
clate_dnn_plot, type(influence) name(clate_example_influence) replace nodraw
clate_dnn_plot, type(tuning) allcandidates name(clate_example_tuning) replace nodraw
