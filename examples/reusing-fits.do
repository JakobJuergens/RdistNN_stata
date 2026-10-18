// Saved numeric fits; safe to run as a standalone example.
version 17
clear
set more off
set seed 42018
set obs 40
gen double x=rnormal()
gen double y=sin(x)+rnormal()/5
np_dnn y, xcols(x) x0(0) tuningmode(grid) sgrid(4 6 8) sdefault(6)
rdistnn_candidate, candidate(3) store(order8)
rdistnn_results
matrix report=r(data)
tempfile fitted
estimates save `fitted', replace
clear
estimates use `fitted'
matrix NEW=(-1\0\1)
rdistnn_predict, x0matrix(NEW) sefit store(predicted)
matrix prediction=r(fit)
matrix covariance=r(vcov)
estimates restore predicted
rdistnn_results, type(weights) targets(1)
np_dnn_plot, type(weights) view(kernel) name(saved_kernel) replace nodraw
