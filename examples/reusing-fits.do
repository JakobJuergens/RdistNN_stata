// Extract a candidate, save its numeric state and predict after clearing data.
// Install rdistnn first, or run dev-setup.do from the repository root.
// This self-contained example clears the active dataset.
version 17
set more off
clear
set seed 42018
set obs 40
gen double x = rnormal()
gen double y = sin(x)+rnormal()/5
np_dnn y, xcols(x) x0(0) tuningmode(grid) sgrid(4 6 8) sdefault(6)
estimates store parent_fit

// Extraction preserves the active parent. Activate order 8 explicitly.
rdistnn_candidate, candidate(3) store(order8)
estimates restore order8
np_dnn
rdistnn_results
matrix report = r(data)
matrix list report

// Saved fits contain observation-level numeric inputs and signals.
// Apply the underlying data's sharing restrictions to the .ster file.
tempfile fitted
estimates save `fitted', replace
clear
estimates use `fitted'
matrix NEW = (-1\0\1)
rdistnn_predict, x0matrix(NEW) sefit store(predicted)
matrix prediction = r(fit)
matrix covariance = r(vcov)
matrix list prediction
matrix list covariance

// The derived fit can itself be reported or plotted, without active data.
estimates restore predicted
rdistnn_results, type(weights) targets(1)
np_dnn_plot, type(weights) view(kernel) name(saved_kernel) replace nodraw
