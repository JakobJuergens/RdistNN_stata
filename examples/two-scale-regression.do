// Two-scale regression with signed weights; replaces the active dataset.
// Install rdistnn first, or run dev-setup.do from the repository root.
version 17
clear
set more off
set seed 20261018
set obs 160
gen double x1=2*runiform()-1
gen double x2=rnormal()
gen double Y=1+x1+x1^2-.4*x2+rnormal()*.3
matrix X0=(-.75,0 \ 0,0 \ .75,0)
np_tdnn Y, xcols(x1 x2) x0matrix(X0) s1(4) s2(12) ///
    distance(l2_scaled) scalerule(sd)
// The combined weights can be negative; they are not probability mass.
matrix summary=r(summary)
matrix list summary
np_tdnn_diagnostics
np_tdnn_plot, type(estimates) name(tdnn_estimates) replace nodraw
np_tdnn_plot, type(weights) name(tdnn_weights) replace nodraw
rdistnn_results, type(vcov)
matrix covariance=r(vcov)
rdistnn_results, type(weights)
matrix observation_weights=r(data)
// New profiles reuse both orders and the fitted distance geometry.
matrix NEW=(-.5,0 \ .5,0)
rdistnn_predict, x0matrix(NEW) sefit store(tdnn_new_profiles)
matrix prediction=r(fit)
matrix prediction_covariance=r(vcov)
weights_tdnn x1 x2, x0matrix(X0) s1(4) s2(12) ///
    distance(l2_scaled) scalerule(sd)
matrix weights=r(weights)
matrix ranks=r(ranks)
matrix distances=r(distances)
