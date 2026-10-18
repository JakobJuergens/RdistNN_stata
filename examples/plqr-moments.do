// Experimental PLQR: supplied coefficients, localized moments and moment
// covariance only. Install rdistnn first, or place src/ on the adopath.
version 17
clear
set more off
set obs 160
set seed 20261017

// The median shape and regressor means are oracle nuisances in this DGP.
// In an application, merge aligned external predictions before this call.
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double md1 = .3*x1 - .2*x2
gen double md2 = -.1*x1 + .4*x2
gen double d1 = md1 + .4*rnormal()
gen double d2 = md2 + .4*rnormal()
gen double g_tau = .25*x1 - .15*x2
gen double Y = d1 - .5*d2 + g_tau + .5*rnormal()

// Independently expose the score and reporting weights for this one target.
// Equality belongs to the indicator's lower side: Y <= threshold.
gen double threshold = d1 - .5*d2 + g_tau
gen double quantile_residual = .5 - (Y <= threshold)
gen double psi1 = quantile_residual*(d1-md1)
gen double psi2 = quantile_residual*(d2-md2)
weights_dnn x1 x2, x0(0 0) k(1) s(12) generate(w)
gen double weighted_psi1 = w*psi1
gen double weighted_psi2 = w*psi2
quietly summarize weighted_psi1, meanonly
scalar manual_moment1 = r(sum)
quietly summarize weighted_psi2, meanonly
scalar manual_moment2 = r(sum)

plqr_dnn Y, x(x1 x2) x0(0 0) d(d1 d2) g(g_tau) md(md1 md2) ///
    thetaeval(1 -.5) tau(.5) k(1) s(12)
matrix moments = e(moment)
matrix moment_covariance = e(moment_vcov)
matrix supplied_theta = e(theta_eval)
assert e(N) == 160 & e(kd) == 2 & e(kx) == 2
assert rowsof(moments) == 1 & colsof(moments) == 2
assert rowsof(moment_covariance) == 2 & colsof(moment_covariance) == 2
assert supplied_theta[1,1] == 1 & supplied_theta[1,2] == -.5
assert abs(moments[1,1]-manual_moment1) < 1e-12
assert abs(moments[1,2]-manual_moment2) < 1e-12
assert abs(moment_covariance[1,2]-moment_covariance[2,1]) < 1e-12
assert moment_covariance[1,1] >= 0 & moment_covariance[1,1] < .
assert moment_covariance[2,2] >= 0 & moment_covariance[2,2] < .
assert "`e(theta_source)'" == "supplied"
assert "`e(inference_target)'" == "moment" & "`e(root_solver)'" == "none"
matrix list moments
matrix list moment_covariance
matrix list supplied_theta

// These SDs describe the localized score components, in regressor units.
// They are not coefficient SEs; no coefficient tests or intervals follow.
scalar moment1_sd = sqrt(moment_covariance[1,1])
scalar moment2_sd = sqrt(moment_covariance[2,2])
display "Moment-component jackknife SDs: " moment1_sd "  " moment2_sd
plqr_dnn
