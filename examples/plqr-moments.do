// Experimental PLQR moments at supplied coefficients; no coefficient fitting.
// Install rdistnn first, or run dev-setup.do from the repository root.
// This self-contained example clears the active dataset.
version 17
set more off
clear
set obs 160
set seed 20261017
// g_tau and md1/md2 are oracle shape and regressor-mean nuisances.
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double md1 = .3*x1 - .2*x2
gen double md2 = -.1*x1 + .4*x2
gen double d1 = md1 + .4*rnormal()
gen double d2 = md2 + .4*rnormal()
gen double g_tau = .25*x1 - .15*x2
gen double Y = d1 - .5*d2 + g_tau + .5*rnormal()

// In an application, merge aligned external predictions before this call.
plqr_dnn Y, x(x1 x2) x0(0 0) d(d1 d2) g(g_tau) md(md1 md2) ///
    thetaeval(1 -.5) tau(.5) k(1) s(12)
matrix moments = e(moment)
matrix moment_covariance = e(moment_vcov)
matrix supplied_theta = e(theta_eval)
matrix list moments
matrix list moment_covariance
matrix list supplied_theta

// These SDs describe score components, not the supplied coefficients.
display "Moment-component SDs: " sqrt(moment_covariance[1,1]) ///
    "  " sqrt(moment_covariance[2,2])

// Saved-fit prediction can evaluate several profiles with joint moment covariance.
matrix NEW = (-.5,0 \ 0,0 \ .5,0)
rdistnn_predict, x0matrix(NEW) sefit store(plqr_profiles)
matrix predicted_moments = r(fit)
matrix joint_moment_covariance = r(vcov)
matrix list predicted_moments
// Supplied theta remains separate; neither command solves for coefficients.
