// CATE input routes, pointwise reporting and a conditional-mean fit.
// Install rdistnn first, or run dev-setup.do from the repository root.
// This self-contained example clears the active dataset.
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
gen double Y = m0_oracle + D*tau_oracle + .3*rnormal()
matrix X0 = (0,0 \ .5,0 \ -.5,0)

// AIPW from supplied oracle predictions. Copy r() tables before another command.
cate_dnn Y, xcols(x1 x2) x0matrix(X0) tcol(D) ///
    ecol(e_oracle) m0col(m0_oracle) m1col(m1_oracle) tuningmode(grid) level(95)
matrix AIPW_summary = r(summary)
matrix AIPW_b = e(b)
matrix AIPW_V = e(V)
estimates store cate_aipw
matrix list AIPW_summary
cate_dnn_diagnostics, targets(1 3)

// A supplied AIPW signal gives the same estimates on the same aligned sample.
gen double eclip = min(max(e_oracle,1e-8),1-1e-8)
gen double rho_aipw = m1_oracle-m0_oracle + D*(Y-m1_oracle)/eclip ///
    - (1-D)*(Y-m0_oracle)/(1-eclip)
cate_dnn rho_aipw, xcols(x1 x2) x0matrix(X0) pseudooutcome
matrix supplied_b = e(b)
matrix route_difference = supplied_b-AIPW_b
matrix list route_difference

// Replay changes the pointwise interval level, without refitting.
estimates restore cate_aipw
cate_dnn, level(90)

// Covariance supplies the arithmetic for a prespecified contrast.
// Valid contrast inference requires the corresponding joint assumptions.
matrix L = (1,-1,0)
matrix difference = AIPW_b*L'
matrix difference_variance = L*AIPW_V*L'
matrix list difference
matrix list difference_variance

// Complete-case NP regression after an explicit sample restriction.
gen double Y_np = Y
replace Y_np = . in 180
np_dnn Y_np if _n > 5, xcols(x1 x2) x0matrix(X0)
matrix mean_report = r(summary)
matrix list mean_report
count if e(sample)
estimates store np_fit

// A supplied-score fit can omit uncertainty; missing SEs are not zero SEs.
cate_dnn rho_aipw, xcols(x1 x2) x0(0 0) pseudooutcome nose
matrix estimate_only = r(summary)
matrix list estimate_only
estimates restore np_fit
