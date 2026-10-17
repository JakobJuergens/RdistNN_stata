*! version 0.1.0  2026-10-17
program define _rdistnn_report_header
    version 17
    local title "DNN Nonparametric Regression"
    if "`e(cmd)'" == "cate_dnn" local title "DNN Conditional Average Treatment Effect Estimator"
    if "`e(cmd)'" == "clate_dnn" local title "DNN Conditional Local Average Treatment Effect (experimental)"
    if "`e(cmd)'" == "plqr_dnn" local title "Experimental DNN PLQR Localized-Moment Evaluator"
    local m = e(m)
    if missing(`m') local m = 1
    local p = e(p)
    if missing(`p') local p = e(kx)
    di _n as text "`title'"
    di as text "Observations: " as result e(N) as text "  Targets: " as result `m' ///
        as text "  Covariates: " as result `p'
    if `m' == 1 & "`e(x0)'" != "" di as text "Evaluation point x0: " as result "`e(x0)'"
    else di as text "Evaluation targets: " as result `m' as text "  (see e(targets))"
    di as text "Covariates: " as result "`e(xvars)'"
    di as text "Method: " as result "`e(method)'" as text "  Distance: " as result "`e(distance)'" ///
        as text "  k = " as result e(k) as text "  s = " as result e(s)
    if `m' == 1 di as text "Effective sample size: " as result %9.3g e(neff)
    else di as text "Effective sample size range: " as result %9.3g e(neff_min) ///
        as text " to " as result %9.3g e(neff_max)
    if "`e(tuning_mode)'" != "" di as text "Tuning: " as result "`e(tuning_mode)'" ///
        as text "  Covariance: " as result "`e(se_method)'"
    if "`e(se_method)'" == "none" di as text "Inference was not requested (nose); uncertainty is unavailable."
end
