*! version 0.1.5  2026-10-25
program define _rdistnn_report_header
    version 17
    local title "DNN Nonparametric Regression"
    if "`e(cmd)'" == "np_tdnn" local title "Two-scale DNN Nonparametric Regression"
    if "`e(cmd)'" == "cate_dnn" local title "DNN Conditional Average Treatment Effect Estimator"
    if "`e(cmd)'" == "clate_dnn" local title "DNN Conditional Local Average Treatment Effect"
    if "`e(cmd)'" == "plqr_dnn" local title "Experimental DNN PLQR Localized-Moment Evaluator"
    local m = e(m)
    local p = e(p)
    di _n as text "`title'"
    di as text "Observations: " as result e(N) as text "  Targets: " as result `m' ///
        as text "  Covariates: " as result `p'
    if `m' == 1 & "`e(x0)'" != "" di as text "Evaluation point x0: " as result "`e(x0)'"
    else di as text "Evaluation targets: " as result `m' as text "  (see e(targets))"
    di as text "Covariates: " as result "`e(xvars)'"
    di as text "Method: " as result "`e(method)'" as text "  Distance: " as result "`e(distance)'"
    local concentration "Effective sample size"
    if "`e(cmd)'" == "np_tdnn" {
        di as text "Orders s1 = " as result e(s1) as text "  s2 = " as result e(s2) ///
            as text "  Localization dimension = " as result e(dim)
        di as text "s1/n = " as result %6.3g (e(s1)/e(N)) ///
            as text "  s2/n = " as result %6.3g (e(s2)/e(N))
        local concentration "Inverse squared-weight concentration"
    }
    else di as text "k = " as result e(k) as text "  s = " as result e(s) ///
        as text "  s/n = " as result %6.3g (e(s)/e(N))
    if `m' == 1 di as text "`concentration': " as result %9.3g e(neff)
    else di as text "`concentration' range: " as result %9.3g e(neff_min) ///
        as text " to " as result %9.3g e(neff_max)
    if "`e(tuning_mode)'" != "" di as text "Tuning: " as result "`e(tuning_mode)'" ///
        as text "  Covariance: " as result "`e(se_method)'"
    if "`e(se_method)'" == "none" di as text "Inference was not requested (nose); uncertainty is unavailable."
end
