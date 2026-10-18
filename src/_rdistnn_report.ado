*! version 0.1.0  2026-10-17
// Full-target reporting, independent of the estimable subset posted in e(b).
program define _rdistnn_report, rclass
    version 17
    syntax [, Level(cilevel)]
    tempname estimates ses covariance table summary status
    local matrices : e(matrices)
    if strpos(" `matrices' "," estimates ") {
        matrix `estimates' = e(estimates)
        matrix `ses' = e(se_targets)
    }
    else {
        // Older CATE/NP estimates remain replayable without altering stored e().
        matrix `estimates' = e(b)
        matrix `ses' = J(1,colsof(`estimates'),.)
        if "`e(se_method)'" != "none" & strpos(" `matrices' "," V ") {
            matrix `covariance' = e(V)
            forvalues j=1/`=colsof(`estimates')' {
                matrix `ses'[1,`j'] = sqrt(el(`covariance',`j',`j'))
            }
        }
    }
    local m = colsof(`estimates')
    local crit = invnormal(.5+`level'/200)
    matrix `table' = J(9,`m',.)
    if "`e(cmd)'" == "clate_dnn" matrix `status' = e(status)
    forvalues j=1/`m' {
        local b = el(`estimates',1,`j')
        local se = el(`ses',1,`j')
        matrix `table'[1,`j'] = `b'
        matrix `table'[2,`j'] = `se'
        matrix `table'[8,`j'] = `crit'
        matrix `table'[9,`j'] = 0
        local wald = !missing(`b',`se')
        if "`e(cmd)'" == "clate_dnn" local wald = (el(`status',`j',5)==1)
        if `wald' {
            matrix `table'[3,`j'] = `b'/`se'
            if `se' > 0 matrix `table'[4,`j'] = 2*normal(-abs(`b'/`se'))
            // R represents a nonzero/zero z statistic as infinity; Stata uses missing.
            else if `b' != 0 matrix `table'[4,`j'] = 0
            matrix `table'[5,`j'] = `b'-`crit'*`se'
            matrix `table'[6,`j'] = `b'+`crit'*`se'
        }
    }
    matrix rownames `table' = b se z pvalue ll ul df crit eform
    matrix colnames `table' = `: colnames `estimates''
    matrix `summary' = `table'[1..6,1..`m']'
    matrix colnames `summary' = estimate std_error z pvalue ci_lower ci_upper
    _rdistnn_report_header
    if "`e(cmd)'"=="plqr_dnn" {
        di as text "Localized score moments at supplied theta_eval; no coefficient inference or root solving."
        di as text "Quantile index tau: " as result %6.3f e(tau)
        matrix list e(theta_eval), noheader format(%12.6g)
    }
    if "`e(cmd)'" == "clate_dnn" {
        di as text "Paired reduced-form / first-stage ratio; pointwise Wald inference."
        matrix list e(components), noheader format(%12.6g)
    }
    di as text "Pointwise confidence level: `level'%"
    matrix list `summary', noheader format(%12.6g)
    if "`e(cmd)'" == "clate_dnn" {
        forvalues j=1/`m' {
            if el(`status',`j',1)!=0 di as text "note: target `j' ratio is unavailable; see e(status) and e(components)."
            if el(`status',`j',2)==1 di as text "note: target `j' has a negative first stage."
            if el(`status',`j',3)==1 di as text "note: target `j' first stage is near zero; the ratio is not regularized."
            if el(`status',`j',5)==0 di as text "note: target `j' has no pointwise Wald interval."
        }
    }
    return matrix table = `table'
    return matrix summary = `summary'
    return scalar level = `level'
end
