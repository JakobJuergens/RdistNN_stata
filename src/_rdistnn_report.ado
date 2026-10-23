*! version 0.1.0  2026-10-17
// Full-target reporting, independent of the estimable subset posted in e(b).
program define _rdistnn_report, rclass
    version 17
    syntax [, Level(cilevel)]
    _rdistnn_model_check
    tempname estimates ses covariance table summary status b se crit
    matrix `estimates' = e(estimates)
    matrix `ses' = e(se_targets)
    local m = colsof(`estimates')
    scalar `crit' = invnormal(.5+`level'/200)
    matrix `table' = J(9,`m',.)
    if "`e(cmd)'" == "clate_dnn" matrix `status' = e(status)
    forvalues j=1/`m' {
        // Numeric scalars preserve double precision. Converting large or
        // tiny values through local text macros rounds the returned tables.
        scalar `b' = el(`estimates',1,`j')
        scalar `se' = el(`ses',1,`j')
        matrix `table'[1,`j'] = scalar(`b')
        matrix `table'[2,`j'] = scalar(`se')
        matrix `table'[8,`j'] = scalar(`crit')
        matrix `table'[9,`j'] = 0
        local wald = !missing(scalar(`b'),scalar(`se'))
        if "`e(cmd)'" == "clate_dnn" local wald = (el(`status',`j',5)==1)
        if `wald' {
            matrix `table'[3,`j'] = scalar(`b')/scalar(`se')
            if scalar(`se') > 0 matrix `table'[4,`j'] = 2*normal(-abs(scalar(`b')/scalar(`se')))
            // R represents a nonzero/zero z statistic as infinity; Stata uses missing.
            else if scalar(`b') != 0 matrix `table'[4,`j'] = 0
            matrix `table'[5,`j'] = scalar(`b')-scalar(`crit')*scalar(`se')
            matrix `table'[6,`j'] = scalar(`b')+scalar(`crit')*scalar(`se')
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
    _rdistnn_report_notices
    return matrix table = `table'
    return matrix summary = `summary'
    return scalar level = `level'
end
