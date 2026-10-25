*! version 0.1.5  2026-10-25
program define _rdistnn_estimation_state, eclass
    version 17
    syntax varlist(numeric)
    tempname estimates ses covariance
    local m = e(m)
    if "`e(cmd)'" == "clate_dnn" {
        matrix `estimates' = e(estimates)
        matrix `ses' = e(se_targets)
        // Retain the paired estimator covariance in the canonical result.
        matrix `covariance' = e(theta_vcov)
    }
    else {
        matrix `estimates' = e(b)
        matrix `ses' = J(1,`m',.)
        matrix `covariance' = J(`m',`m',.)
        if "`e(se_method)'" != "none" {
            matrix `covariance' = e(V)
            forvalues j=1/`m' {
                matrix `ses'[1,`j'] = sqrt(el(`covariance',`j',`j'))
            }
        }
        matrix colnames `ses' = `: colnames `estimates''
        matrix colnames `covariance' = `: colnames `estimates''
        matrix rownames `covariance' = `: colnames `estimates''
    }
    ereturn matrix estimates = `estimates', copy
    ereturn matrix se_targets = `ses', copy
    ereturn matrix estimate_vcov = `covariance', copy
    ereturn local tie_rule "exact_distance_block_averaging_v1"
    quietly signestimationsample `varlist'
    quietly _rdistnn_sample_signature `varlist' if e(sample)
    ereturn local ordered_signature "`r(signature)'"
    ereturn local signature_vars "`varlist'"
end
