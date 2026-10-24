*! version 0.1.0  2026-10-18
// Isolate all derived-fit changes, including failures, from the caller's e().
program define _rdistnn_reuse, rclass
    version 17
    syntax, TARGETS(name) S(integer) COVARIANCE(string) [POST STORE(name) Level(cilevel) CANDIDATE(integer 0) VALIDATED]
    local notice_noisily ""
    if c(noisily) local notice_noisily "noisily"
    tempname original derived
    _estimates hold `original', copy restore
    capture noisily {
        quietly _rdistnn_model_fit, targets(`targets') s(`s') covariance(`covariance') candidate(`candidate') `validated'
        estimates store `derived'
        quietly rdistnn_results, level(`level')
        local quantity "`r(quantity)'"
        foreach field in fit se_fit vcov table summary index data {
            tempname out_`field'
            matrix `out_`field''=r(`field')
        }
        if "`quantity'"=="moment" {
            tempname theta
            matrix `theta'=r(theta_eval)
        }
        if "`store'"!="" estimates store `store'
    }
    local rc=_rc
    if !`rc' {
        capture `notice_noisily' _rdistnn_report_notices
        local rc=_rc
        if !`rc' {
            local n_tied_targets=r(n_tied_targets)
            local n_weak_targets=r(n_weak_targets)
        }
    }
    _estimates unhold `original'
    if `rc' {
        capture estimates drop `derived'
        exit `rc'
    }
    if "`post'"!="" estimates restore `derived'
    estimates drop `derived'
    foreach field in fit se_fit vcov table summary index data {
        return matrix `field'=`out_`field''
    }
    if "`quantity'"=="moment" return matrix theta_eval=`theta'
    return local quantity "`quantity'"
    return scalar n_tied_targets=`n_tied_targets'
    return scalar n_weak_targets=`n_weak_targets'
    return scalar level=`level'
    return scalar s=`s'
    if `candidate'>0 return scalar candidate=`candidate'
    if "`store'"!="" return local stored "`store'"
end
