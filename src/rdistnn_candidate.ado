*! version 0.1.5  2026-10-25
program define rdistnn_candidate, rclass
    version 17
    syntax, CANDIDATE(integer) [POST STORE(name) Level(cilevel)]
    if !inlist("`e(cmd)'","np_dnn","np_tdnn","cate_dnn","clate_dnn","plqr_dnn") error 301
    _rdistnn_model_check
    local matrices : e(matrices)
    local field "tuning_candidates"
    if !`:list field in matrices' {
        _rdistnn_messages error, rc(498) message(`"This fit has no retained candidate grid."')
    }
    if `candidate'<1 | `candidate'>rowsof(e(tuning_candidates)) error 198
    if el(e(tuning_candidates),`candidate',4)!=1 {
        _rdistnn_messages error, rc(198) message(`"Requested candidate is infeasible."')
    }
    local s=el(e(tuning_candidates),`candidate',8)
    local covariance "compute"
    if "`e(se_method)'"=="none" local covariance "none"
    tempname targets
    matrix `targets'=e(model_targets)
    if "`store'"!="" local storeopt "store(`store')"
    _rdistnn_reuse, targets(`targets') s(`s') covariance(`covariance') candidate(`candidate') level(`level') `post' `storeopt' validated
    return add
end
