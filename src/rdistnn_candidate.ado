*! version 0.1.0  2026-10-18
program define rdistnn_candidate, rclass
    version 17
    syntax, CANDIDATE(integer) [POST STORE(name) Level(cilevel)]
    if !inlist("`e(cmd)'","np_dnn","cate_dnn","clate_dnn","plqr_dnn") error 301
    local matrices : e(matrices)
    local field "tuning_candidates"
    if !`:list field in matrices' {
        di as error "This fit has no retained candidate grid."
        exit 498
    }
    if `candidate'<1 | `candidate'>rowsof(e(tuning_candidates)) error 198
    if el(e(tuning_candidates),`candidate',4)!=1 {
        di as error "Requested candidate is infeasible."
        exit 198
    }
    local s=el(e(tuning_candidates),`candidate',8)
    local covariance "compute"
    if "`e(se_method)'"=="none" local covariance "none"
    tempname targets
    matrix `targets'=e(model_targets)
    if "`store'"!="" local storeopt "store(`store')"
    _rdistnn_reuse, targets(`targets') s(`s') covariance(`covariance') candidate(`candidate') level(`level') `post' `storeopt'
    return add
end
