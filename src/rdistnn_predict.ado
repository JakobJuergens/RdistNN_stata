*! version 0.1.5  2026-10-25
program define rdistnn_predict, rclass
    version 17
    syntax [, X0(numlist) X0Matrix(name) SEFIT INTERVAL(string) Level(cilevel) POST STORE(name)]
    if !inlist("`e(cmd)'","np_dnn","np_tdnn","cate_dnn","clate_dnn","plqr_dnn") error 301
    if ("`x0'"=="")==("`x0matrix'"=="") {
        _rdistnn_messages error, rc(198) message(`"Specify exactly one of x0() and x0matrix()."')
    }
    if "`interval'"!="" & "`interval'"!="confidence" error 198
    _rdistnn_model_check
    tempname targets
    if "`x0matrix'"!="" {
        confirm matrix `x0matrix'
        matrix `targets'=`x0matrix'
    }
    else mata: st_matrix("`targets'",strtoreal(tokens("`x0'")))
    local covariance "none"
    if "`sefit'`interval'"!="" local covariance "compute"
    if "`store'"!="" local storeopt "store(`store')"
    _rdistnn_reuse, targets(`targets') s(`=e(s)') covariance(`covariance') level(`level') `post' `storeopt' validated
    return add
end
