*! version 0.1.0  2026-10-18
// Post a derived fit from frozen numeric state; callers control e() restoration.
program define _rdistnn_model_fit, eclass
    version 17
    syntax, TARGETS(name) S(integer) [COVariance(string) CANDIDATE(integer 0)]
    if !inlist("`e(cmd)'","np_dnn","np_tdnn","cate_dnn","clate_dnn","plqr_dnn") error 301
    if "`e(cmd)'"=="cate_dnn" & !inlist("`e(score_method)'","aipw","pseudo_outcome") {
        di as error "Current CATE saved results require AIPW or supplied pseudo-outcomes; refit cate_dnn."
        exit 498
    }
    if e(model_version)!=1 {
        di as error "Active estimates have no supported saved model; refit with the current package version."
        exit 498
    }
    local matrices : e(matrices)
    foreach field in model_X model_signals model_inputs model_ranks model_distances model_targets model_row_ids {
        if !`: list field in matrices' {
            di as error "Saved model is missing e(`field'); refit with the current package version."
            exit 498
        }
    }
    if "`covariance'"=="" local covariance "compute"
    if !inlist("`covariance'","compute","none") error 198
    local active=("`covariance'"=="compute")
    local cmd "`e(cmd)'"
    local lambda=e(lambda)
    if `candidate'>0 local lambda=el(e(tuning_candidates),`candidate',7)
    local p=colsof(`targets')
    local m=rowsof(`targets')
    if `p'!=colsof(e(model_X)) {
        di as error "Targets must have one column per fitted localization covariate."
        exit 198
    }
    foreach file in _rdistnn_mata.ado _rdistnn_clate_mata.ado _rdistnn_model_mata.ado {
        findfile "`file'"
        quietly run "`r(fn)'"
    }
    tempname b V diagnostics components componentV status neff ranks DX W fit ses
    mata: rdistnn_model_evaluate("`targets'",`s',`active',tokens("`b' `V' `diagnostics' `components' `componentV' `status' `neff' `ranks' `DX' `W' `fit'"))
    mata: st_matrix("`ses'",sqrt(diagonal(st_matrix("`V'")))')
    _rdistnn_model_post, targets(`targets') s(`s') active(`active') candidate(`candidate') ///
        b(`b') v(`V') diagnostics(`diagnostics') components(`components') componentv(`componentV') ///
        status(`status') neff(`neff') ranks(`ranks') dx(`DX') w(`W') fit(`fit') ses(`ses') ///
        lambda(`lambda') p(`p') m(`m')
end
