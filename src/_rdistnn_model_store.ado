*! version 0.1.0  2026-10-18
// Numeric snapshots survive estimates save/use and do not depend on Mata globals.
program define _rdistnn_model_store, eclass
    version 17
    syntax varlist(numeric), SIGNALS(name) RANKS(name) DISTANCES(name) GROUPS(name) [GEOMETRY(name)]
    local inputs : list uniq varlist
    local xvars "`e(xvars)'"
    tempname data X row_ids targets
    tempvar sample
    quietly gen byte `sample'=e(sample)
    mata: st_matrix("`data'",st_data(.,tokens("`inputs'"),"`sample'"))
    mata: st_matrix("`X'",st_data(.,tokens("`xvars'"),"`sample'"))
    mata: st_matrix("`row_ids'",select((1::st_nobs()),st_data(.,"`sample'")))
    matrix colnames `data'=`inputs'
    matrix colnames `X'=`xvars'
    matrix colnames `row_ids'=observation_id
    matrix `targets'=e(targets)
    local n=e(N)
    foreach value in `data' `X' `row_ids' `signals' `ranks' `distances' {
        if rowsof(`value')!=`n' {
            di as error "Cannot save fitted model: snapshot rows do not match e(N)."
            exit 498
        }
    }
    ereturn matrix model_inputs=`data', copy
    ereturn matrix model_X=`X', copy
    ereturn matrix model_signals=`signals', copy
    ereturn matrix model_ranks=`ranks', copy
    ereturn matrix model_distances=`distances', copy
    ereturn matrix model_row_ids=`row_ids', copy
    ereturn matrix model_targets=`targets', copy
    if "`geometry'"!="" {
        if rowsof(`geometry')>0 & colsof(`geometry')>0 {
            ereturn matrix model_geometry=`geometry', copy
        }
    }
    else if inlist("`e(distance)'","mahalanobis","l1_scaled","l2_scaled") {
        di as error "Fitted scaled/Mahalanobis models require explicit frozen geometry()."
        exit 498
    }
    ereturn matrix model_groups=`groups', copy
    ereturn scalar model_version=3
    ereturn local tie_rule "exact_distance_block_averaging_v1"
    ereturn local model_vars "`inputs'"
end
