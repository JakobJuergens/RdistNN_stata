*! version 0.1.5  2026-10-25
program define weights_tdnn, rclass
    version 17
    syntax varlist(numeric) [if] [in], [X0(numlist) X0MATRIX(name) ///
        S1(integer 1) S2(integer 2) DIM(string) METHOD(string) ///
        noNORMALIZE DISTance(string) SCALERule(string) SCALEVec(numlist) ///
        SCALEEPS(real 1e-12) GENerate(name)]
    marksample touse
    markout `touse' `varlist'
    quietly count if `touse'
    local n=r(N)
    if `n'==0 error 2000
    local p : word count `varlist'
    if "`dim'"=="" local dim=`p'
    confirm integer number `dim'
    if `s1'<1 | `s1'>=`s2' | `s2'>`n' | `dim'<1 {
        _rdistnn_messages error, rc(198) message(`"Require 1 <= s1 < s2 <= n and positive dim()."')
    }
    if "`method'"=="" local method "exact"
    if !inlist("`method'","exact","asymp") error 198
    if "`distance'"=="" local distance "l2"
    if "`scalerule'"=="" local scalerule "sd"
    if ("`x0'"=="")==("`x0matrix'"=="") {
        _rdistnn_messages error, rc(198) message(`"Specify exactly one of x0() or x0matrix()."')
    }
    tempname targets W R DX groups table diag geom row_ids target_notices notice_counts
    if "`x0matrix'"!="" {
        confirm matrix `x0matrix'
        matrix `targets'=`x0matrix'
    }
    else {
        local row : subinstr local x0 " " ",", all
        matrix `targets'=(`row')
    }
    local m=rowsof(`targets')
    if `m'<1 | colsof(`targets')!=`p' {
        local target_rc = cond("`x0matrix'"=="",198,503)
        _rdistnn_messages error, rc(`target_rc') message(`"Targets must have one column per localization covariate."')
    }
    if "`generate'"!="" {
        confirm new variable `generate'
        if `m'!=1 {
            _rdistnn_messages error, rc(198) message(`"generate() requires a single target; use r(weights) for batches."')
        }
    }
    local normalize=("`normalize'"=="")
    foreach file in _rdistnn_mata.ado _rdistnn_tdnn_mata.ado {
        findfile "`file'"
        quietly run "`r(fn)'"
    }
    mata: rdistnn_tdnn_weights_ado()
    capture confirm matrix `geom'
    local geometry_rc=_rc
    if !inlist(`geometry_rc',0,111) exit `geometry_rc'
    mata: st_matrix("`target_notices'",rdistnn_target_notices(st_matrix("`groups'"),st_matrix("`diag'")[,5]))
    matrix colnames `target_notices'=distance_ties weak_localization
    mata: st_matrix("`notice_counts'",colsum(st_matrix("`target_notices'")))
    if "`generate'"!="" {
        quietly gen double `generate'=.
        mata: st_store(.,"`generate'","`touse'",st_matrix("`W'"))
    }
    matrix colnames `diag'=signed_sum positive_mass negative_mass abs_mass ///
        inverse_l2_concentration max_abs_weight n_nonzero_weight abs_mean_distance
    matrix colnames `row_ids'=observation_id
    return matrix weights=`W'
    return matrix target_notices=`target_notices'
    return scalar n_tied_targets=el(`notice_counts',1,1)
    return scalar n_weak_targets=el(`notice_counts',1,2)
    return matrix ranks=`R'
    return matrix distances=`DX'
    return matrix groups=`groups'
    return matrix rank_weights=`table'
    return matrix targets=`targets'
    return matrix diagnostics=`diag'
    return matrix row_ids=`row_ids'
    if `geometry_rc'==0 return matrix geometry=`geom'
    return scalar n=`n'
    return scalar m=`m'
    return scalar p=`p'
    return scalar s1=`s1'
    return scalar s2=`s2'
    return scalar dim=`dim'
    return scalar normalize=`normalize'
    return local method "`method'"
    return local distance "`distance'"
    return local tie_rule "exact_distance_block_averaging_v1"
end
