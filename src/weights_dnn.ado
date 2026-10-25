*! version 0.1.5  2026-10-25  Jakob Juergens
*! Exact/asymptotic DNN weights for one target or a target matrix
program define weights_dnn, rclass
    version 17
    syntax varlist(numeric min=1) [if] [in], [ ///
        X0(numlist) X0MATRIX(name) K(integer 1) S(integer 0) ///
        Generate(name) DISTance(string) METHod(string) ///
        SCALERule(string) SCALEVec(numlist) SCALEEPS(real 1e-12)]
    marksample touse
    markout `touse' `varlist'
    local kx : word count `varlist'
    quietly count if `touse'
    local n=r(N)
    if `n'==0 error 2000
    if `k'<1 {
        _rdistnn_messages error, rc(198) message("k() must be a positive integer; got `k'.") hint("Set k() to at least 1.")
    }
    if ("`x0'"=="")==("`x0matrix'"=="") {
        _rdistnn_messages error, rc(198) message(`"Specify exactly one of x0() and x0matrix()."')
    }
    tempname targets weights ranks distances groups row_ids neff geometry target_notices notice_counts
    if "`x0matrix'"!="" {
        confirm matrix `x0matrix'
        matrix `targets'=`x0matrix'
    }
    else {
        local nx0 : word count `x0'
        if `nx0'!=`kx' {
            _rdistnn_messages error, rc(198) message("x0() requires `kx' values; got `nx0'.") hint("Supply one value per localization covariate.")
        }
        mata: st_matrix("`targets'",strtoreal(tokens("`x0'")))
    }
    local m=rowsof(`targets')
    if colsof(`targets')!=`kx' | `m'<1 {
        _rdistnn_messages error, rc(503) message(`"x0matrix() must have one target per row and one column per covariate."')
    }
    if "`generate'"!="" {
        confirm new variable `generate'
        if `m'!=1 {
            _rdistnn_messages error, rc(198) message(`"generate() requires a single target; use r(weights) for a batch."')
        }
    }
    if "`distance'"=="" local distance "l2"
    if "`method'"=="" local method "exact"
    if "`scalerule'"=="" local scalerule "sd"
    if !inlist("`method'","exact","asymp") {
        _rdistnn_messages error, rc(198) message(`"method() must be exact or asymp."')
    }
    if "`distance'"=="mahalanobis" & `kx'<2 {
        _rdistnn_messages error, rc(198) message(`"distance(mahalanobis) requires at least 2 covariates."')
    }
    if `s'==0 local s=floor(.4*`n'^((`kx'+1)/(`kx'+2)))
    if `s'<`k' | `s'>`n' {
        _rdistnn_messages error, rc(198) message(`"s() must satisfy k <= s <= n."')
    }
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"
    mata: rdistnn_weights_fit_ado()
    mata: st_matrix("`target_notices'",rdistnn_target_notices(st_matrix("`groups'"),(1:/colsum(st_matrix("`weights'"):^2))'))
    matrix colnames `target_notices'=distance_ties weak_localization
    mata: st_matrix("`notice_counts'",colsum(st_matrix("`target_notices'")))
    if "`generate'"!="" {
        quietly gen double `generate'=.
        mata: st_store(.,"`generate'","`touse'",st_matrix("`weights'"))
    }
    if inlist("`distance'","mahalanobis","l1_scaled","l2_scaled") {
        return matrix geometry=`geometry'
    }
    return matrix weights=`weights'
    return matrix target_notices=`target_notices'
    return scalar n_tied_targets=el(`notice_counts',1,1)
    return scalar n_weak_targets=el(`notice_counts',1,2)
    return matrix ranks=`ranks'
    return matrix distances=`distances'
    return matrix groups=`groups'
    return matrix row_ids=`row_ids'
    return matrix targets=`targets'
    return scalar neff=`neff'[1,1]
    return matrix neff_targets=`neff'
    return scalar n=`n'
    return scalar m=`m'
    return scalar k=`k'
    return scalar s=`s'
    return scalar kx=`kx'
    return scalar scaleeps=`scaleeps'
    return local distance "`distance'"
    return local scalerule "`scalerule'"
    return local method "`method'"
    return local tie_rule "exact_distance_block_averaging_v1"
    di as text "DNN weights: " as result `n' as text " observations, " as result `m' as text " target(s)."
end
