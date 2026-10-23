*! version 0.1.0  2026-10-18  Jakob Juergens
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
        di as error "Specify exactly one of x0() and x0matrix()."
        exit 198
    }
    tempname targets weights ranks distances groups row_ids neff geometry
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
        di as error "x0matrix() must have one target per row and one column per covariate."
        exit 198
    }
    if "`generate'"!="" {
        confirm new variable `generate'
        if `m'!=1 {
            di as error "generate() requires a single target; use r(weights) for a batch."
            exit 198
        }
    }
    if "`distance'"=="" local distance "l2"
    if "`method'"=="" local method "exact"
    if "`scalerule'"=="" local scalerule "sd"
    if !inlist("`method'","exact","asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }
    if "`distance'"=="mahalanobis" & `kx'<2 {
        di as error "distance(mahalanobis) requires at least 2 covariates."
        exit 198
    }
    if `s'==0 local s=floor(.4*`n'^((`kx'+1)/(`kx'+2)))
    if `s'<`k' | `s'>`n' {
        di as error "s() must satisfy k <= s <= n."
        exit 198
    }
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"
    mata {
        _wd_X=st_data(.,tokens("`varlist'"),"`touse'")
        _wd_X0=st_matrix("`targets'")
        _wd_G=rdistnn_fit_geometry(_wd_X,"`distance'","`scalerule'",strtoreal(tokens("`scalevec'")),`scaleeps')
        _wd_D=rdistnn_distance_matrix(_wd_X,_wd_X0,"`distance'",_wd_G)
        _wd_R=rdistnn_ranks_from_dmat(_wd_D)
        _wd_groups=rdistnn_distance_groups(_wd_R,_wd_D,1)
        rdistnn_distance_warning(_wd_groups)
        _wd_W=rdistnn_weights_from_rmat(_wd_R,`k',`s',"`method'",_wd_D,_wd_groups)
    }
    mata {
        st_matrix("`weights'",_wd_W); st_matrix("`ranks'",_wd_R)
        st_matrix("`distances'",_wd_D)
        st_matrix("`groups'",_wd_groups)
        st_matrix("`row_ids'",select((1::st_nobs()),st_data(.,"`touse'")))
        st_matrix("`neff'",1:/colsum(_wd_W:^2))
    }
    if "`generate'"!="" {
        quietly gen double `generate'=.
        mata: st_store(.,"`generate'","`touse'",_wd_W)
    }
    if inlist("`distance'","mahalanobis","l1_scaled","l2_scaled") {
        mata: st_matrix("`geometry'",_wd_G)
        return matrix geometry=`geometry'
    }
    return matrix weights=`weights'
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
