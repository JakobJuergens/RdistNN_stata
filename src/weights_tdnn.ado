*! version 0.1.0  2026-10-18
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
    capture confirm integer number `dim'
    if _rc error 198
    if `s1'<1 | `s1'>=`s2' | `s2'>`n' | `dim'<1 {
        di as error "Require 1 <= s1 < s2 <= n and positive dim()."
        exit 198
    }
    if "`method'"=="" local method "exact"
    if !inlist("`method'","exact","asymp") error 198
    if "`distance'"=="" local distance "l2"
    if "`scalerule'"=="" local scalerule "sd"
    if ("`x0'"=="")==("`x0matrix'"=="") {
        di as error "Specify exactly one of x0() or x0matrix()."
        exit 198
    }
    tempname targets W R DX table diag geom row_ids
    if "`x0matrix'"!="" {
        capture confirm matrix `x0matrix'
        if _rc error 198
        matrix `targets'=`x0matrix'
    }
    else {
        local row : subinstr local x0 " " ",", all
        matrix `targets'=(`row')
    }
    local m=rowsof(`targets')
    if `m'<1 | colsof(`targets')!=`p' {
        di as error "Targets must have one column per localization covariate."
        exit 198
    }
    if "`generate'"!="" {
        confirm new variable `generate'
        if `m'!=1 {
            di as error "generate() requires a single target; use r(weights) for batches."
            exit 198
        }
    }
    local normalize=("`normalize'"=="")
    foreach file in _rdistnn_mata.ado _rdistnn_tdnn_mata.ado {
        findfile "`file'"
        quietly run "`r(fn)'"
    }
    mata {
        _td_X=st_data(.,tokens("`varlist'"),"`touse'")
        _td_X0=st_matrix("`targets'")
        if(hasmissing(_td_X0)) exit(error(198))
        _td_geom=rdistnn_fit_geometry(_td_X,"`distance'","`scalerule'",strtoreal(tokens("`scalevec'")),`scaleeps')
        _td_DX=rdistnn_distance_matrix(_td_X,_td_X0,"`distance'",_td_geom)
        _td_R=rdistnn_ranks_from_dmat(_td_DX,_td_X)
        _td_table=rdistnn_tdnn_table(`n',`s1',`s2',`dim',"`method'",`normalize')
        _td_W=rdistnn_tdnn_weights(_td_R,`s1',`s2',`dim',"`method'",`normalize')
        st_matrix("`geom'",_td_geom)
        st_matrix("`row_ids'",select((1::st_nobs()),st_data(.,"`touse'")))
        st_matrix("`W'",_td_W); st_matrix("`R'",_td_R); st_matrix("`DX'",_td_DX)
        st_matrix("`table'",_td_table)
        st_matrix("`diag'",rdistnn_tdnn_diagnostics(_td_W,_td_DX))
    }
    if "`generate'"!="" {
        quietly gen double `generate'=.
        mata: st_store(.,"`generate'","`touse'",st_matrix("`W'"))
    }
    matrix colnames `diag'=signed_sum positive_mass negative_mass abs_mass ///
        inverse_l2_concentration max_abs_weight n_nonzero_weight abs_mean_distance
    matrix colnames `row_ids'=observation_id
    return matrix weights=`W'
    return matrix ranks=`R'
    return matrix distances=`DX'
    return matrix rank_weights=`table'
    return matrix targets=`targets'
    return matrix diagnostics=`diag'
    return matrix row_ids=`row_ids'
    capture confirm matrix `geom'
    if !_rc return matrix geometry=`geom'
    return scalar n=`n'
    return scalar m=`m'
    return scalar p=`p'
    return scalar s1=`s1'
    return scalar s2=`s2'
    return scalar dim=`dim'
    return scalar normalize=`normalize'
    return local method "`method'"
    return local distance "`distance'"
    return local tie_rule "distance_then_covariates_then_input_position"
end
