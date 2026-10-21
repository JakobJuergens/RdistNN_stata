*! version 0.1.0  2026-10-18
program define np_tdnn, eclass
    version 17
    if replay() {
        if "`e(cmd)'"!="np_tdnn" error 301
        _np_tdnn_Display `0'
        exit
    }
    syntax varname(numeric) [if] [in], [X(varlist numeric) XCOLS(varlist numeric) ///
        X0(numlist) X0MATRIX(name) S1(integer 1) S2(integer 2) DIM(string) ///
        METHOD(string) noNORMALIZE DISTance(string) SCALERule(string) ///
        SCALEVec(numlist) SCALEEPS(real 1e-12) noSE Level(cilevel)]
    if "`x'"!="" & "`xcols'"!="" error 198
    local xvars `x' `xcols'
    if "`xvars'"=="" {
        di as error "xcols() is required (legacy alias: x())."
        exit 198
    }
    local outcome `varlist'
    marksample touse
    markout `touse' `xvars'
    quietly count if `touse'
    local n=r(N)
    if `n'==0 error 2000
    if "`se'"=="" & `s2'>=`n' {
        di as error "Delete-one jackknife requires s2 < n; reduce s2() or specify nose."
        exit 198
    }
    local p : word count `xvars'
    if "`dim'"=="" local dim=`p'
    capture confirm integer number `dim'
    if _rc error 198
    if "`method'"=="" local method "exact"
    if "`distance'"=="" local distance "l2"
    if "`scalerule'"=="" local scalerule "sd"
    local targetopt ""
    if "`x0'"!="" local targetopt "x0(`x0')"
    if "`x0matrix'"!="" local targetopt "`targetopt' x0matrix(`x0matrix')"
    local scaleopt ""
    if "`scalevec'"!="" local scaleopt "scalevec(`scalevec')"
    local normopt "`normalize'"
    local normalized=("`normalize'"=="")
    quietly weights_tdnn `xvars' if `touse', `targetopt' s1(`s1') s2(`s2') ///
        dim(`dim') method(`method') distance(`distance') scalerule(`scalerule') ///
        `scaleopt' scaleeps(`scaleeps') `normopt'
    tempname W R DX T diag geom b V signals concentration
    matrix `W'=r(weights)
    matrix `R'=r(ranks)
    matrix `DX'=r(distances)
    matrix `T'=r(targets)
    matrix `diag'=r(diagnostics)
    local m=rowsof(`T')
    local geopt ""
    capture matrix `geom'=r(geometry)
    if !_rc local geopt "geometry(`geom')"
    mata {
        _td_y=st_data(.,"`outcome'","`touse'")
        st_matrix("`signals'",_td_y)
        st_matrix("`b'",rdistnn_weighted_signals(_td_y,st_matrix("`W'")))
    }
    local labels ""
    forvalues j=1/`m' {
        local labels "`labels' target`j'"
    }
    matrix colnames `b'=`labels'
    matrix rownames `T'=`labels'
    matrix colnames `T'=`xvars'
    matrix rownames `diag'=`labels'
    mata: st_matrix("`concentration'",st_matrix("`diag'")[,5]')
    matrix colnames `concentration'=`labels'
    if "`se'"=="" {
        mata: st_matrix("`V'",rdistnn_tdnn_cov(st_matrix("`R'"),st_matrix("`signals'"),`s1',`s2',`dim',"`method'",`normalized'))
        matrix rownames `V'=`labels'
        matrix colnames `V'=`labels'
        ereturn post `b' `V', esample(`touse')
    }
    else ereturn post `b', esample(`touse')
    ereturn scalar N=`n'
    ereturn scalar m=`m'
    ereturn scalar p=`p'
    ereturn scalar kx=`p'
    ereturn scalar k=1
    ereturn scalar s=`s2'
    ereturn scalar s1=`s1'
    ereturn scalar s2=`s2'
    ereturn scalar dim=`dim'
    ereturn scalar normalize=`normalized'
    mata: st_local("concentration_min",strofreal(min(st_matrix("`concentration'")),"%21.0g")); st_local("concentration_max",strofreal(max(st_matrix("`concentration'")),"%21.0g"))
    // Compatibility names retain signed squared-weight concentration only.
    ereturn scalar neff_min=`concentration_min'
    ereturn scalar neff_max=`concentration_max'
    if `m'==1 ereturn scalar neff=el(`concentration',1,1)
    ereturn matrix neff_targets=`concentration', copy
    ereturn scalar scaleeps=`scaleeps'
    ereturn local scalerule "`scalerule'"
    ereturn local method "`method'"
    ereturn local distance "`distance'"
    ereturn local xvars "`xvars'"
    ereturn local x_cols "`xvars'"
    ereturn local outcome "`outcome'"
    ereturn local signal "`outcome'"
    ereturn local x0 "`x0'"
    ereturn local x0matrix "`x0matrix'"
    ereturn local target_source=cond("`x0matrix'"!="","X0","x0")
    ereturn matrix targets=`T', copy
    ereturn matrix x0_matrix=`T', copy
    ereturn matrix diagnostics=`diag', copy
    ereturn matrix target_diagnostics=`diag', copy
    ereturn local se_method=cond("`se'"=="","signed_rank_shift","none")
    ereturn local vce=cond("`se'"=="","signed_rank_shift","none")
    ereturn local vcetype=cond("`se'"=="","Jackknife","")
    ereturn local depvar "`outcome'"
    ereturn local jackknife=cond("`se'"=="","fast","none")
    ereturn local tuning_mode "default"
    ereturn local cmd "np_tdnn"
    ereturn local cmdline "np_tdnn `0'"
    local vars `xvars' `outcome'
    local vars : list uniq vars
    quietly _rdistnn_estimation_state `vars'
    quietly _rdistnn_model_store `vars', signals(`signals') ranks(`R') distances(`DX') `geopt'
    _np_tdnn_Display, level(`level')
end

program define _np_tdnn_Display, rclass
    version 17
    _rdistnn_report `0'
    return add
end
