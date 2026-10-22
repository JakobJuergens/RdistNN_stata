*! version 0.1.0  2026-10-18
program define np_tdnn_plot, rclass
    version 17
    syntax [, TYPE(string) STYLE(string) XINDEX(integer 0) ///
        TARGETS(numlist integer >=1) Level(cilevel) NAME(name) REPLACE ///
        NODRAW SAVING(string asis)]
    if "`e(cmd)'"!="np_tdnn" error 301
    if "`type'"=="" local type "estimates"
    if !inlist("`type'","estimates","weights","diagnostics") {
        di as error "type() must be estimates, weights or diagnostics for signed TDNN."
        exit 198
    }
    if "`type'"!="estimates" & ("`style'"!="" | `xindex'!=0) error 198
    if `xindex'<0 | `xindex'>e(p) error 198
    foreach j of local targets {
        if `j'>e(m) error 198
    }
    if "`name'"=="" local name "np_tdnn_`type'"
    local nameopt "name(`name')"
    if "`replace'"!="" local nameopt "name(`name', replace)"
    tempname data
    if "`type'"=="estimates" {
        capture program list _cate_dnn_plot_estimates
        if _rc {
            findfile "_rdistnn_plot_helpers.ado"
            quietly run "`r(fn)'"
        }
        findfile "_rdistnn_mata.ado"
        quietly run "`r(fn)'"
        local targetopt ""
        if "`targets'"!="" local targetopt "targets(`targets')"
        local styleopt ""
        if "`style'"!="" local styleopt "style(`style')"
        _cate_dnn_plot_estimates, name(`name') level(`level') xindex(`xindex') ///
            `targetopt' `styleopt' `replace' `nodraw'
        matrix `data'=r(plot_data)
    }
    else {
        if "`type'"=="weights" {
            foreach file in _rdistnn_mata.ado _rdistnn_tdnn_mata.ado _rdistnn_model_mata.ado {
                findfile "`file'"
                quietly run "`r(fn)'"
            }
            if e(model_version)!=2 error 498
            if "`targets'"=="" {
                numlist "1/`=e(m)'"
                local targets "`r(numlist)'"
            }
            noisily mata: _td_plot=rdistnn_model_weights(strtoreal(tokens("`targets'")))
            mata: _td_plot=_td_plot[order(_td_plot,(2,3)),]
            mata: st_matrix("`data'",_td_plot[,(3,5,2)])
            matrix colnames `data'=rank weight target_id
        }
        else {
            local targetopt ""
            if "`targets'"!="" local targetopt "targets(`targets')"
            quietly np_tdnn_diagnostics, `targetopt'
            matrix `data'=r(compact)
        }
        preserve
        clear
        quietly svmat double `data', names(col)
        if "`type'"=="weights" {
            twoway line weight rank, by(target_id, note("")) yline(0) xtitle("Distance rank") ///
                ytitle("Signed weight") title("TDNN rank weights") `nameopt' `nodraw'
        }
        else {
            twoway scatter abs_mass target_id, xtitle("Target") ///
                ytitle("Sum of absolute weights") title("TDNN weight concentration") ///
                yscale(range(0)) ylabel(, format(%5.2f)) `nameopt' `nodraw'
        }
        restore
    }
    if `"`saving'"'!="" {
        capture program list _rdistnn_plot_export
        if _rc {
            findfile "_rdistnn_plot_dispatch.ado"
            quietly run "`r(fn)'"
        }
        quietly _rdistnn_plot_export, name(`name') saving(`saving') `replace'
    }
    local n_rows=rowsof(`data')
    return matrix plot_data=`data'
    return scalar n_rows=`n_rows'
    return local type "`type'"
    return local name "`name'"
    return local saved `saving'
end
