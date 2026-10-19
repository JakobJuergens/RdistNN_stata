*! version 0.1.0  2026-10-18
program define np_tdnn_diagnostics, rclass
    version 17
    syntax [, TARGETS(numlist integer >=1) DETAIL SAVING(string asis) REPLACE]
    if "`e(cmd)'"!="np_tdnn" error 301
    local selected "`targets'"
    if "`selected'"=="" {
        forvalues j=1/`=e(m)' {
            local selected "`selected' `j'"
        }
    }
    local selected : list uniq selected
    foreach j of local selected {
        if `j'>e(m) error 198
    }
    tempname out compact
    mata: st_matrix("`out'",(strtoreal(tokens("`selected'"))',st_matrix("e(targets)")[strtoreal(tokens("`selected'")),],st_matrix("e(diagnostics)")[strtoreal(tokens("`selected'")),]))
    local xnames ""
    forvalues j=1/`=e(p)' {
        local xnames "`xnames' x0_`j'"
    }
    matrix colnames `out'=target_id `xnames' signed_sum positive_mass ///
        negative_mass abs_mass inverse_l2_concentration max_abs_weight ///
        n_nonzero_weight abs_mean_distance
    mata: st_matrix("`compact'",st_matrix("`out'")[,(1,`=e(p)+2',`=e(p)+4',`=e(p)+5',`=e(p)+7')])
    matrix colnames `compact'=target_id signed_sum negative_mass abs_mass max_abs_weight
    di as text "Signed weights: positive/negative mass and absolute concentration; no probability quantiles."
    if "`detail'"!="" matrix list `out', noheader format(%12.6g)
    else matrix list `compact', noheader format(%12.6g)
    if `"`saving'"'!="" _rdistnn_diagnostics_save, matrix(`out') saving(`saving') `replace'
    local n_rows=rowsof(`out')
    return matrix diagnostics=`out'
    return matrix compact=`compact'
    return scalar n_rows=`n_rows'
    return local selected "`selected'"
    return local scope "targets"
    return local saved `saving'
end
