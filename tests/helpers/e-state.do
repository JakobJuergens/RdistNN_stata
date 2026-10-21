version 17
// Complete e() snapshots do not create estimation stores or dataset markers.
capture program drop _rdistnn_test_e_snapshot
program define _rdistnn_test_e_snapshot
    syntax, NAME(name)
    global `name'_matrices : e(matrices)
    global `name'_scalars : e(scalars)
    global `name'_macros : e(macros)
    local fields : e(matrices)
    local i 0
    foreach field of local fields {
        local ++i
        matrix `name'_m`i'=e(`field')
    }
    local fields : e(scalars)
    local i 0
    foreach field of local fields {
        local ++i
        scalar `name'_s`i'=e(`field')
    }
    local fields : e(macros)
    local i 0
    foreach field of local fields {
        local ++i
        global `name'_g`i' `"`e(`field')'"'
    }
    global `name'_nobs=_N
    if _N>0 {
        tempvar sample
        quietly gen byte `sample'=e(sample)
        quietly mkmat `sample', matrix(`name'_sample)
    }
end

capture program drop _rdistnn_test_e_snapshot_drop
program define _rdistnn_test_e_snapshot_drop
    syntax, NAME(name)
    local fields : copy global `name'_matrices
    local n : word count `fields'
    if `n'>0 {
        forvalues i=1/`n' {
            capture matrix drop `name'_m`i'
        }
    }
    local fields : copy global `name'_scalars
    local n : word count `fields'
    if `n'>0 {
        forvalues i=1/`n' {
            capture scalar drop `name'_s`i'
        }
    }
    capture matrix drop `name'_sample
    macro drop `name'_*
end

capture program drop _rdistnn_test_e_snapshot_compare
program define _rdistnn_test_e_snapshot_compare, rclass
    syntax, NAME(name) [NUMERIC]
    local ok 1
    local matrices : copy global `name'_matrices
    local current_matrices : e(matrices)
    // Serialization may reverse field enumeration without changing state.
    local extra : list current_matrices - matrices
    local absent : list matrices - current_matrices
    if "`extra'`absent'"!="" {
        di as error "Matrix fields differ: extra `extra'; absent `absent'"
        local ok 0
    }
    local i 0
    foreach field of local matrices {
        local ++i
        if !`: list field in current_matrices' local ok 0
        else {
            tempname current
            matrix `current'=e(`field')
            if rowsof(`current')!=rowsof(`name'_m`i') | colsof(`current')!=colsof(`name'_m`i') local ok 0
            else {
                mata: st_numscalar("r(_test_e_equal)",all(st_matrix("`current'"):==st_matrix("`name'_m`i'")) & all(st_matrixrowstripe("`current'"):==st_matrixrowstripe("`name'_m`i'")) & all(st_matrixcolstripe("`current'"):==st_matrixcolstripe("`name'_m`i'")))
                if !r(_test_e_equal) {
                    di as error "Numeric values or stripes differ: e(`field')"
                    local ok 0
                }
            }
        }
    }
    if "`numeric'"=="" {
        local macros : copy global `name'_macros
        local scalars : copy global `name'_scalars
        local current_macros : e(macros)
        local current_scalars : e(scalars)
        local extra : list current_macros - macros
        local absent : list macros - current_macros
        local extra_s : list current_scalars - scalars
        local absent_s : list scalars - current_scalars
        if "`extra'`absent'`extra_s'`absent_s'"!="" {
            di as error "Macro/scalar fields differ"
            local ok 0
        }
        local i 0
        foreach field of local macros {
            local ++i
            local original : copy global `name'_g`i'
            if `"`e(`field')'"'!=`"`original'"' {
                di as error "Macro differs: e(`field')"
                local ok 0
            }
        }
        local i 0
        foreach field of local scalars {
            local ++i
            if e(`field')!=scalar(`name'_s`i') {
                di as error "Scalar differs: e(`field')"
                local ok 0
            }
        }
        local nobs : copy global `name'_nobs
        if _N!=`nobs' local ok 0
        else if _N>0 {
            tempvar sample
            tempname current
            quietly gen byte `sample'=e(sample)
            quietly mkmat `sample', matrix(`current')
            mata: st_numscalar("r(_test_e_equal)",all(st_matrix("`current'"):==st_matrix("`name'_sample")))
            if !r(_test_e_equal) {
                di as error "e(sample) differs"
                local ok 0
            }
        }
    }
    return scalar equal=`ok'
end

// Compare with an existing named estimate or an inert complete-e snapshot.
// Full mode includes every macro, scalar and e(sample); numeric mode is matrices.
capture program drop _rdistnn_test_e_compare
program define _rdistnn_test_e_compare, rclass
    syntax [name(name=reference)] [, NUMERIC SNAPSHOT(name)]
    if "`snapshot'"!="" {
        if "`reference'"!="" error 198
        _rdistnn_test_e_snapshot_compare, name(`snapshot') `numeric'
        return add
        exit
    }
    if "`reference'"=="" error 198
    tempname held snapshot_name
    // Leading underscores are reserved for Stata local macros, not globals.
    local snapshot_name "RDN`snapshot_name'"
    _estimates hold `held', copy restore
    capture noisily {
        _rdistnn_test_e_snapshot, name(`snapshot_name')
        quietly estimates restore `reference'
        _rdistnn_test_e_snapshot_compare, name(`snapshot_name') `numeric'
        local ok=r(equal)
    }
    local rc=_rc
    _estimates unhold `held'
    _rdistnn_test_e_snapshot_drop, name(`snapshot_name')
    if `rc' exit `rc'
    return scalar equal=`ok'
end
