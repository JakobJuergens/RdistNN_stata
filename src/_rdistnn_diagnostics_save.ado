*! version 0.1.0  2026-10-17  Jakob Juergens
*! Save diagnostic tables while preserving the caller's data on errors

program define _rdistnn_diagnostics_save
    version 17
    syntax, MATRIX(name) [SAVING(string asis) REPLACE]

    if `"`saving'"' == "" exit
    preserve
    clear
    capture noisily {
        quietly svmat double `matrix', names(col)
        if "`replace'" != "" save `saving', replace
        else save `saving'
    }
    local rc = _rc
    restore
    if `rc' exit `rc'
end
