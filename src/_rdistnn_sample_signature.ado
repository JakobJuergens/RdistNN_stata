*! version 0.1.5  2026-10-25
// Pair exact input values with row indices: ordinary data signatures ignore order.
program define _rdistnn_sample_signature, rclass
    version 17
    syntax varlist(numeric) [if] [in]
    marksample touse, novarlist
    tempvar encoded
    quietly gen str64 `encoded' = "" if `touse'
    local signature ""
    foreach v of local varlist {
        quietly replace `encoded' = string(_n,"%12.0f") + ":" + ///
            string(`v',"%21x") if `touse'
        quietly _datasignature `encoded' if `touse', nonames
        local signature "`signature'|`r(datasignature)'"
    }
    return local signature "`signature'"
end
