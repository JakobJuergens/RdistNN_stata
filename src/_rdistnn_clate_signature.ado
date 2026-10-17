*! version 0.1.0  2026-10-17
// Unlike Stata's ordinary data signature, this records input row order.
// Each checksum pairs a variable's exact hexadecimal values with row indices.
program define _rdistnn_clate_signature, rclass
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
