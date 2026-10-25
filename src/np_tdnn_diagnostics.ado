*! version 0.1.5  2026-10-25
program define np_tdnn_diagnostics, rclass
    version 17
    syntax [, TARGETS(numlist integer >=1) DETAIL SAVING(string asis) REPLACE]
    local targetopt ""
    if "`targets'"!="" local targetopt "targets(`targets')"
    local savingopt ""
    if `"`saving'"'!="" local savingopt `"saving(`saving')"'
    _rdistnn_diagnostics, command(np_tdnn) `targetopt' `detail' `savingopt' `replace'
    return add
end
