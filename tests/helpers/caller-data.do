version 17
// Caller-data snapshot metadata is independent of estimation state comparison.
capture program drop _rdistnn_test_data_snapshot
program define _rdistnn_test_data_snapshot
    global PMPE_datalabel : data label
    unab variables : _all
    global PMPE_variables "`variables'"
    local i 0
    foreach variable of local variables {
        local ++i
        global PMPE_type`i' : type `variable'
        global PMPE_format`i' : format `variable'
        global PMPE_label`i' : variable label `variable'
        global PMPE_valuelabel`i' : value label `variable'
    }
end
capture program drop _rdistnn_test_data_compare
program define _rdistnn_test_data_compare, rclass
    syntax, DATA(string)
    local ok 1
    capture cf _all using `"`data'"'
    if _rc local ok 0
    local datalabel : data label
    if `"`datalabel'"'!=`"$PMPE_datalabel"' local ok 0
    unab variables : _all
    if "`variables'"!="$PMPE_variables" local ok 0
    local i 0
    foreach variable of local variables {
        local ++i
        foreach property in type format label valuelabel {
            local original : copy global PMPE_`property'`i'
            if "`property'"=="type" local current : type `variable'
            if "`property'"=="format" local current : format `variable'
            if "`property'"=="label" local current : variable label `variable'
            if "`property'"=="valuelabel" local current : value label `variable'
            if `"`current'"'!=`"`original'"' local ok 0
        }
    }
    return scalar equal=`ok'
end
