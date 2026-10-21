version 17
do "tests/helpers/e-state.do"
do "tests/helpers/caller-data.do"
capture program drop _pmpe_snapshot
program define _pmpe_snapshot
    _rdistnn_test_data_snapshot
    _rdistnn_test_e_snapshot, name(PMPE_ESTATE)
end
capture program drop _pmpe_unchanged
program define _pmpe_unchanged, rclass
    syntax, DATA(string)
    _rdistnn_test_data_compare, data(`"`data'"')
    local ok=r(equal)
    _rdistnn_test_e_compare, snapshot(PMPE_ESTATE)
    return scalar ok=`ok' & r(equal)
end
