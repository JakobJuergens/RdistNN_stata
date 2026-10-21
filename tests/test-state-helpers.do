version 17
do "tests/helpers/fitted-data.do"
do "tests/helpers/e-state.do"
test_section "shared complete e() comparator owns estimation state"
quietly np_dnn fr_y, xcols(fr_x) x0(0) s(3)
estimates store STATE_REFERENCE
tempfile serialized
quietly estimates save "`serialized'", replace
quietly estimates use "`serialized'"
estimates store STATE_CURRENT
_rdistnn_test_e_compare STATE_REFERENCE, numeric
local equal=r(equal)
test_assert "numeric e() comparison ignores field enumeration order after serialization" `equal'
_rdistnn_test_e_compare STATE_CURRENT
local equal=r(equal)
test_assert "complete e() comparison preserves the active result" `equal'
capture program drop _state_test_mutate
program define _state_test_mutate, eclass
    syntax, KIND(string)
    if "`kind'"=="macro" ereturn local extra "descriptive metadata"
    if "`kind'"=="scalar" ereturn scalar extra=42
    if "`kind'"=="sample" {
        tempvar sample
        quietly gen byte `sample'=e(sample) & _n>1
        ereturn repost, esample(`sample')
    }
    if "`kind'"=="stripe" {
        tempname b
        matrix `b'=e(estimates)
        matrix colnames `b'=different_stripe
        ereturn matrix estimates=`b'
    }
end
foreach kind in macro scalar sample stripe {
    quietly estimates restore STATE_REFERENCE
    _state_test_mutate, kind(`kind')
    _rdistnn_test_e_compare STATE_REFERENCE
    local equal=r(equal)
    test_assert "complete e() comparison detects changed `kind'" !`equal'
    estimates store STATE_MUTATED
    _rdistnn_test_e_compare STATE_MUTATED
    local equal=r(equal)
    test_assert "comparison preserves changed `kind' state on return" `equal'
    if "`kind'"!="stripe" {
        _rdistnn_test_e_compare STATE_REFERENCE, numeric
        local equal=r(equal)
        test_assert "numeric mode compares matrices independently of descriptive `kind'" `equal'
    }
    estimates drop STATE_MUTATED
}
quietly estimates restore STATE_REFERENCE
gen str8 caller_text=cond(mod(_n,2),"odd","even")
gen byte caller_class=mod(_n,2)
label data "Complete e helper caller data"
label variable fr_x "caller covariate"
format fr_y %12.4f
label define STATE_CLASS 0 "control" 1 "treated"
label values caller_class STATE_CLASS
do "tests/helpers/caller-data.do"
_rdistnn_test_data_snapshot
tempfile caller
quietly save "`caller'", replace
_rdistnn_test_e_compare STATE_REFERENCE
capture cf _all using "`caller'"
local rc=_rc
test_assert "complete e() comparator leaves caller data and sample marker ownership intact" `rc'==0
_rdistnn_test_e_snapshot, name(STATE_INERT)
capture cf _all using "`caller'"
local rc=_rc
test_assert "inert complete e() snapshot creates no dataset marker variables" `rc'==0
_rdistnn_test_data_compare, data("`caller'")
local equal=r(equal)
test_assert "inert capture preserves caller values, variable order, storage, labels and formats" `equal'
_rdistnn_test_e_compare STATE_REFERENCE
local equal=r(equal)
test_assert "inert capture preserves every active e() field and e(sample)" `equal'
_rdistnn_test_e_compare, snapshot(STATE_INERT)
local equal=r(equal)
test_assert "inert snapshot compares every e() field and the sample exactly" `equal'
capture noisily _rdistnn_test_e_compare STATE_NO_SUCH_RESULT
local rc=_rc
test_assert "missing named reference returns a failure" `rc'!=0
_rdistnn_test_e_compare, snapshot(STATE_INERT)
local equal=r(equal)
_rdistnn_test_data_compare, data("`caller'")
local same_data=r(equal)
test_assert "failed named comparison preserves complete caller e() and data" `equal' & `same_data'
_rdistnn_test_e_snapshot_drop, name(STATE_INERT)
estimates drop STATE_CURRENT STATE_REFERENCE
program drop _state_test_mutate
