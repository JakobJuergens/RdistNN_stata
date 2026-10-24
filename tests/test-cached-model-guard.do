version 17
do "tests/test_framework.do"
do "tests/helpers/e-state.do"
test_section "all cached consumers validate distance-group snapshots"
clear
set obs 8
gen double guard_x=_n-4.5
gen double guard_y=_n^2
matrix GUARD_TARGETS=(0\1)
quietly np_dnn guard_y, xcols(guard_x) x0matrix(GUARD_TARGETS) s(3)
capture program drop _distance_guard_corrupt
program define _distance_guard_corrupt, eclass
    args fault
    tempname groups notices
    if "`fault'"!="groups" {
        matrix `notices'=e(target_notices)
        if "`fault'"=="tie" matrix `notices'[1,1]=1-el(`notices',1,1)
        if "`fault'"=="weak" matrix `notices'[1,2]=1-el(`notices',1,2)
        if "`fault'"=="missing" matrix `notices'=J(1,1,.)
        ereturn matrix target_notices=`notices', copy
        // Keep the aggregate consistent with the false flag, so validation
        // must check the actual geometry/diagnostic facts as well as totals.
        if "`fault'"=="tie" ereturn scalar n_tied_targets=e(n_tied_targets)+2*el(`notices',1,1)-1
        if "`fault'"=="weak" ereturn scalar n_weak_targets=e(n_weak_targets)+2*el(`notices',1,2)-1
        if "`fault'"=="count" ereturn scalar n_tied_targets=e(n_tied_targets)+1
        if "`fault'"=="weakcount" ereturn scalar n_weak_targets=.
        if "`fault'"=="schema" ereturn scalar model_version=3
        exit
    }
    matrix `groups'=e(model_groups)
    matrix `groups'[1,1]=99
    ereturn matrix model_groups=`groups', copy
end
estimates store DISTANCE_GUARD_GOOD
foreach fault in groups tie weak missing count weakcount schema {
quietly estimates restore DISTANCE_GUARD_GOOD
_distance_guard_corrupt `fault'
estimates store DISTANCE_GUARD_BAD
preserve
clear
mata: mata clear
_rdistnn_test_e_snapshot, name(DISTANCE_GUARD_STATE)
foreach route in estimates vcov resultdiag diagnostics estimateplot diagplot {
    quietly estimates restore DISTANCE_GUARD_BAD
    if "`route'"=="estimates" capture quietly rdistnn_results
    if "`route'"=="vcov" capture quietly rdistnn_results, type(vcov)
    if "`route'"=="resultdiag" capture quietly rdistnn_results, type(diagnostics)
    if "`route'"=="diagnostics" capture quietly np_dnn_diagnostics
    if "`route'"=="estimateplot" capture quietly np_dnn_plot, type(estimates) name(distance_guard) replace nodraw
    if "`route'"=="diagplot" capture quietly np_dnn_plot, type(diagnostics) name(distance_guard) replace nodraw
    local rc=_rc
    test_assert "cold `route' rejects inconsistent cached `fault' state" `rc'==498
    _rdistnn_test_e_compare, snapshot(DISTANCE_GUARD_STATE)
    local unchanged=r(equal)
    test_assert "cold `route' `fault' preserves all fitted state and empty data" `unchanged' & _N==0
}
_rdistnn_test_e_snapshot_drop, name(DISTANCE_GUARD_STATE)
restore
estimates drop DISTANCE_GUARD_BAD
}
estimates drop DISTANCE_GUARD_GOOD
program drop _distance_guard_corrupt
