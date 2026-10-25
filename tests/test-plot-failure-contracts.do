// Plot-bound probes and required dashboard panels propagate genuine failures.
// Use public plot consumers with actual helper source, independent caller-state
// snapshots and deterministic data; restore the replaced helper providers.
version 17
do "tests/test_framework.do"
clear
set obs 12
gen double x=_n
gen double y=2*x
test_section "plot bound checks propagate native failures"
quietly np_dnn y, xcols(x) x0(0) tuningmode(grid) sgrid(3 4) sdefault(4)
do "tests/helpers/reporting-instrument.do"
do "tests/helpers/e-state.do"
do "tests/helpers/caller-data.do"
capture program drop MR_bounds_assert
program MR_bounds_assert
    if "$MR_BOUND_RC"!="0" exit $MR_BOUND_RC
    assert `0'
end
capture program drop MR_bound_fixture
program MR_bound_fixture
    tempfile injected caller_data
    quietly save "`caller_data'", replace
    _rdistnn_test_data_snapshot
    _rdistnn_test_e_snapshot, name(MR_BOUND_STATE)
    mata: rdistnn_test_instrument("src/_rdistnn_plot_helpers.ado",st_local("injected"),"capture assert missing(ci_lower)","capture MR_bounds_assert missing(ci_lower)")
    local helpers _cate_dnn_plot_estimates _cate_dnn_plot_diagnostics _cate_dnn_plot_tuning _cate_dnn_plot_weights _cate_dnn_plot_overlap _cate_dnn_plot_influence _cate_dnn_plot_dashboard _np_dnn_plot_outcome _np_dnn_plot_dashboard _rdistnn_plot_wald_bounds
    foreach helper of local helpers {
        capture program drop `helper'
    }
    quietly run "`injected'"
    capture noisily {
        foreach code in 1 909 498 {
            global MR_BOUND_RC `code'
            foreach style in profile forest pointrange tuning {
                local opts "type(estimates) style(`style') xindex(1)"
                if "`style'"=="tuning" local opts "type(tuning)"
                capture quietly np_dnn_plot, `opts' name(mr_bounds) replace nodraw
                local result=_rc
                test_assert "`style': bounds detection preserves native`code'" `result'==`code'
                _rdistnn_test_e_compare, snapshot(MR_BOUND_STATE)
                test_assert "`style': bounds failure preserves complete estimates" r(equal)
                _rdistnn_test_data_compare, data("`caller_data'")
                test_assert "`style': bounds failure preserves caller data" r(equal)
            }
        }
    }
    local fixture_rc=_rc
    foreach helper of local helpers {
        capture program drop `helper'
    }
    quietly run "src/_rdistnn_plot_helpers.ado"
    _rdistnn_test_e_snapshot_drop, name(MR_BOUND_STATE)
    macro drop MR_BOUND_RC
    if `fixture_rc' exit `fixture_rc'
end
MR_bound_fixture
program drop MR_bound_fixture MR_bounds_assert

clear
set obs 24
gen double x=mod(_n,7)
gen double y=_n
quietly cate_dnn y, xcols(x) x0(0) pseudooutcome s(4)
// Capability checks already selected this required weights panel. Its resource
// failure must not turn into a dashboard that silently drops the panel.
capture program drop _cate_dnn_plot_weights
program _cate_dnn_plot_weights
    exit 909
end
capture quietly cate_dnn_plot, type(dashboard) name(mr_fault) replace nodraw
local MR_PANEL_RC=_rc
test_assert "dashboard propagates an unexpected native/resource panel failure" `MR_PANEL_RC'==909
capture program drop _cate_dnn_plot_weights
// Reload actual helpers for later suites without redefining existing programs.
program drop _cate_dnn_plot_estimates _cate_dnn_plot_diagnostics _cate_dnn_plot_tuning _cate_dnn_plot_overlap _cate_dnn_plot_influence _cate_dnn_plot_dashboard _np_dnn_plot_outcome _np_dnn_plot_dashboard _rdistnn_plot_wald_bounds
quietly run "src/_rdistnn_plot_helpers.ado"
