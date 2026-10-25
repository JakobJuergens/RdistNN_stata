*! version 0.1.5  2026-10-25  Jakob Juergens
*! Stata-native postestimation plots for clate_dnn
program define clate_dnn_plot, rclass
    version 17
    capture program list _rdistnn_plot_dispatch
    local dispatch_rc=_rc
    if !inlist(`dispatch_rc',0,111) exit `dispatch_rc'
    if `dispatch_rc'==111 {
        findfile "_rdistnn_plot_dispatch.ado"
        quietly run "`r(fn)'"
    }
    _rdistnn_plot_dispatch clate `0'
    return add
end
