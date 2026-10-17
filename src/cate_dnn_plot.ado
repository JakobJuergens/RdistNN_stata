*! version 0.2.0  2026-10-17  Jakob Juergens
*! Stata-native postestimation plots for cate_dnn
program define cate_dnn_plot, rclass
    version 17
    capture program list _rdistnn_plot_dispatch
    if _rc {
        findfile "_rdistnn_plot_dispatch.ado"
        quietly run "`r(fn)'"
    }
    _rdistnn_plot_dispatch cate `0'
    return add
end
