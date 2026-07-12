*! version 0.1.0  2026-07-12  Jakob Juergens
*! Stata-native postestimation plots for np_dnn

program define np_dnn_plot, rclass
    version 17

    syntax [, TYPE(string) VIEW(string) STYLE(string) XINDEX(integer 0) ///
        TARGETS(numlist integer >=1) CANDIDATES(numlist integer >=1) ///
        ALLCandidates Level(cilevel) LOGWeights NAME(name) REPLACE ///
        NODRAW SAVING(string asis)]

    if "`e(cmd)'" != "np_dnn" {
        di as error "np_dnn_plot requires active np_dnn estimates."
        exit 301
    }
    if "`type'" == "" local type "estimates"
    local type = lower("`type'")
    if !inlist("`type'", "estimates", "diagnostics", "tuning", ///
        "weights", "outcome", "dashboard") {
        di as error "type() must be estimates, diagnostics, tuning, weights, outcome, or dashboard."
        exit 198
    }
    if "`candidates'" != "" & "`allcandidates'" != "" {
        di as error "Specify either candidates() or allcandidates, not both."
        exit 198
    }
    if `xindex' < 0 | `xindex' > e(p) {
        di as error "xindex() must be 0 or a target-coordinate index in 1..`=e(p)'."
        exit 198
    }
    foreach j of local targets {
        if `j' > e(m) {
            di as error "targets() entries must be in 1..`=e(m)'."
            exit 198
        }
    }

    capture program list _cate_dnn_plot_estimates
    if _rc {
        findfile "cate_dnn_plot.ado"
        quietly run "`r(fn)'"
    }
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    if "`name'" == "" local name "np_dnn_`type'"
    local common "name(`name')"
    if "`replace'" != "" local common "`common' replace"
    if "`nodraw'" != "" local common "`common' nodraw"
    local candidate_opts ""
    if "`candidates'" != "" local candidate_opts "candidates(`candidates')"
    if "`allcandidates'" != "" local candidate_opts "allcandidates"
    local target_opts ""
    if "`targets'" != "" local target_opts "targets(`targets')"
    local style_opt ""
    if "`style'" != "" local style_opt "style(`style')"
    local view_opt ""
    if "`view'" != "" local view_opt "view(`view')"

    if "`type'" == "estimates" {
        _cate_dnn_plot_estimates, `common' level(`level') ///
            `style_opt' xindex(`xindex') `target_opts' `candidate_opts'
    }
    else if "`type'" == "diagnostics" {
        _cate_dnn_plot_diagnostics, `common' `candidate_opts'
    }
    else if "`type'" == "tuning" {
        _cate_dnn_plot_tuning, `common' level(`level') ///
            `target_opts' `candidate_opts'
    }
    else if "`type'" == "weights" {
        _cate_dnn_plot_weights, `common' `view_opt' `target_opts' ///
            `candidate_opts' `logweights'
    }
    else if "`type'" == "outcome" {
        if "`candidates'`allcandidates'" != "" {
            di as error "Candidate selection is not available for type(outcome)."
            exit 198
        }
        _np_dnn_plot_outcome, `common' xindex(`xindex') `target_opts'
    }
    else {
        if "`candidates'`allcandidates'" != "" {
            di as error "Candidate selection is not available for type(dashboard)."
            exit 198
        }
        _np_dnn_plot_dashboard, `common' level(`level') ///
            `style_opt' xindex(`xindex')
    }

    local built_view "`r(view)'"
    local built_n = r(n_plotted)

    if `"`saving'"' != "" {
        local saving_lower = lower(`"`saving'"')
        local export_type ""
        foreach ext in gph png pdf svg tif jpg eps {
            if regexm(`"`saving_lower'"', "[.]`ext'") local export_type "`ext'"
        }
        if "`export_type'" == "gph" {
            if "`replace'" != "" graph save `name' `saving', replace
            else graph save `name' `saving'
        }
        else {
            local asopt ""
            if inlist("`export_type'", "png", "pdf", "svg", "tif", "jpg", "eps") ///
                local asopt "as(`export_type')"
            quietly graph display `name'
            if "`replace'" != "" graph export `saving', `asopt' replace
            else graph export `saving', `asopt'
        }
    }

    return local type "`type'"
    return local view "`built_view'"
    return local graph "`name'"
    return scalar n_plotted = `built_n'
    return local saved `saving'
end

program define _np_dnn_plot_outcome, rclass
    version 17
    syntax , NAME(name) XINDEX(integer) ///
        [TARGETS(numlist integer >=1) REPLACE NODRAW]

    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"
    tempname targets_mat diag plot
    matrix `targets_mat' = e(targets)
    matrix `diag' = e(diagnostics)
    local p = e(p)
    mata: st_matrix("`plot'", (st_matrix("`targets_mat'"), ///
        (1::rows(st_matrix("`targets_mat'"))), ///
        st_matrix("`diag'")[,13..17]))
    local xnames ""
    forvalues j = 1/`p' {
        local xnames "`xnames' x0_`j'"
    }
    matrix colnames `plot' = `xnames' target_id weighted_mean_y ///
        weighted_sd_y weighted_y_q10 weighted_y_q50 weighted_y_q90

    preserve
    clear
    quietly svmat double `plot', names(col)
    if "`targets'" != "" {
        quietly gen byte __keep = 0
        foreach j of local targets {
            quietly replace __keep = 1 if target_id == `j'
        }
        quietly keep if __keep
        quietly drop __keep
    }
    quietly count
    local nplot = r(N)

    local axis "target_id"
    local axis_label "Target"
    local varying ""
    if `xindex' > 0 local varying "x0_`xindex'"
    else {
        forvalues j = 1/`p' {
            quietly levelsof x0_`j', local(vals)
            local nv : word count `vals'
            if `nv' > 1 local varying "`varying' x0_`j'"
        }
        local nvary : word count `varying'
        if `nvary' != 1 local varying ""
    }
    if "`varying'" != "" {
        local axis "`varying'"
        local j = subinstr("`varying'", "x0_", "", .)
        local axis_label : word `j' of `e(xvars)'
    }

    if "`varying'" != "" & `nplot' >= 4 {
        twoway (rarea weighted_y_q10 weighted_y_q90 `axis', sort ///
                color(gs12%45) lcolor(none)) ///
            (line weighted_mean_y `axis', sort lcolor(navy) lwidth(medthick)), ///
            xtitle("`axis_label'") ytitle("Localized outcome") ///
            title("DNN localized outcome distribution") ///
            legend(order(2 "weighted mean" 1 "weighted 10th-90th percentiles")) ///
            `nameopt' `drawopt'
        local view "profile"
    }
    else {
        twoway (rcap weighted_y_q10 weighted_y_q90 `axis', lcolor(gs8)) ///
            (scatter weighted_mean_y `axis', mcolor(navy) msymbol(O)), ///
            xtitle("`axis_label'") ytitle("Localized outcome") ///
            title("DNN localized outcome distribution") ///
            legend(order(2 "weighted mean" 1 "weighted 10th-90th percentiles")) ///
            `nameopt' `drawopt'
        local view "pointrange"
    }
    restore
    return scalar n_plotted = `nplot'
    return local view "`view'"
end

program define _np_dnn_plot_dashboard, rclass
    version 17
    syntax , NAME(name) Level(cilevel) XINDEX(integer) ///
        [STYLE(string) REPLACE NODRAW]

    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"
    local style_opt ""
    if "`style'" != "" local style_opt "style(`style')"
    tempname gest gdiag gout gweights gtune

    quietly _cate_dnn_plot_estimates, name(`gest') replace nodraw ///
        level(`level') xindex(`xindex') `style_opt'
    quietly _cate_dnn_plot_diagnostics, name(`gdiag') replace nodraw
    quietly _np_dnn_plot_outcome, name(`gout') replace nodraw xindex(`xindex')
    quietly _cate_dnn_plot_weights, name(`gweights') replace nodraw view(cumulative)
    local graphs "`gest' `gdiag' `gout' `gweights'"
    local npanels 4
    if "`e(tuning_mode)'" == "grid" {
        quietly _cate_dnn_plot_tuning, name(`gtune') replace nodraw level(`level')
        local graphs "`graphs' `gtune'"
        local npanels 5
    }
    graph combine `graphs', cols(2) title("DNN nonparametric-regression dashboard") ///
        `nameopt' `drawopt'
    return scalar n_plotted = `npanels'
    return local view "dashboard"
end
