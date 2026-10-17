*! version 0.1.0  2026-10-17  Jakob Juergens
*! Stata-native postestimation plots for clate_dnn

program define clate_dnn_plot, rclass
    version 17

    syntax [, TYPE(string) VIEW(string) STYLE(string) XINDEX(integer 0) ///
        TARGETS(numlist integer >=1) CANDIDATES(numlist integer >=1) ///
        ALLCandidates Level(cilevel) SHARELimits(numlist min=2 max=2) ///
        NLABEL(integer 3) LOGWeights NAME(name) REPLACE NODRAW ///
        SAVING(string asis)]

    if "`e(cmd)'" != "clate_dnn" {
        di as error "clate_dnn_plot requires active clate_dnn estimates."
        exit 301
    }
    if "`type'" == "" local type "estimates"
    local type = lower("`type'")
    if !inlist("`type'", "estimates", "diagnostics", "tuning", ///
        "overlap", "weights", "influence", "dashboard") {
        di as error "type() must be estimates, diagnostics, tuning, overlap, weights, influence, or dashboard."
        exit 198
    }
    if "`candidates'" != "" & "`allcandidates'" != "" {
        di as error "Specify either candidates() or allcandidates, not both."
        exit 198
    }
    capture program list _cate_dnn_plot_estimates
    if _rc {
        findfile "_rdistnn_plot_helpers.ado"
        quietly run "`r(fn)'"
    }
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"
    if `xindex' < 0 | `xindex' > e(p) {
        di as error "xindex() must be 0 or a target-coordinate index in 1..`=e(p)'."
        exit 198
    }
    if `nlabel' < 0 {
        di as error "nlabel() must be a nonnegative integer."
        exit 198
    }
    foreach j of local targets {
        if `j' > e(m) {
            di as error "targets() entries must be in 1..`=e(m)'."
            exit 198
        }
    }

    if "`type'" == "influence" | ("`type'" == "weights" & lower("`view'")=="kernel") {
        capture quietly checkestimationsample
        if _rc | "`e(datasignature)'"=="" {
            di as error "The estimation data have changed or cannot be verified; refit clate_dnn before plotting observation-level results."
            exit 498
        }
        capture quietly _rdistnn_clate_signature `e(signature_vars)' if e(sample)
        if _rc | "`r(signature)'"!="`e(ordered_signature)'" | "`e(ordered_signature)'"=="" {
            di as error "The estimation data or their order have changed; refit clate_dnn before plotting observation-level results."
            exit 498
        }
    }
    if "`name'" == "" local name "clate_dnn_`type'"
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
    local share_opt ""
    if "`sharelimits'" != "" local share_opt "sharelimits(`sharelimits')"

    if "`type'" == "estimates" {
        _cate_dnn_plot_estimates, `common' level(`level') ///
            `style_opt' xindex(`xindex') `target_opts' `candidate_opts'
    }
    else if "`type'" == "diagnostics" {
        _cate_dnn_plot_diagnostics, `common' `target_opts' `candidate_opts'
    }
    else if "`type'" == "tuning" {
        _cate_dnn_plot_tuning, `common' level(`level') ///
            `target_opts' `candidate_opts'
    }
    else if "`type'" == "overlap" {
        _cate_dnn_plot_overlap, `common' `view_opt' xindex(`xindex') ///
            `share_opt' `target_opts'
    }
    else if "`type'" == "weights" {
        _cate_dnn_plot_weights, `common' `view_opt' `target_opts' ///
            `candidate_opts' `logweights'
    }
    else if "`type'" == "influence" {
        _cate_dnn_plot_influence, `common' `target_opts' nlabel(`nlabel')
    }
    else {
        _cate_dnn_plot_dashboard, `common' level(`level') ///
            `style_opt' xindex(`xindex')
    }

    tempname plot_data
    capture matrix `plot_data' = r(plot_data)
    local has_plot_data = (_rc==0)
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

    if `has_plot_data' return matrix plot_data = `plot_data'
    return local type "`type'"
    return local view "`built_view'"
    return local graph "`name'"
    return scalar n_plotted = `built_n'
    return local saved `saving'
end
