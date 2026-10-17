*! version 0.1.0  2026-10-17  Jakob Juergens
// One public option, sample-guard, graph export and return policy for scalar fits.
program define _rdistnn_plot_dispatch, rclass
    version 17
    syntax name(name=family) [, TYPE(string) VIEW(string) STYLE(string) ///
        XINDEX(string) TARGETS(numlist integer >=1) ///
        CANDIDATES(numlist integer >=1) ALLCandidates Level(string) ///
        SHARELimits(numlist min=2 max=2) NLABEL(string) LOGWeights ///
        NAME(name) REPLACE NODRAW SAVING(string asis)]

    if "`e(cmd)'" != "`family'_dnn" {
        di as error "`family'_dnn_plot requires active `family'_dnn estimates."
        exit 301
    }
    if "`type'" == "" local type "estimates"
    local type = lower("`type'")
    local view = lower("`view'")
    local style = lower("`style'")
    if "`style'"!="" & !inlist("`style'", "auto", "profile", "pointrange", "forest") {
        di as error "style() must be auto, profile, pointrange, or forest."
        exit 198
    }
    if !inlist("`type'", "estimates", "diagnostics", "tuning", "weights", "dashboard") & ///
        !("`family'"=="np" & "`type'"=="outcome") & ///
        !("`family'"!="np" & inlist("`type'", "overlap", "influence")) {
        di as error "Unsupported plot type() for `family'_dnn_plot."
        exit 198
    }
    if "`candidates'"!="" & "`allcandidates'"!="" {
        di as error "Specify either candidates() or allcandidates, not both."
        exit 198
    }
    local candidate_mode = ("`candidates'`allcandidates'"!="")
    if `candidate_mode' & !inlist("`type'", "estimates", "diagnostics", "tuning", "weights") {
        di as error "Candidate selection is not available for type(`type')."
        exit 198
    }
    if "`targets'"!="" & "`type'"=="dashboard" {
        di as error "targets() is not available for type(dashboard); select a single plot type."
        exit 198
    }
    if "`targets'"!="" & `candidate_mode' & "`type'"=="diagnostics" {
        di as error "Specify targets() or candidate selection for diagnostics, not both."
        exit 198
    }
    if "`style'"!="" & !inlist("`type'", "estimates", "dashboard") {
        di as error "style() applies only to estimates and dashboard estimate panels."
        exit 198
    }
    if "`view'"!="" & !inlist("`type'", "weights", "overlap") {
        di as error "view() applies only to weights and overlap plots."
        exit 198
    }
    if "`xindex'"!="" & (!inlist("`type'", "estimates", "overlap", "outcome", "dashboard") | ///
        ("`type'"=="overlap" & "`view'"=="armess")) {
        di as error "xindex() does not apply to this plot."
        exit 198
    }
    if "`xindex'"=="" local xindex 0
    capture confirm integer number `xindex'
    if _rc error 198
    if `xindex'<0 | `xindex'>e(p) {
        di as error "xindex() must be 0 or a target-coordinate index in 1..`=e(p)'."
        exit 198
    }
    if "`level'"!="" & !inlist("`type'", "estimates", "tuning", "dashboard") {
        di as error "level() applies only to estimates, tuning, or dashboard."
        exit 198
    }
    if "`level'"=="" local level = c(level)
    quietly _rdistnn_plot_level, level(`level')
    if "`nlabel'"!="" & "`type'"!="influence" {
        di as error "nlabel() applies only to influence plots."
        exit 198
    }
    if "`nlabel'"=="" local nlabel 3
    capture confirm integer number `nlabel'
    if _rc error 198
    if `nlabel'<0 error 198
    if "`sharelimits'"!="" & !("`type'"=="overlap" & inlist("`view'", "", "mass")) {
        di as error "sharelimits() applies only to overlap view(mass)."
        exit 198
    }
    if "`type'"=="weights" {
        if "`view'"=="" local view "profile"
        if "`targets'"!="" & "`view'"!="kernel" {
            di as error "targets() applies to weight view(kernel); rank weights are common to every target."
            exit 198
        }
        if `candidate_mode' & "`view'"!="profile" {
            di as error "Candidate weights require view(profile)."
            exit 198
        }
    }
    if "`logweights'"!="" & !("`type'"=="weights" & "`view'"=="profile") {
        di as error "logweights applies only to weight view(profile)."
        exit 198
    }
    foreach t of local targets {
        if `t'>e(m) {
            di as error "targets() entries must be in 1..`=e(m)'."
            exit 198
        }
    }
    local targets : list uniq targets
    local candidates : list uniq candidates
    if `candidate_mode' | "`type'"=="tuning" {
        quietly _rdistnn_plot_candidates, `=cond("`candidates'"!="","candidates(`candidates')", "")' `allcandidates'
        local selected_candidates "`r(candidates)'"
    }
    if "`type'"=="influence" | ("`type'"=="weights" & "`view'"=="kernel") {
        quietly _rdistnn_plot_sample_guard
    }

    capture program list _cate_dnn_plot_estimates
    if _rc {
        findfile "_rdistnn_plot_helpers.ado"
        quietly run "`r(fn)'"
    }
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"
    if "`name'"=="" local name "`family'_dnn_`type'"
    local common "name(`name') `replace' `nodraw'"
    local topt ""
    if "`targets'"!="" local topt "targets(`targets')"
    local copt ""
    if `candidate_mode' local copt "candidates(`selected_candidates')"
    local sopt ""
    if "`style'"!="" local sopt "style(`style')"
    local vopt ""
    if "`view'"!="" local vopt "view(`view')"
    local shareopt ""
    if "`sharelimits'"!="" local shareopt "sharelimits(`sharelimits')"
    return clear
    if "`type'"=="estimates" _cate_dnn_plot_estimates, `common' level(`level') xindex(`xindex') `topt' `copt' `sopt'
    else if "`type'"=="diagnostics" _cate_dnn_plot_diagnostics, `common' `topt' `copt'
    else if "`type'"=="tuning" _cate_dnn_plot_tuning, `common' level(`level') `topt' candidates(`selected_candidates')
    else if "`type'"=="overlap" _cate_dnn_plot_overlap, `common' xindex(`xindex') `topt' `vopt' `shareopt'
    else if "`type'"=="weights" _cate_dnn_plot_weights, `common' `topt' `copt' `vopt' `logweights'
    else if "`type'"=="influence" _cate_dnn_plot_influence, `common' `topt' nlabel(`nlabel')
    else if "`type'"=="outcome" _np_dnn_plot_outcome, `common' xindex(`xindex') `topt'
    else if "`family'"=="np" _np_dnn_plot_dashboard, `common' level(`level') xindex(`xindex') `sopt'
    else _cate_dnn_plot_dashboard, `common' level(`level') xindex(`xindex') `sopt'

    tempname plot_data panel_info
    local built_view "`r(view)'"
    local built_n = r(n_plotted)
    local panels "`r(panels)'"
    local built_rows = r(n_rows)
    local omitted = r(n_omitted)
    local selected_targets "`r(selected_targets)'"
    if "`type'"!="dashboard" & `built_rows'!=0 matrix `plot_data' = r(plot_data)
    else if "`type'"=="dashboard" matrix `panel_info' = r(panel_info)
    if `"`saving'"'!="" quietly _rdistnn_plot_export, name(`name') saving(`saving') `replace'
    return clear
    if "`type'"!="dashboard" {
        if `built_rows'!=0 {
        return scalar n_rows = rowsof(`plot_data')
        local columns : colnames `plot_data'
        local tid : list posof "target_id" in columns
        local cid : list posof "candidate_id" in columns
        if `tid' {
            mata: _rdnn_ids=st_matrix("`plot_data'")[,`tid']; st_local("plotted_targets", invtokens(strofreal(uniqrows(select(_rdnn_ids,_rdnn_ids:<.)))'))
            if "`plotted_targets'"!="" return local targets "`plotted_targets'"
        }
        if `cid' {
            mata: st_local("plotted_candidates", invtokens(strofreal(uniqrows(st_matrix("`plot_data'")[,`cid']))'))
            return local candidates "`plotted_candidates'"
        }
        return matrix plot_data = `plot_data'
        }
        else return scalar n_rows = 0
    }
    else {
        return matrix panel_info = `panel_info'
        return scalar n_rows = .
        return scalar n_panels = `built_n'
        return local panels "`panels'"
    }
    if inlist("`type'", "estimates", "tuning") {
        return scalar n_omitted = `omitted'
        return local selected_targets "`selected_targets'"
    }
    return local type "`type'"
    return local view "`built_view'"
    return local graph "`name'"
    return scalar n_plotted = `built_n'
    return local saved `saving'
end

program define _rdistnn_plot_level
    version 17
    syntax , Level(cilevel)
end

program define _rdistnn_plot_candidates, rclass
    version 17
    syntax [, CANDIDATES(numlist integer >=1) ALLCandidates]
    if "`e(tuning_mode)'"!="grid" {
        di as error "Candidate plots require tuningmode(grid)."
        exit 498
    }
    tempname meta
    matrix `meta' = e(tuning_candidates)
    local selected "`candidates'"
    if "`selected'"=="" {
        forvalues g=1/`=rowsof(`meta')' {
            if el(`meta',`g',6)==1 local selected "`selected' `g'"
        }
    }
    foreach g of local selected {
        if `g'>rowsof(`meta') {
            di as error "candidates() entries must identify rows of e(tuning_candidates)."
            exit 198
        }
        if el(`meta',`g',6)!=1 {
            di as error "Candidate `g' is unavailable; select status=1 candidates."
            exit 498
        }
    }
    if "`selected'"=="" error 498
    return local candidates "`selected'"
end

program define _rdistnn_plot_sample_guard
    version 17
    capture quietly checkestimationsample
    if _rc | "`e(datasignature)'"=="" {
        di as error "The estimation data have changed or cannot be verified; refit `e(cmd)' before plotting observation-level results."
        exit 498
    }
    if "`e(ordered_signature)'"=="" | "`e(signature_vars)'"=="" {
        di as error "This saved fit cannot verify observation order; refit `e(cmd)' before plotting observation-level results."
        exit 498
    }
    capture quietly _rdistnn_sample_signature `e(signature_vars)' if e(sample)
    if _rc | "`r(signature)'"!="`e(ordered_signature)'" {
        di as error "The estimation data or their order have changed; refit `e(cmd)' before plotting observation-level results."
        exit 498
    }
end

program define _rdistnn_plot_export
    version 17
    syntax , NAME(name) SAVING(string asis) [REPLACE]
    local path `saving'
    local path = lower(`"`path'"')
    if regexm(`"`path'"', "[.]gph$") graph save `name' `saving', `replace'
    else {
        local asopt ""
        foreach ext in png pdf svg tif jpg eps {
            if regexm(`"`path'"', "[.]`ext'$") local asopt "as(`ext')"
        }
        if "`asopt'"=="" {
            di as error "saving() requires .gph, .png, .pdf, .svg, .tif, .jpg, or .eps."
            exit 198
        }
        graph display `name'
        graph export `saving', `asopt' `replace'
    }
end

// Record each dashboard panel without claiming a common data table.
program define _rdistnn_plot_add_panel
    version 17
    syntax , MATRIX(name) ID(integer)
    tempname data
    local n = r(n_plotted)
    local nr = r(n_rows)
    if `nr'!=0 {
        matrix `data' = r(plot_data)
        local nr = rowsof(`data')
    }
    if `id'==1 matrix `matrix' = (`id', `nr', `n')
    else matrix `matrix' = `matrix' \ (`id', `nr', `n')
    matrix colnames `matrix' = panel_id n_rows n_plotted
end
