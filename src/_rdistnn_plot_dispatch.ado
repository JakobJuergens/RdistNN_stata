*! version 0.1.5  2026-10-25  Jakob Juergens
// One public option, sample-guard, graph export and return policy for scalar fits.
program define _rdistnn_plot_dispatch, rclass
    version 17
    syntax name(name=family) [, TYPE(string) VIEW(string) STYLE(string) ///
        XINDEX(string) TARGETS(numlist integer >=1) ///
        CANDIDATES(numlist integer >=1) ALLCandidates Level(string) ///
        SHARELimits(numlist min=2 max=2) NLABEL(string) LOGWeights ///
        NAME(name) REPLACE NODRAW SAVING(string asis)]

    if "`e(cmd)'" != "`family'_dnn" {
        _rdistnn_messages error, rc(301) message(`"`family'_dnn_plot requires active `family'_dnn estimates."')
    }
    _rdistnn_model_check
    if "`type'" == "" local type "estimates"
    local type = lower("`type'")
    local view = lower("`view'")
    local style = lower("`style'")
    if "`style'"!="" & !inlist("`style'", "auto", "profile", "pointrange", "forest") {
        _rdistnn_messages error, rc(198) message(`"style() must be auto, profile, pointrange, or forest."')
    }
    if !inlist("`type'", "estimates", "diagnostics", "tuning", "weights", "dashboard") & ///
        !("`family'"=="np" & "`type'"=="outcome") & ///
        !("`family'"!="np" & inlist("`type'", "overlap", "influence")) {
        _rdistnn_messages error, rc(198) message(`"Unsupported plot type() for `family'_dnn_plot."')
    }
    if "`candidates'"!="" & "`allcandidates'"!="" {
        _rdistnn_messages error, rc(198) message(`"Specify either candidates() or allcandidates, not both."')
    }
    local candidate_mode = ("`candidates'`allcandidates'"!="")
    if `candidate_mode' & !inlist("`type'", "estimates", "diagnostics", "tuning", "weights") {
        _rdistnn_messages error, rc(198) message(`"Candidate selection is not available for type(`type')."')
    }
    if "`targets'"!="" & "`type'"=="dashboard" {
        _rdistnn_messages error, rc(198) message(`"targets() is not available for type(dashboard); select a single plot type."')
    }
    if "`targets'"!="" & `candidate_mode' & "`type'"=="diagnostics" {
        _rdistnn_messages error, rc(198) message(`"Specify targets() or candidate selection for diagnostics, not both."')
    }
    if "`style'"!="" & !inlist("`type'", "estimates", "dashboard") {
        _rdistnn_messages error, rc(198) message(`"style() applies only to estimates and dashboard estimate panels."')
    }
    if "`view'"!="" & !inlist("`type'", "weights", "overlap") {
        _rdistnn_messages error, rc(198) message(`"view() applies only to weights and overlap plots."')
    }
    if "`xindex'"!="" & (!inlist("`type'", "estimates", "overlap", "outcome", "dashboard") | ///
        ("`type'"=="overlap" & "`view'"=="armess")) {
        _rdistnn_messages error, rc(198) message(`"xindex() does not apply to this plot."')
    }
    if "`xindex'"=="" local xindex 0
    confirm integer number `xindex'
    if `xindex'<0 | `xindex'>e(p) {
        _rdistnn_messages error, rc(198) message(`"xindex() must be 0 or a target-coordinate index in 1..`=e(p)'."')
    }
    if "`level'"!="" & !inlist("`type'", "estimates", "tuning", "dashboard") {
        _rdistnn_messages error, rc(198) message(`"level() applies only to estimates, tuning, or dashboard."')
    }
    if "`level'"=="" local level = c(level)
    quietly _rdistnn_plot_level, level(`level')
    if "`nlabel'"!="" & "`type'"!="influence" {
        _rdistnn_messages error, rc(198) message(`"nlabel() applies only to influence plots."')
    }
    if "`nlabel'"=="" local nlabel 3
    confirm integer number `nlabel'
    if `nlabel'<0 error 198
    if "`sharelimits'"!="" & !("`type'"=="overlap" & inlist("`view'", "", "mass")) {
        _rdistnn_messages error, rc(198) message(`"sharelimits() applies only to overlap view(mass)."')
    }
    if "`type'"=="weights" {
        if "`view'"=="" local view "profile"

        if `candidate_mode' & "`view'"!="profile" {
            _rdistnn_messages error, rc(198) message(`"Candidate weights require view(profile)."')
        }
    }
    if "`logweights'"!="" & !("`type'"=="weights" & "`view'"=="profile") {
        _rdistnn_messages error, rc(198) message(`"logweights applies only to weight view(profile)."')
    }
    foreach t of local targets {
        if `t'>e(m) {
            _rdistnn_messages error, rc(198) message(`"targets() entries must be in 1..`=e(m)'."')
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
    local helper_rc=_rc
    if !inlist(`helper_rc',0,111) exit `helper_rc'
    if `helper_rc'==111 {
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
    else if "`type'"=="weights" _cate_dnn_plot_weights, `common' `topt' `copt' `vopt' `logweights' validated
    else if "`type'"=="influence" _cate_dnn_plot_influence, `common' `topt' nlabel(`nlabel') validated
    else if "`type'"=="outcome" _np_dnn_plot_outcome, `common' xindex(`xindex') `topt'
    else if "`family'"=="np" _np_dnn_plot_dashboard, `common' level(`level') xindex(`xindex') `sopt' validated
    else _cate_dnn_plot_dashboard, `common' level(`level') xindex(`xindex') `sopt' validated

    tempname plot_data panel_info
    local built_view "`r(view)'"
    local built_panels = r(n_panels)
    local panels "`r(panels)'"
    local built_rows = r(n_rows)
    local omitted = r(n_omitted)
    local selected_targets "`r(selected_targets)'"
    if "`type'"!="dashboard" & `built_rows'!=0 matrix `plot_data' = r(plot_data)
    else if "`type'"=="dashboard" matrix `panel_info' = r(panel_info)
    if `"`saving'"'!="" quietly _rdistnn_plot_export, name(`name') saving(`saving') `replace'
    _rdistnn_report_notices, `topt'
    return clear
    if "`type'"!="dashboard" {
        if `built_rows'!=0 {
        return scalar n_rows = rowsof(`plot_data')
        local columns : colnames `plot_data'
        local tid : list posof "target_id" in columns
        local cid : list posof "candidate_id" in columns
        if `tid' {
            mata: rdistnn_plot_target_ids("`plot_data'",`tid',"plotted_targets")
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
        return scalar n_panels = `built_panels'
        return local panels "`panels'"
    }
    if inlist("`type'", "estimates", "tuning") {
        return scalar n_omitted = `omitted'
        return local selected_targets "`selected_targets'"
    }
    return local type "`type'"
    return local view "`built_view'"
    return local graph "`name'"
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
        _rdistnn_messages error, rc(498) message(`"Candidate plots require tuningmode(grid)."')
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
            _rdistnn_messages error, rc(198) message(`"candidates() entries must identify rows of e(tuning_candidates)."')
        }
        if el(`meta',`g',6)!=1 {
            _rdistnn_messages error, rc(498) message(`"Candidate `g' is unavailable; select status=1 candidates."')
        }
    }
    if "`selected'"=="" _rdistnn_messages error, rc(498) message("No usable tuning candidate is available; inspect the retained candidate status.")
    return local candidates "`selected'"
end

program define _rdistnn_plot_sample_guard
    version 17
    if e(model_version)!=4 {
        _rdistnn_messages error, rc(498) message(`"Observation-level plots require current fitted model state; refit `e(cmd)' with the current package."')
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
            _rdistnn_messages error, rc(198) message(`"saving() requires .gph, .png, .pdf, .svg, .tif, .jpg, or .eps."')
        }
        graph display `name'
        graph export `saving', `asopt' `replace'
    }
end

// Record each dashboard panel without claiming a common data table.
program define _rdistnn_plot_add_panel
    version 17
    syntax , MATRIX(name) ID(integer)
    local nr = r(n_rows)
    if `id'==1 matrix `matrix' = (`id', `nr')
    else matrix `matrix' = `matrix' \ (`id', `nr')
    matrix colnames `matrix' = panel_id n_rows
end
