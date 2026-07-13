*! version 0.1.1  2026-07-12  Jakob Juergens
*! Stata-native postestimation plots for cate_dnn

program define cate_dnn_plot, rclass
    version 17

    syntax [, TYPE(string) VIEW(string) STYLE(string) XINDEX(integer 0) ///
        TARGETS(numlist integer >=1) CANDIDATES(numlist integer >=1) ///
        ALLCandidates Level(cilevel) SHARELimits(numlist min=2 max=2) ///
        NLABEL(integer 3) LOGWeights NAME(name) REPLACE NODRAW ///
        SAVING(string asis)]

    if "`e(cmd)'" != "cate_dnn" {
        di as error "cate_dnn_plot requires active cate_dnn estimates."
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

    if "`name'" == "" local name "cate_dnn_`type'"
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
        _cate_dnn_plot_diagnostics, `common' `candidate_opts'
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

program define _cate_dnn_plot_estimates, rclass
    version 17
    syntax , NAME(name) Level(cilevel) XINDEX(integer) ///
        [STYLE(string) TARGETS(numlist integer >=1) ///
        CANDIDATES(numlist integer >=1) ALLCandidates REPLACE NODRAW]

    if "`style'" == "" local style "auto"
    local style = lower("`style'")
    if !inlist("`style'", "auto", "profile", "pointrange", "forest") {
        di as error "style() must be auto, profile, pointrange, or forest."
        exit 198
    }
    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"

    tempname targets_mat bmat vmat semat plotmat candidate_meta grid_est grid_se
    matrix `targets_mat' = e(targets)
    local m = e(m)
    local p = e(p)
    local xvars "`e(xvars)'"
    local z = invnormal(1 - (1 - `level'/100) / 2)
    local use_candidates = ("`candidates'`allcandidates'" != "")
    local estimand "CATE"
    local plot_title "DNN CATE estimates"
    if "`e(cmd)'" == "np_dnn" {
        local estimand "Conditional mean"
        local plot_title "DNN nonparametric regression estimates"
    }

    if !`use_candidates' {
        matrix `bmat' = e(b)'
        capture matrix `vmat' = e(V)
        if _rc matrix `semat' = J(`m', 1, .)
        else if rowsof(`vmat') != `m' | colsof(`vmat') != `m' {
            matrix `semat' = J(`m', 1, .)
        }
        else mata: st_matrix("`semat'", sqrt(diagonal(st_matrix("`vmat'"))))
        mata: st_matrix("`plotmat'", rdistnn_plot_est_default(st_matrix("`targets_mat'"), st_matrix("`bmat'"), st_matrix("`semat'"), `=e(s)', `z'))
        local selected "0"
    }
    else {
        if "`e(tuning_mode)'" != "grid" {
            di as error "Candidate estimate plots require the active estimator to use tuningmode(grid)."
            exit 498
        }
        matrix `candidate_meta' = e(tuning_candidates)
        matrix `grid_est' = e(tuning_estimates)
        matrix `grid_se' = e(tuning_se)
        local gmax = rowsof(`candidate_meta')
        local selected "`candidates'"
        if "`allcandidates'" != "" {
            local selected ""
            forvalues g = 1/`gmax' {
                if el(`candidate_meta', `g', 6) == 1 local selected "`selected' `g'"
            }
        }
        foreach g of local selected {
            if `g' > `gmax' {
                di as error "candidates() entries must be in 1..`gmax'."
                exit 198
            }
            if el(`candidate_meta', `g', 6) != 1 {
                di as error "Candidate `g' is unavailable; select a candidate with status=1."
                exit 498
            }
        }
        mata: st_matrix("`plotmat'", rdistnn_plot_est_grid(st_matrix("`targets_mat'"), st_matrix("`candidate_meta'"), st_matrix("`grid_est'"), st_matrix("`grid_se'"), strtoreal(tokens("`selected'"))', `z'))
    }
    local ng : word count `selected'

    local xnames ""
    forvalues j = 1/`p' {
        local xnames "`xnames' x0_`j'"
    }
    matrix colnames `plotmat' = `xnames' candidate_id s is_default ///
        target_id estimate se ci_lower ci_upper

    preserve
    clear
    quietly svmat double `plotmat', names(col)
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
    local mt = 0
    quietly levelsof target_id, local(target_levels)
    local mt : word count `target_levels'

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
        local axis_label : word `j' of `xvars'
    }
    if "`style'" == "auto" {
        if "`varying'" != "" & `mt' >= 4 local style "profile"
        else if "`varying'" == "" & `mt' > 1 & `mt' <= 12 local style "forest"
        else local style "pointrange"
    }
    if "`style'" == "profile" & "`varying'" == "" {
        restore
        di as error "style(profile) requires targets varying in exactly one coordinate or xindex()."
        exit 198
    }

    local plots ""
    local legend_order ""
    local layercount 0
    local plotnum 0
    foreach g of local selected {
        local ++plotnum
        local ps = mod(`plotnum' - 1, 15) + 1
        local width "medthin"
        local symbol "Oh"
        quietly summarize is_default if candidate_id == `g', meanonly
        if r(max) == 1 {
            local width "medthick"
            local symbol "O"
        }
        if "`style'" == "profile" {
            capture assert missing(se) if candidate_id == `g'
            if _rc {
                local plots `"`plots' (rarea ci_lower ci_upper `axis' if candidate_id==`g', sort color(gs12%45) lcolor(none))"'
                local ++layercount
            }
            local plots `"`plots' (line estimate `axis' if candidate_id==`g', sort pstyle(p`ps') lwidth(`width'))"'
            local ++layercount
        }
        else if "`style'" == "forest" {
            quietly gen double __y`plotnum' = target_id + ///
                (`plotnum' - (`ng' + 1)/2) * .12
            capture assert missing(se) if candidate_id == `g'
            if _rc {
                local plots `"`plots' (rcap ci_lower ci_upper __y`plotnum' if candidate_id==`g', horizontal pstyle(p`ps'))"'
                local ++layercount
            }
            local plots `"`plots' (scatter __y`plotnum' estimate if candidate_id==`g', msymbol(`symbol') pstyle(p`ps'))"'
            local ++layercount
        }
        else {
            capture assert missing(se) if candidate_id == `g'
            if _rc {
                local plots `"`plots' (rcap ci_lower ci_upper `axis' if candidate_id==`g', pstyle(p`ps'))"'
                local ++layercount
            }
            local plots `"`plots' (scatter estimate `axis' if candidate_id==`g', msymbol(`symbol') pstyle(p`ps'))"'
            local ++layercount
        }
        if `use_candidates' {
            local sg = el(`candidate_meta', `g', 8)
            local legend_order `"`legend_order' `layercount' "s = `sg'""'
        }
    }
    local legendopt "legend(off)"
    if `ng' > 1 local legendopt `"legend(order(`legend_order') cols(1))"'

    local subtitle "n = `=e(N)', k = `=e(k)', s = `=e(s)'"
    if "`style'" == "forest" {
        twoway `plots', xline(0, lpattern(dash) lcolor(gs8)) ///
            xtitle("`estimand' estimate") ytitle("Target") ///
            title("`plot_title'") subtitle("`subtitle'") ///
            `legendopt' `nameopt' `drawopt'
    }
    else {
        twoway `plots', yline(0, lpattern(dash) lcolor(gs8)) ///
            xtitle("`axis_label'") ytitle("`estimand' estimate") ///
            title("`plot_title'") subtitle("`subtitle'") ///
            `legendopt' `nameopt' `drawopt'
    }
    restore
    return scalar n_plotted = `nplot'
    return local view "`style'"
end

program define _cate_dnn_plot_tuning, rclass
    version 17
    syntax , NAME(name) Level(cilevel) [TARGETS(numlist integer >=1) ///
        CANDIDATES(numlist integer >=1) ALLCandidates REPLACE NODRAW]

    if "`e(tuning_mode)'" != "grid" {
        di as error "type(tuning) requires the active estimator to use tuningmode(grid)."
        exit 498
    }
    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"
    tempname meta est ses plot
    matrix `meta' = e(tuning_candidates)
    matrix `est' = e(tuning_estimates)
    matrix `ses' = e(tuning_se)
    local gmax = rowsof(`meta')
    local m = e(m)
    local z = invnormal(1 - (1 - `level'/100) / 2)
    local estimand "CATE"
    local plot_title "CATE tuning sensitivity"
    if "`e(cmd)'" == "np_dnn" {
        local estimand "Conditional mean"
        local plot_title "Nonparametric-regression tuning sensitivity"
    }
    local selected "`candidates'"
    if "`allcandidates'" != "" | "`selected'" == "" {
        local selected ""
        forvalues g = 1/`gmax' {
            if el(`meta', `g', 6) == 1 local selected "`selected' `g'"
        }
    }
    foreach g of local selected {
        if `g' > `gmax' | el(`meta', `g', 6) != 1 {
            di as error "Selected tuning candidates must be available rows in 1..`gmax'."
            exit 498
        }
    }
    local ng : word count `selected'
    local target_list "`targets'"
    if "`target_list'" == "" {
        if `m' <= 3 {
            forvalues j = 1/`m' {
                local target_list "`target_list' `j'"
            }
        }
        else local target_list "1 `=ceil(`m'/2)' `m'"
    }
    mata: st_matrix("`plot'", rdistnn_plot_tuning(st_matrix("`meta'"), st_matrix("`est'"), st_matrix("`ses'"), strtoreal(tokens("`selected'"))', strtoreal(tokens("`target_list'"))', `z'))
    matrix colnames `plot' = candidate_id s lambda is_default target_id ///
        estimate se ci_lower ci_upper
    preserve
    clear
    quietly svmat double `plot', names(col)
    capture assert missing(se)
    local hasse = (_rc != 0)
    if `hasse' {
        twoway (rcap ci_lower ci_upper s, lcolor(gs9)) ///
            (connected estimate s, sort mcolor(black) lcolor(black)) ///
            (scatter estimate s if is_default, msymbol(O) mcolor(navy)), ///
            by(target_id, note("") title("`plot_title'")) ///
            xline(`=e(s)', lpattern(dot) lcolor(gs8)) ///
            xtitle("Kernel order s") ytitle("`estimand' estimate") ///
            legend(order(2 "candidate path" 3 "fitted s")) ///
            `nameopt' `drawopt'
    }
    else {
        twoway (connected estimate s, sort mcolor(black) lcolor(black)) ///
            (scatter estimate s if is_default, msymbol(O) mcolor(navy)), ///
            by(target_id, note("") title("`plot_title'")) ///
            xline(`=e(s)', lpattern(dot) lcolor(gs8)) ///
            xtitle("Kernel order s") ytitle("`estimand' estimate") ///
            legend(order(1 "candidate path" 2 "fitted s")) ///
            `nameopt' `drawopt'
    }
    restore
    return scalar n_plotted = `ng'
    return local view "tuning"
end

program define _cate_dnn_plot_diagnostics, rclass
    version 17
    syntax , NAME(name) [CANDIDATES(numlist integer >=1) ///
        ALLCandidates REPLACE NODRAW]

    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"
    tempname plot raw meta diag g1 g2 g3
    local candidate_mode = ("`candidates'`allcandidates'" != "")

    if !`candidate_mode' {
        matrix `raw' = e(diagnostics)
        local m = rowsof(`raw')
        if "`e(cmd)'" == "cate_dnn" {
            mata: st_matrix("`plot'", ((1::rows(st_matrix("`raw'"))), st_matrix("`raw'")[,(1,2,4,20,25)], J(rows(st_matrix("`raw'")),1,0)))
            matrix colnames `plot' = plot_id ess max_weight top5_weight_share ///
                physical_radius_90 top5_score_share is_default
        }
        else {
            mata: st_matrix("`plot'", ((1::rows(st_matrix("`raw'"))), st_matrix("`raw'")[,(1,2,4)], J(rows(st_matrix("`raw'")),1,0)))
            matrix colnames `plot' = plot_id ess max_weight top5_weight_share is_default
        }
        local xvar "plot_id"
        local xlabel "Target"
        local nplot = `m'
        local guide = .9 * e(N)
    }
    else {
        if "`e(tuning_mode)'" != "grid" {
            di as error "Candidate diagnostic plots require the active estimator to use tuningmode(grid)."
            exit 498
        }
        matrix `meta' = e(tuning_candidates)
        matrix `diag' = e(tuning_diagnostics)
        local gmax = rowsof(`meta')
        local selected "`candidates'"
        if "`allcandidates'" != "" {
            local selected ""
            forvalues g = 1/`gmax' {
                if el(`meta', `g', 6) == 1 local selected "`selected' `g'"
            }
        }
        foreach g of local selected {
            if `g' > `gmax' | el(`meta', `g', 6) != 1 {
                di as error "Selected diagnostic candidates must be available rows in 1..`gmax'."
                exit 498
            }
        }
        mata: st_matrix("`plot'", (st_matrix("`meta'")[strtoreal(tokens("`selected'"))',8], st_matrix("`diag'")[strtoreal(tokens("`selected'"))',(1,2,4)], st_matrix("`meta'")[strtoreal(tokens("`selected'"))',3]))
        matrix colnames `plot' = plot_id ess max_weight top5_weight_share is_default
        local xvar "plot_id"
        local xlabel "Kernel order s"
        local nplot : word count `selected'
        local guide = .9 * e(N)
    }

    preserve
    clear
    quietly svmat double `plot', names(col)
    twoway (connected ess `xvar', sort mcolor(black) lcolor(black)) ///
        (scatter ess `xvar' if is_default, mcolor(navy) msymbol(O)), ///
        yscale(range(0 `guide')) ///
        ylabel(0 `guide', format(%9.0g)) ///
        yline(`guide', lpattern(dot) lcolor(gs8)) ///
        xtitle("`xlabel'") ytitle("Effective sample size") ///
        title("Effective sample size") legend(off) name(`g1', replace) nodraw
    twoway (connected max_weight `xvar', sort mcolor(black) lcolor(black)) ///
        (scatter max_weight `xvar' if is_default, mcolor(navy) msymbol(O)), ///
        yscale(range(0 1)) ylabel(0(.2)1) ///
        xtitle("`xlabel'") ytitle("Maximum weight") ///
        title("Maximum weight") legend(off) name(`g2', replace) nodraw
    twoway (connected top5_weight_share `xvar', sort mcolor(black) lcolor(black)) ///
        (scatter top5_weight_share `xvar' if is_default, mcolor(navy) msymbol(O)), ///
        yscale(range(0 1)) ylabel(0(.2)1) ///
        xtitle("`xlabel'") ytitle("Top-5 weight share") ///
        title("Weight concentration") legend(off) name(`g3', replace) nodraw
    local plot_title "DNN CATE diagnostics"
    if "`e(cmd)'" == "np_dnn" local plot_title "DNN nonparametric-regression diagnostics"
    graph combine `g1' `g2' `g3', cols(3) ///
        title("`plot_title'") `nameopt' `drawopt'
    if !`candidate_mode' & "`e(cmd)'" == "cate_dnn" {
        quietly summarize physical_radius_90, meanonly
        di as text "90% physical-radius range: " as result %9.4g r(min) ///
            as text " to " as result %9.4g r(max)
        quietly summarize top5_score_share, meanonly
        di as text "Top-5 absolute score-contribution share range: " ///
            as result %6.3f r(min) as text " to " as result %6.3f r(max)
    }
    restore
    return scalar n_plotted = `nplot'
    if `candidate_mode' return local view "candidates"
    else return local view "targets"
end

program define _cate_dnn_plot_overlap, rclass
    version 17
    syntax , NAME(name) XINDEX(integer) [VIEW(string) ///
        SHARELimits(numlist min=2 max=2) TARGETS(numlist integer >=1) ///
        REPLACE NODRAW]

    if "`view'" == "" local view "mass"
    local view = lower("`view'")
    if !inlist("`view'", "mass", "armess", "propensity") {
        di as error "For type(overlap), view() must be mass, armess, or propensity."
        exit 198
    }
    if "`sharelimits'" == "" local sharelimits ".1 .9"
    local lower : word 1 of `sharelimits'
    local upper : word 2 of `sharelimits'
    if `lower' < 0 | `upper' > 1 | `lower' >= `upper' {
        di as error "sharelimits() must contain two increasing values in [0,1]."
        exit 198
    }
    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"

    tempname diag targets_mat plot
    matrix `diag' = e(diagnostics)
    matrix `targets_mat' = e(targets)
    local m = e(m)
    local p = e(p)
    local xvars "`e(xvars)'"
    if "`view'" == "mass" {
        mata: st_matrix("`plot'", rdistnn_plot_overlap(st_matrix("`diag'"), st_matrix("`targets_mat'"), (12,13,14,15)))
        matrix colnames `plot' = `: colnames `targets_mat'' target_id ///
            overlap_available treated_mass control_mass treated_share
    }
    else if "`view'" == "armess" {
        mata: st_matrix("`plot'", rdistnn_plot_overlap(st_matrix("`diag'"), st_matrix("`targets_mat'"), (21,22,23)))
        matrix colnames `plot' = `: colnames `targets_mat'' target_id ///
            arm_ess_available treated_ess control_ess
    }
    else {
        mata: st_matrix("`plot'", rdistnn_plot_overlap(st_matrix("`diag'"), st_matrix("`targets_mat'"), (16,17,18,19)))
        matrix colnames `plot' = `: colnames `targets_mat'' target_id ///
            propensity_available propensity_mean propensity_q10 propensity_q90
    }
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
    if `xindex' > 0 {
        local axis : word `xindex' of `xvars'
        local axis_label "`axis'"
    }
    else {
        local j = 0
        foreach xv of local xvars {
            local ++j
            quietly levelsof `xv', local(vals)
            local nv : word count `vals'
            if `nv' > 1 local varying "`varying' `xv'"
        }
        local nvary : word count `varying'
        if `nvary' == 1 {
            local axis "`varying'"
            local axis_label "`axis'"
        }
    }

    if "`view'" == "mass" {
        quietly count if overlap_available == 1
        if r(N) == 0 {
            restore
            di as error "Overlap mass is unavailable for a precomputed pseudo-outcome fit."
            exit 498
        }
        if "`axis'" != "target_id" & `nplot' >= 4 {
            twoway (line treated_mass `axis', sort lcolor(navy)) ///
                (line control_mass `axis', sort lcolor(orange)), ///
                yscale(range(0 1)) ylabel(0(.2)1) ///
                yline(`lower' `upper', lpattern(dot) lcolor(gs8)) ///
                xtitle("`axis_label'") ytitle("Weighted arm mass") ///
                title("Local treatment overlap") ///
                legend(order(1 "treated" 2 "control")) `nameopt' `drawopt'
        }
        else {
            graph bar treated_mass control_mass, over(target_id) ///
                yscale(range(0 1)) ylabel(0(.2)1) ///
                yline(`lower' `upper', lpattern(dot) lcolor(gs8)) ///
                ytitle("Weighted arm mass") title("Local treatment overlap") ///
                legend(order(1 "treated" 2 "control")) `nameopt' `drawopt'
        }
        quietly count if treated_share < `lower' | treated_share > `upper'
        if r(N) > 0 di as text "note: " r(N) ///
            " target(s) have treated share outside [`lower', `upper']."
    }
    else if "`view'" == "armess" {
        quietly count if arm_ess_available == 1
        if r(N) == 0 {
            restore
            di as error "Arm-specific ESS is unavailable for a precomputed pseudo-outcome fit."
            exit 498
        }
        graph bar treated_ess control_ess, over(target_id) ///
            ytitle("Normalized within-arm ESS") ///
            title("Local treatment-arm effective sample sizes") ///
            legend(order(1 "treated ESS" 2 "control ESS")) `nameopt' `drawopt'
    }
    else {
        quietly count if propensity_available == 1
        if r(N) == 0 {
            restore
            di as error "Propensity diagnostics are unavailable for a precomputed pseudo-outcome fit."
            exit 498
        }
        twoway (rcap propensity_q10 propensity_q90 `axis', lcolor(navy)) ///
            (scatter propensity_mean `axis', mcolor(navy) msymbol(O)), ///
            yscale(range(0 1)) ylabel(0(.2)1) yline(.5, lpattern(dot) lcolor(gs8)) ///
            xtitle("`axis_label'") ytitle("Weighted propensity (q10-q90)") ///
            title("Local propensity overlap") legend(off) `nameopt' `drawopt'
    }
    restore
    return scalar n_plotted = `nplot'
    return local view "`view'"
end

program define _cate_dnn_plot_weights, rclass
    version 17
    syntax , NAME(name) [VIEW(string) TARGETS(numlist integer >=1) ///
        CANDIDATES(numlist integer >=1) ALLCandidates LOGWeights ///
        REPLACE NODRAW]

    if "`view'" == "" local view "profile"
    local view = lower("`view'")
    if !inlist("`view'", "profile", "cumulative", "kernel") {
        di as error "For type(weights), view() must be profile, cumulative, or kernel."
        exit 198
    }
    if "`view'" != "profile" & "`candidates'`allcandidates'" != "" {
        di as error "Candidate comparisons are available only with view(profile)."
        exit 198
    }
    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"

    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"
    tempname plot meta targets_mat
    local n = e(N)
    local k = e(k)
    local s = e(s)
    local method "`e(method)'"
    local candidate_mode = ("`candidates'`allcandidates'" != "")

    if "`view'" != "kernel" {
        local selected "0"
        if `candidate_mode' {
            if "`e(tuning_mode)'" != "grid" {
                di as error "Candidate weight plots require the active estimator to use tuningmode(grid)."
                exit 498
            }
            matrix `meta' = e(tuning_candidates)
            local gmax = rowsof(`meta')
            local selected "`candidates'"
            if "`allcandidates'" != "" {
                local selected ""
                forvalues g = 1/`gmax' {
                    if el(`meta', `g', 6) == 1 local selected "`selected' `g'"
                }
            }
            foreach g of local selected {
                if `g' > `gmax' | el(`meta', `g', 6) != 1 {
                    di as error "Selected weight candidates must be available rows in 1..`gmax'."
                    exit 498
                }
            }
        }
        if `candidate_mode' mata: st_matrix("`plot'", rdistnn_plot_rank_weights(st_matrix("`meta'"), strtoreal(tokens("`selected'"))', `n', `k', `s', "`method'", 1))
        else mata: st_matrix("`plot'", rdistnn_plot_rank_weights(J(1,8,.), strtoreal(tokens("`selected'"))', `n', `k', `s', "`method'", 0))
        local nselected : word count `selected'
        matrix colnames `plot' = rank weight cum_share candidate_id s is_default
        preserve
        clear
        quietly svmat double `plot', names(col)
        if "`logweights'" != "" quietly drop if weight <= 0
        local plots ""
        local pi = 0
        foreach g of local selected {
            local ++pi
            local ps = mod(`pi' - 1, 15) + 1
            local width "medthin"
            quietly summarize is_default if candidate_id == `g', meanonly
            if r(max) == 1 local width "medthick"
            if "`view'" == "profile" {
                local plots `"`plots' (line weight rank if candidate_id==`g', sort pstyle(p`ps') lwidth(`width'))"'
            }
            else {
                local plots `"`plots' (line cum_share rank if candidate_id==`g', sort pstyle(p`ps') lwidth(`width'))"'
            }
        }
        if "`view'" == "profile" {
            local yscale ""
            if "`logweights'" != "" local yscale "yscale(log)"
            local edgeopt ""
            if !`candidate_mode' {
                local edge = `n' - `s' + `k'
                local edgeopt "xline(`edge', lpattern(dot) lcolor(gs8))"
            }
            twoway `plots', `edgeopt' ///
                `yscale' xtitle("Rank") ytitle("DNN weight") ///
                title("DNN rank-weight profile") `nameopt' `drawopt'
        }
        else {
            twoway `plots', yline(.5 .9, lpattern(dot) lcolor(gs8)) ///
                xtitle("Rank") ytitle("Cumulative weight share") ///
                title("DNN weight concentration") `nameopt' `drawopt'
        }
        restore
        return scalar n_plotted = `nselected'
        return local view "`view'"
        exit
    }

    local xvars "`e(xvars)'"
    foreach v of local xvars {
        capture confirm numeric variable `v'
        if _rc {
            di as error "Current data do not contain numeric covariate `v'; refit `e(cmd)' before plotting kernels."
            exit 498
        }
    }
    tempvar esample
    quietly gen byte `esample' = e(sample)
    quietly count if `esample'
    if r(N) != e(N) {
        di as error "The current estimation sample no longer matches the active `e(cmd)' fit."
        exit 498
    }
    matrix `targets_mat' = e(targets)
    local target_list "`targets'"
    if "`target_list'" == "" {
        if e(m) <= 8 {
            forvalues j = 1/`=e(m)' {
                local target_list "`target_list' `j'"
            }
        }
        else local target_list "1 `=ceil(e(m)/2)' `=e(m)'"
    }
    local ntargets : word count `target_list'
    mata: st_matrix("`plot'", rdistnn_plot_kernel(st_data(., "`xvars'", "`esample'"), st_matrix("`targets_mat'"), strtoreal(tokens("`target_list'"))', `k', `s', "`method'", "`e(distance)'"))
    matrix colnames `plot' = target_id rank distance weight cum_share
    preserve
    clear
    quietly svmat double `plot', names(col)
    twoway (connected weight distance, sort msize(vsmall)), ///
        by(target_id, note("") title("DNN effective kernel")) ///
        xtitle("Distance") ytitle("DNN weight") `nameopt' `drawopt'
    restore
    return scalar n_plotted = `ntargets'
    return local view "kernel"
end

program define _cate_dnn_plot_influence, rclass
    version 17
    syntax , NAME(name) NLABEL(integer) [TARGETS(numlist integer >=1) ///
        REPLACE NODRAW]

    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"
    local xvars "`e(xvars)'"
    foreach v of local xvars {
        capture confirm numeric variable `v'
        if _rc {
            di as error "Current data do not contain numeric covariate `v'; refit cate_dnn before plotting influence."
            exit 498
        }
    }
    local score_method "`e(score_method)'"
    local signal_var "`e(signal)'"
    local outcome_var "`e(outcome)'"
    local treat_var "`e(t_col)'"
    local propensity_var "`e(e_col)'"
    local m0_var "`e(m0_col)'"
    local m1_var "`e(m1_col)'"
    local m_var "`e(m_col)'"
    local needed "`signal_var'"
    if "`score_method'" != "pseudo_outcome" {
        local needed "`outcome_var' `treat_var' `propensity_var'"
        if "`score_method'" == "aipw" local needed "`needed' `m0_var' `m1_var'"
        else local needed "`needed' `m_var'"
    }
    foreach v of local needed {
        capture confirm numeric variable `v'
        if _rc {
            di as error "Current data do not contain numeric score variable `v'; refit cate_dnn before plotting influence."
            exit 498
        }
    }
    tempvar esample rho eclip obsid
    quietly gen byte `esample' = e(sample)
    quietly count if `esample'
    if r(N) != e(N) {
        di as error "The current estimation sample no longer matches the active cate_dnn fit."
        exit 498
    }
    quietly gen long `obsid' = _n
    if "`score_method'" == "pseudo_outcome" {
        quietly gen double `rho' = `signal_var' if `esample'
    }
    else {
        quietly gen double `eclip' = min(max(`propensity_var', e(clip)), 1-e(clip)) if `esample'
        if "`score_method'" == "aipw" {
            quietly gen double `rho' = `treat_var' * (`outcome_var' - `m1_var') / `eclip' ///
                - (1 - `treat_var') * (`outcome_var' - `m0_var') / (1 - `eclip') ///
                + `m1_var' - `m0_var' if `esample'
        }
        else {
            quietly gen double `rho' = (`treat_var' - `eclip') * ///
                (`outcome_var' - `m_var') / (`eclip' * (1 - `eclip')) if `esample'
        }
    }

    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"
    tempname targets_mat plot bmat
    matrix `targets_mat' = e(targets)
    matrix `bmat' = e(b)
    local target_list "`targets'"
    if "`target_list'" == "" {
        if e(m) <= 8 {
            forvalues j = 1/`=e(m)' {
                local target_list "`target_list' `j'"
            }
        }
        else local target_list "1 `=ceil(e(m)/2)' `=e(m)'"
    }
    local ntargets : word count `target_list'
    mata: st_matrix("`plot'", rdistnn_plot_influence(st_data(., "`xvars'", "`esample'"), st_matrix("`targets_mat'"), st_data(., "`rho'", "`esample'"), st_data(., "`obsid'", "`esample'"), st_matrix("`bmat'"), strtoreal(tokens("`target_list'"))', `=e(k)', `=e(s)', "`e(method)'", "`e(distance)'"))
    matrix colnames `plot' = target_id obs_id rank distance weight score contribution
    preserve
    clear
    quietly svmat double `plot', names(col)
    quietly gen double abs_contribution = abs(contribution)
    quietly gsort target_id -abs_contribution
    quietly by target_id: gen str12 obs_label = ///
        string(obs_id, "%9.0g") if _n <= `nlabel' & `nlabel' > 0
    quietly by target_id: gen byte obs_labpos = ///
        cond(_n==1, 12, cond(_n==2, 3, cond(_n==3, 9, 6)))
    quietly sort target_id distance
    twoway (scatter contribution distance, mcolor(black%65) msize(small) ///
        mlabel(obs_label) mlabsize(vsmall) mlabposition(obs_labpos)), ///
        by(target_id, note("") title("CATE influence contributions")) ///
        yline(0, lpattern(dash) lcolor(gs8)) xtitle("Distance") ///
        ytitle("w(i,j) * (score(i) - tau(j))") `nameopt' `drawopt'
    restore
    return scalar n_plotted = `ntargets'
    return local view "influence"
end

program define _cate_dnn_plot_dashboard, rclass
    version 17
    syntax , NAME(name) Level(cilevel) XINDEX(integer) ///
        [STYLE(string) REPLACE NODRAW]

    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"
    local styleopt ""
    if "`style'" != "" local styleopt "style(`style')"
    tempname gest gover gprop gdiag gweight gtune diagcheck
    local graphs ""
    local npanels 0
    local has_overlap 0

    capture quietly _cate_dnn_plot_estimates, name(`gest') replace nodraw ///
        level(`level') xindex(`xindex') `styleopt'
    if !_rc {
        local graphs "`graphs' `gest'"
        local ++npanels
    }
    capture matrix `diagcheck' = e(diagnostics)
    if !_rc {
        if el(`diagcheck', 1, 12) == 1 {
            capture quietly _cate_dnn_plot_overlap, name(`gover') replace nodraw ///
                view(mass) xindex(`xindex')
            if !_rc {
                local graphs "`graphs' `gover'"
                local ++npanels
                local has_overlap 1
            }
        }
        if el(`diagcheck', 1, 16) == 1 {
            capture quietly _cate_dnn_plot_overlap, name(`gprop') replace nodraw ///
                view(propensity) xindex(`xindex')
            if !_rc {
                local graphs "`graphs' `gprop'"
                local ++npanels
                local has_overlap 1
            }
        }
    }
    if !`has_overlap' {
        capture quietly _cate_dnn_plot_diagnostics, name(`gdiag') replace nodraw
        if !_rc {
            local graphs "`graphs' `gdiag'"
            local ++npanels
        }
    }
    capture quietly _cate_dnn_plot_weights, name(`gweight') replace nodraw ///
        view(cumulative)
    if !_rc {
        local graphs "`graphs' `gweight'"
        local ++npanels
    }
    if "`e(tuning_mode)'" == "grid" {
        capture quietly _cate_dnn_plot_tuning, name(`gtune') replace nodraw ///
            level(`level') allcandidates
        if !_rc {
            local graphs "`graphs' `gtune'"
            local ++npanels
        }
    }
    if `npanels' == 0 {
        di as error "No dashboard panel could be constructed from the active fit."
        exit 498
    }
    graph combine `graphs', cols(2) title("DNN CATE dashboard") ///
        note("n=`=e(N)', m=`=e(m)', k=`=e(k)', s=`=e(s)', distance=`e(distance)', jackknife=`e(jackknife)'") ///
        `nameopt' `drawopt'
    return scalar n_plotted = `npanels'
    return local view "dashboard"
end
