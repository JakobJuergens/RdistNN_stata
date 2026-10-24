*! version 0.1.0  2026-10-14  Jakob Juergens
// Shared CATE/CLATE/NP plot builders, loaded explicitly to give both commands access.
// Keep these routines outside public ado-local subroutine namespaces.
// Loading is guarded by program list in each public plot command.

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

    tempname targets_mat bmat vmat semat plotmat candidate_meta grid_est grid_se statusmask
    matrix `targets_mat' = e(targets)
    local m = e(m)
    local p = e(p)
    local xvars "`e(xvars)'"
    local z = invnormal(1 - (1 - `level'/100) / 2)
    local use_candidates = ("`candidates'`allcandidates'" != "")
    local estimand "CATE"
    local plot_title "DNN CATE estimates"
    if "`e(cmd)'" == "clate_dnn" {
        local estimand "CLATE"
        local plot_title "DNN CLATE estimates"
    }
    if inlist("`e(cmd)'", "np_dnn", "np_tdnn") {
        local estimand "Conditional mean"
        local plot_title "DNN nonparametric regression estimates"
        if "`e(cmd)'"=="np_tdnn" local plot_title "Two-scale regression estimates"
    }

    if !`use_candidates' {
        local emats : e(matrices)
        foreach field in estimates se_targets {
            if !`: list field in emats' {
                di as error "Current fitted results require e(`field'); refit `e(cmd)' with the current package."
                exit 498
            }
        }
        matrix `bmat' = e(estimates)'
        matrix `semat' = e(se_targets)'
        if "`e(jackknife)'"=="none" | "`e(se_method)'"=="none" matrix `semat' = J(`m',1,.)
        if "`e(cmd)'"=="clate_dnn" {
            matrix `statusmask' = e(status)
        }
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
        if "`e(cmd)'" == "clate_dnn" {
            matrix `statusmask' = e(tuning_target_status)
        }
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
    if "`e(cmd)'"=="clate_dnn" {
        local cidx 0
        if `use_candidates' local cidx = `p'+1
        _rdistnn_plot_wald_bounds, matrix(`plotmat') status(`statusmask') ///
            target(`=`p'+4') lower(`=`p'+7') upper(`=`p'+8') candidate(`cidx')
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

    // Keep the confidence fill visible in StataNow 19; lcolor(none)
    // collapses rarea polygons in SVG and raster export on that version.
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
            capture assert missing(ci_lower) | missing(ci_upper) if candidate_id == `g'
            if _rc {
                local plots `"`plots' (rarea ci_lower ci_upper `axis' if candidate_id==`g', sort color(gs12%45))"'
                local ++layercount
            }
            local plots `"`plots' (line estimate `axis' if candidate_id==`g', sort pstyle(p`ps') lwidth(`width'))"'
            local ++layercount
        }
        else if "`style'" == "forest" {
            quietly gen double __y`plotnum' = target_id + ///
                (`plotnum' - (`ng' + 1)/2) * .12
            capture assert missing(ci_lower) | missing(ci_upper) if candidate_id == `g'
            if _rc {
                local plots `"`plots' (rcap ci_lower ci_upper __y`plotnum' if candidate_id==`g', horizontal pstyle(p`ps'))"'
                local ++layercount
            }
            local plots `"`plots' (scatter __y`plotnum' estimate if candidate_id==`g', msymbol(`symbol') pstyle(p`ps'))"'
            local ++layercount
        }
        else {
            capture assert missing(ci_lower) | missing(ci_upper) if candidate_id == `g'
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
    if "`e(cmd)'"=="np_tdnn" local subtitle "n = `=e(N)', s1 = `=e(s1)', s2 = `=e(s2)', dim = `=e(dim)'"
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
    quietly keep if !missing(estimate)
    quietly count
    local nrows = r(N)
    if `nrows'>0 quietly mkmat `xnames' candidate_id s ///
        is_default target_id estimate se ci_lower ci_upper, matrix(`plotmat')
    restore
    return scalar n_rows = `nrows'
    return scalar n_omitted = `nplot'-`nrows'
    return local selected_targets "`target_levels'"
    if `nrows'>0 return matrix plot_data = `plotmat'
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
    tempname meta est ses plot statusmask
    matrix `meta' = e(tuning_candidates)
    matrix `est' = e(tuning_estimates)
    matrix `ses' = e(tuning_se)
    if "`e(cmd)'" == "clate_dnn" {
        matrix `statusmask' = e(tuning_target_status)
    }
    local gmax = rowsof(`meta')
    local m = e(m)
    local z = invnormal(1 - (1 - `level'/100) / 2)
    local estimand "CATE"
    local plot_title "CATE tuning sensitivity"
    if "`e(cmd)'" == "clate_dnn" {
        local estimand "CLATE"
        local plot_title "CLATE tuning sensitivity"
    }
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
    if "`e(cmd)'"=="clate_dnn" _rdistnn_plot_wald_bounds, ///
        matrix(`plot') status(`statusmask') target(5) lower(8) upper(9) candidate(1)
    matrix colnames `plot' = candidate_id s lambda is_default target_id ///
        estimate se ci_lower ci_upper
    preserve
    clear
    quietly svmat double `plot', names(col)
    capture assert missing(ci_lower) | missing(ci_upper)
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
    quietly keep if !missing(estimate)
    quietly count
    local nrows = r(N)
    if `nrows'>0 quietly mkmat candidate_id s lambda is_default target_id estimate se ci_lower ci_upper, matrix(`plot')
    restore
    return scalar n_rows = `nrows'
    local nt : word count `target_list'
    return scalar n_omitted = `ng'*`nt'-`nrows'
    return local selected_targets "`target_list'"
    if `nrows'>0 return matrix plot_data = `plot'
    return local view "tuning"
end

program define _cate_dnn_plot_diagnostics, rclass
    version 17
    syntax , NAME(name) [TARGETS(numlist integer >=1) CANDIDATES(numlist integer >=1) ///
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
        if inlist("`e(cmd)'", "cate_dnn", "clate_dnn") {
            mata: st_matrix("`plot'", ((1::rows(st_matrix("`raw'"))), st_matrix("`raw'")[,(1,2,4,20,25)], J(rows(st_matrix("`raw'")),1,0)))
            matrix colnames `plot' = plot_id ess max_weight top5_weight_share ///
                physical_radius_90 top5_score_share is_default
        }
        else {
            mata: st_matrix("`plot'", ((1::rows(st_matrix("`raw'"))), st_matrix("`raw'")[,(1,2,4)], J(rows(st_matrix("`raw'")),1,0)))
            matrix colnames `plot' = plot_id ess max_weight top5_weight_share is_default
        }
        if "`targets'" != "" {
            foreach t of local targets {
                if `t'>`m' error 198
            }
            local plotcols : colnames `plot'
            mata: st_matrix("`plot'", st_matrix("`plot'")[strtoreal(tokens("`targets'"))',])
            matrix colnames `plot' = `plotcols'
            local m : word count `targets'
        }
        local xvar "plot_id"
        local xlabel "Target"
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
        local guide = .9 * e(N)
    }

    local oldnames : colnames `plot'
    if !`candidate_mode' {
        mata: st_matrix("`plot'", (st_matrix("`plot'"), st_matrix("`plot'")[,1], J(rows(st_matrix("`plot'")),1,0)))
    }
    else {
        mata: st_matrix("`plot'", (st_matrix("`plot'"), J(rows(st_matrix("`plot'")),1,.), strtoreal(tokens("`selected'"))'))
    }
    matrix colnames `plot' = `oldnames' target_id candidate_id
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
    if "`e(cmd)'" == "clate_dnn" local plot_title "DNN CLATE diagnostics"
    if "`e(cmd)'" == "np_dnn" local plot_title "DNN nonparametric-regression diagnostics"
    graph combine `g1' `g2' `g3', cols(3) ///
        title("`plot_title'") `nameopt' `drawopt'
    if !`candidate_mode' & inlist("`e(cmd)'", "cate_dnn", "clate_dnn") {
        quietly summarize physical_radius_90, meanonly
        di as text "90% physical-radius range: " as result %9.4g r(min) ///
            as text " to " as result %9.4g r(max)
        quietly summarize top5_score_share, meanonly
        di as text "Top-5 absolute score-contribution share range: " ///
            as result %6.3f r(min) as text " to " as result %6.3f r(max)
    }
    restore
    local nrows = rowsof(`plot')
    return matrix plot_data = `plot'
    return scalar n_rows = `nrows'
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

    local armrole "treatment"
    local propensity_title "Local propensity overlap"
    local arm1 "treated"
    local arm0 "control"
    if "`e(cmd)'" == "clate_dnn" {
        local armrole "instrument"
        local propensity_title "Local instrument propensity overlap"
        local arm1 "instrument = 1"
        local arm0 "instrument = 0"
    }
    tempname diag targets_mat plot
    matrix `diag' = e(diagnostics)
    matrix `targets_mat' = e(targets)
    local m = e(m)
    local p = e(p)
    local xvars "`e(xvars)'"
    local xnames ""
    forvalues j=1/`p' {
        local xnames "`xnames' x0_`j'"
    }
    if "`view'" == "mass" {
        mata: st_matrix("`plot'", rdistnn_plot_overlap(st_matrix("`diag'"), st_matrix("`targets_mat'"), (12,13,14,15)))
        matrix colnames `plot' = `xnames' target_id ///
            overlap_available treated_mass control_mass treated_share
    }
    else if "`view'" == "armess" {
        mata: st_matrix("`plot'", rdistnn_plot_overlap(st_matrix("`diag'"), st_matrix("`targets_mat'"), (21,22,23)))
        matrix colnames `plot' = `xnames' target_id ///
            arm_ess_available treated_ess control_ess
    }
    else {
        mata: st_matrix("`plot'", rdistnn_plot_overlap(st_matrix("`diag'"), st_matrix("`targets_mat'"), (16,17,18,19)))
        matrix colnames `plot' = `xnames' target_id ///
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
    local available "overlap_available"
    if "`view'"=="armess" local available "arm_ess_available"
    if "`view'"=="propensity" local available "propensity_available"
    quietly keep if `available'==1
    quietly count
    if r(N)==0 {
        restore
        di as error "Selected overlap diagnostics are unavailable for this fit."
        exit 498
    }
    local nplot = r(N)

    local axis "target_id"
    local axis_label "Target"
    local varying ""
    if `xindex'>0 {
        local axis "x0_`xindex'"
        local axis_label : word `xindex' of `xvars'
    }
    else {
        forvalues j=1/`p' {
            quietly levelsof x0_`j', local(vals)
            local nv : word count `vals'
            if `nv'>1 local varying "`varying' `j'"
        }
        local nvary : word count `varying'
        if `nvary'==1 {
            local coordinate : word 1 of `varying'
            local axis "x0_`coordinate'"
            local axis_label : word `coordinate' of `xvars'
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
                title("Local `armrole' overlap") ///
                legend(order(1 "`arm1'" 2 "`arm0'")) `nameopt' `drawopt'
        }
        else {
            graph bar treated_mass control_mass, over(target_id) ///
                yscale(range(0 1)) ylabel(0(.2)1) ///
                yline(`lower' `upper', lpattern(dot) lcolor(gs8)) ///
                ytitle("Weighted arm mass") title("Local `armrole' overlap") ///
                legend(order(1 "`arm1'" 2 "`arm0'")) `nameopt' `drawopt'
        }
        quietly count if treated_share < `lower' | treated_share > `upper'
        if r(N) > 0 di as text "note: " r(N) ///
            " target(s) have `arm1' share outside [`lower', `upper']."
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
            title("Local `armrole'-arm effective sample sizes") ///
            legend(order(1 "`arm1' ESS" 2 "`arm0' ESS")) `nameopt' `drawopt'
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
            title("`propensity_title'") legend(off) `nameopt' `drawopt'
    }
    local cols : colnames `plot'
    quietly mkmat `cols', matrix(`plot')
    if "`e(cmd)'"=="clate_dnn" {
        local cols : subinstr local cols "treated_mass" "instrument_1_mass", all
        local cols : subinstr local cols "control_mass" "instrument_0_mass", all
        local cols : subinstr local cols "treated_share" "instrument_1_share", all
        local cols : subinstr local cols "treated_ess" "instrument_1_ess", all
        local cols : subinstr local cols "control_ess" "instrument_0_ess", all
        matrix colnames `plot' = `cols'
    }
    restore
    local nrows = rowsof(`plot')
    return matrix plot_data = `plot'
    return scalar n_rows = `nrows'
    return local view "`view'"
end

program define _cate_dnn_plot_weights, rclass
    version 17
    syntax , NAME(name) [VIEW(string) TARGETS(numlist integer >=1) ///
        CANDIDATES(numlist integer >=1) ALLCandidates LOGWeights ///
        REPLACE NODRAW VALIDATED]

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

    if e(model_version)!=4 {
        di as error "Weight reconstruction requires a current saved model; refit with the current package."
        exit 498
    }
    findfile "_rdistnn_model_decode.ado"
    quietly run "`r(fn)'"
    if "`validated'"=="" mata: rdistnn_model_validate()
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
        if `candidate_mode' noisily mata: st_matrix("`plot'", rdistnn_plot_rank_weights(st_matrix("`meta'"), strtoreal(tokens("`selected'"))', `n', `k', `s', "`method'", 1,st_matrix("e(model_ranks)"),rdistnn_model_target_index(strtoreal(tokens("`targets'")),st_numscalar("e(m)")),st_matrix("e(model_distances)")))
        else noisily mata: st_matrix("`plot'", rdistnn_plot_rank_weights(J(1,8,.), strtoreal(tokens("`selected'"))', `n', `k', `s', "`method'", 0,st_matrix("e(model_ranks)"),rdistnn_model_target_index(strtoreal(tokens("`targets'")),st_numscalar("e(m)")),st_matrix("e(model_distances)")))
        matrix colnames `plot' = rank weight cum_share candidate_id s is_default target_id
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
            if !`candidate_mode' & "`method'"=="exact" {
                local edge = `n' - `s' + `k'
                local edgeopt "xline(`edge', lpattern(dot) lcolor(gs8)) note(Nominal exact-support edge: `edge')"
            }
            twoway `plots', `edgeopt' by(target_id, note("")) ///
                `yscale' xtitle("Rank") ytitle("DNN weight") ///
                title("DNN rank-weight profile") `nameopt' `drawopt'
        }
        else {
            twoway `plots', by(target_id, note("")) yline(.5 .9, lpattern(dot) lcolor(gs8)) ///
                xtitle("Rank") ytitle("Cumulative weight share") ///
                title("DNN weight concentration") `nameopt' `drawopt'
        }
        quietly mkmat rank weight cum_share candidate_id s is_default target_id, matrix(`plot')
        restore
        local nrows = rowsof(`plot')
        return matrix plot_data = `plot'
        return scalar n_rows = `nrows'
        return local view "`view'"
        exit
    }

    _rdistnn_plot_sample_guard
    _rdistnn_model_plot, type(kernel) name(`name') `replace' `nodraw' targets(`targets') validated
    return add
end

program define _cate_dnn_plot_influence, rclass
    version 17
    syntax , NAME(name) NLABEL(integer) [TARGETS(numlist integer >=1) ///
        REPLACE NODRAW VALIDATED]
    _rdistnn_plot_sample_guard
    _rdistnn_model_plot, type(influence) name(`name') nlabel(`nlabel') `replace' `nodraw' targets(`targets') `validated'
    return add
end

program define _cate_dnn_plot_dashboard, rclass
    version 17
    syntax , NAME(name) Level(cilevel) XINDEX(integer) ///
        [STYLE(string) REPLACE NODRAW VALIDATED]
    if "`validated'"=="" _rdistnn_model_check

    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"
    local styleopt ""
    if "`style'" != "" local styleopt "style(`style')"
    tempname gest gover gprop gdiag gweight gtune diagcheck
    local graphs ""
    local panels ""
    tempname panel_info
    local npanels 0
    local has_overlap 0

    quietly _cate_dnn_plot_estimates, name(`gest') replace nodraw ///
        level(`level') xindex(`xindex') `styleopt'
    if !_rc {
        _rdistnn_plot_add_panel, matrix(`panel_info') id(`=`npanels'+1')
                local panels "`panels' estimates"
                local graphs "`graphs' `gest'"
        local ++npanels
    }
    matrix `diagcheck' = e(diagnostics)
    if !_rc {
        if el(`diagcheck', 1, 12) == 1 {
            quietly _cate_dnn_plot_overlap, name(`gover') replace nodraw ///
                view(mass) xindex(`xindex')
            if !_rc {
                _rdistnn_plot_add_panel, matrix(`panel_info') id(`=`npanels'+1')
                local panels "`panels' overlap_mass"
                local graphs "`graphs' `gover'"
                local ++npanels
                local has_overlap 1
            }
        }
        if el(`diagcheck', 1, 16) == 1 {
            quietly _cate_dnn_plot_overlap, name(`gprop') replace nodraw ///
                view(propensity) xindex(`xindex')
            if !_rc {
                _rdistnn_plot_add_panel, matrix(`panel_info') id(`=`npanels'+1')
                local panels "`panels' overlap_propensity"
                local graphs "`graphs' `gprop'"
                local ++npanels
                local has_overlap 1
            }
        }
    }
    if !`has_overlap' {
        quietly _cate_dnn_plot_diagnostics, name(`gdiag') replace nodraw
        if !_rc {
            _rdistnn_plot_add_panel, matrix(`panel_info') id(`=`npanels'+1')
                local panels "`panels' diagnostics"
                local graphs "`graphs' `gdiag'"
            local ++npanels
        }
    }
    quietly _cate_dnn_plot_weights, name(`gweight') replace nodraw ///
        view(cumulative) validated
    if !_rc {
        _rdistnn_plot_add_panel, matrix(`panel_info') id(`=`npanels'+1')
                local panels "`panels' weights_cumulative"
                local graphs "`graphs' `gweight'"
        local ++npanels
    }
    if "`e(tuning_mode)'" == "grid" {
        quietly _cate_dnn_plot_tuning, name(`gtune') replace nodraw ///
            level(`level') allcandidates
        if !_rc {
            _rdistnn_plot_add_panel, matrix(`panel_info') id(`=`npanels'+1')
                local panels "`panels' tuning"
                local graphs "`graphs' `gtune'"
            local ++npanels
        }
    }
    if `npanels' == 0 {
        di as error "No dashboard panel could be constructed from the active fit."
        exit 498
    }
    local dashboard_title "DNN CATE dashboard"
    if "`e(cmd)'" == "clate_dnn" local dashboard_title "DNN CLATE dashboard"
    graph combine `graphs', cols(2) title("`dashboard_title'") ///
        note("n=`=e(N)', m=`=e(m)', k=`=e(k)', s=`=e(s)', distance=`e(distance)', jackknife=`e(jackknife)'") ///
        `nameopt' `drawopt'
    return matrix panel_info = `panel_info'
    return local panels "`panels'"
    return scalar n_panels = `npanels'
    return local view "dashboard"
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
                color(gs12%45)) ///
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
    quietly mkmat `xnames' target_id weighted_mean_y weighted_sd_y ///
        weighted_y_q10 weighted_y_q50 weighted_y_q90, matrix(`plot')
    restore
    local nrows = rowsof(`plot')
    return matrix plot_data = `plot'
    return scalar n_rows = `nrows'
    return local view "`view'"
end

program define _np_dnn_plot_dashboard, rclass
    version 17
    syntax , NAME(name) Level(cilevel) XINDEX(integer) ///
        [STYLE(string) REPLACE NODRAW VALIDATED]
    if "`validated'"=="" _rdistnn_model_check

    local nameopt "name(`name')"
    if "`replace'" != "" local nameopt "name(`name', replace)"
    local drawopt ""
    if "`nodraw'" != "" local drawopt "nodraw"
    local style_opt ""
    if "`style'" != "" local style_opt "style(`style')"
    tempname gest gdiag gout gweights gtune panel_info

    quietly _cate_dnn_plot_estimates, name(`gest') replace nodraw ///
        level(`level') xindex(`xindex') `style_opt'
    _rdistnn_plot_add_panel, matrix(`panel_info') id(1)
    quietly _cate_dnn_plot_diagnostics, name(`gdiag') replace nodraw
    _rdistnn_plot_add_panel, matrix(`panel_info') id(2)
    quietly _np_dnn_plot_outcome, name(`gout') replace nodraw xindex(`xindex')
    _rdistnn_plot_add_panel, matrix(`panel_info') id(3)
    quietly _cate_dnn_plot_weights, name(`gweights') replace nodraw view(cumulative) validated
    _rdistnn_plot_add_panel, matrix(`panel_info') id(4)
    local graphs "`gest' `gdiag' `gout' `gweights'"
    local panels "estimates diagnostics outcome weights_cumulative"
    local npanels 4
    if "`e(tuning_mode)'" == "grid" {
        quietly _cate_dnn_plot_tuning, name(`gtune') replace nodraw level(`level')
        _rdistnn_plot_add_panel, matrix(`panel_info') id(5)
        local panels "`panels' tuning"
        local graphs "`graphs' `gtune'"
        local npanels 5
    }
    graph combine `graphs', cols(2) title("DNN nonparametric-regression dashboard") ///
        `nameopt' `drawopt'
    return matrix panel_info = `panel_info'
    return local panels "`panels'"
    return scalar n_panels = `npanels'
    return local view "dashboard"
end

// R retains the stored CLATE SE even when a status rules out Wald bounds.
program define _rdistnn_plot_wald_bounds
    version 17
    syntax , MATRIX(name) STATUS(name) TARGET(integer) LOWER(integer) ///
        UPPER(integer) [CANDIDATE(integer 0)]
    tempname work statuses available
    matrix `work' = `matrix'
    matrix `statuses' = `status'
    mata: st_matrix("`available'",rdistnn_plot_wald_available(st_matrix("`work'"),st_matrix("`statuses'"),`target',`candidate'))
    forvalues i=1/`=rowsof(`work')' {
        if el(`available',`i',1)!=1 {
            matrix `matrix'[`i',`lower'] = .
            matrix `matrix'[`i',`upper'] = .
        }
    }
end
