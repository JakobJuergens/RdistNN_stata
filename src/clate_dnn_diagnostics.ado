*! version 0.1.0  2026-10-17  Jakob Juergens
*! Postestimation diagnostics tables for clate_dnn

program define clate_dnn_diagnostics, rclass
    version 17

    syntax [, TARGETS(numlist integer >=1) CANDIDATES(numlist integer >=1) ALLCandidates DETAIL ///
        SAVING(string asis) REPLACE]

    if "`e(cmd)'" != "clate_dnn" {
        di as error "clate_dnn_diagnostics requires active clate_dnn estimates."
        exit 301
    }
    if "`candidates'" != "" & "`allcandidates'" != "" {
        di as error "Specify either candidates() or allcandidates, not both."
        exit 198
    }
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    if "`targets'" != "" & "`candidates'`allcandidates'" != "" {
        di as error "targets() cannot be combined with candidate diagnostics."
        exit 198
    }
    local target_selection "`targets'"
    tempname out compact raw targets candidate_meta candidate_diag ///
        components status component_raw component_vcov estimates report
    local scope "targets"

    if "`candidates'`allcandidates'" == "" {
        capture matrix `raw' = e(diagnostics)
        if _rc {
            di as error "Active clate_dnn estimates do not contain e(diagnostics); refit with the current package version."
            exit 498
        }
        matrix `targets' = e(targets)
        local m = rowsof(`raw')
        local p = colsof(`targets')
        mata: st_matrix("`out'", rdistnn_plot_target_diag(st_matrix("`raw'"), st_matrix("`targets'")))
        mata: st_matrix("`compact'", st_matrix("`out'")[,(1, `p' + 2, `p' + 5, `p' + 21, `p' + 23, `p' + 24, `p' + 26)])
        local x0_names ""
        forvalues j = 1/`p' {
            local x0_names "`x0_names' x0_`j'"
        }
        matrix colnames `out' = target_id `x0_names' ///
            ess max_weight top1_weight_share top5_weight_share ///
            top10_weight_share n_nonzero_weight distance_available ///
            radius_rank_k radius_rank_s weighted_mean_distance ///
            weighted_median_distance overlap_available instrument_1_mass ///
            instrument_0_mass instrument_1_share propensity_available ///
            weighted_mean_propensity weighted_propensity_q10 ///
            weighted_propensity_q90 physical_radius_90 ///
            arm_ess_available instrument_1_ess instrument_0_ess ///
            score_contribution_available top5_score_share
        matrix rownames `out' = `: rownames `raw''
        matrix colnames `compact' = target_id ess top5_weight_share ///
            physical_radius_90 instrument_1_ess instrument_0_ess top5_score_share
        matrix rownames `compact' = `: rownames `raw''
        matrix `component_raw' = e(components)
        matrix `component_vcov' = e(component_vcov)
        matrix `estimates' = e(estimates)
        matrix `status' = e(status)
        matrix `components' = J(`m',9,.)
        forvalues j=1/`m' {
            matrix `components'[`j',1] = `j'
            matrix `components'[`j',2] = el(`component_raw',`j',1)
            matrix `components'[`j',3] = el(`component_raw',`j',2)
            matrix `components'[`j',4] = sqrt(el(`component_vcov',2*`j'-1,2*`j'-1))
            matrix `components'[`j',5] = sqrt(el(`component_vcov',2*`j',2*`j'))
            matrix `components'[`j',6] = el(`component_vcov',2*`j'-1,2*`j')
            if el(`components',`j',4)>0 & el(`components',`j',4)<. & ///
                el(`components',`j',5)>0 & el(`components',`j',5)<. ///
                matrix `components'[`j',7] = el(`components',`j',6) / ///
                    (el(`components',`j',4)*el(`components',`j',5))
            if el(`components',`j',5)>0 & el(`components',`j',5)<. ///
                matrix `components'[`j',8] = el(`components',`j',3)/el(`components',`j',5)
            matrix `components'[`j',9] = el(`estimates',1,`j')
        }
        matrix colnames `components' = target_id reduced_form first_stage ///
            numerator_se denominator_se component_covariance component_correlation ///
            denominator_stat estimate
        matrix rownames `components' = `: rownames `raw''
        if "`target_selection'" != "" {
            foreach mat in out compact components status {
                local cn_`mat' : colnames ``mat''
            }
            foreach j of local target_selection {
                if `j'>`m' {
                    di as error "targets() entries must be in 1..`m'."
                    exit 198
                }
            }
            mata: st_matrix("`out'", st_matrix("`out'")[strtoreal(tokens("`target_selection'"))',]); st_matrix("`compact'", st_matrix("`compact'")[strtoreal(tokens("`target_selection'"))',]); st_matrix("`components'", st_matrix("`components'")[strtoreal(tokens("`target_selection'"))',]); st_matrix("`status'", st_matrix("`status'")[strtoreal(tokens("`target_selection'"))',])
            local selected_rows ""
            foreach j of local target_selection {
                local selected_rows "`selected_rows' target`j'"
            }
            foreach mat in out compact components status {
                matrix rownames ``mat'' = `selected_rows'
                matrix colnames ``mat'' = `cn_`mat''
            }
            local m : word count `target_selection'
        }
        matrix `report' = `out', `components'[1...,2...], `status'
    }
    else {
        if "`e(tuning_mode)'" != "grid" {
            di as error "Candidate diagnostics require clate_dnn, tuningmode(grid)."
            exit 498
        }
        capture matrix `candidate_meta' = e(tuning_candidates)
        if _rc {
            di as error "Active estimates do not contain e(tuning_candidates)."
            exit 498
        }
        matrix `candidate_diag' = e(tuning_diagnostics)
        local gmax = rowsof(`candidate_meta')
        local selected "`candidates'"
        if "`allcandidates'" != "" {
            local selected ""
            forvalues g = 1/`gmax' {
                local selected "`selected' `g'"
            }
        }
        foreach g of local selected {
            if `g' > `gmax' {
                di as error "candidates() entries must be in 1..`gmax'."
                exit 198
            }
        }
        local nselected : word count `selected'
        mata: st_matrix("`out'", rdistnn_plot_candidate_diag(st_matrix("`candidate_meta'"), st_matrix("`candidate_diag'"), strtoreal(tokens("`selected'"))'))
        mata: st_matrix("`compact'", st_matrix("`out'")[,(1, 3, 4, 6, 7, 8, 9, 10, 12, 13)])
        matrix colnames `out' = candidate_id computed_id is_default ///
            feasible reason status lambda s ess max_weight top_n ///
            top5_weight_share degenerate
        local candidate_rows ""
        foreach g of local selected {
            local candidate_rows "`candidate_rows' candidate`g'"
        }
        matrix rownames `out' = `candidate_rows'
        matrix colnames `compact' = candidate_id is_default feasible ///
            status lambda s ess max_weight top5_weight_share degenerate
        matrix rownames `compact' = `candidate_rows'
        local scope "candidates"
        local m = `nselected'
        matrix `report' = `out'
    }

    if "`detail'" != "" matrix list `out', format(%12.6g)
    else matrix list `compact', format(%12.6g)
    if "`scope'" == "targets" {
        matrix list `components', format(%12.6g)
        matrix list `status', format(%9.0g)
        di as text "Arms refer to the instrument A; these are pointwise Wald diagnostics, without weak-IV-uniform guarantees."
    }

    if `"`saving'"' != "" {
        preserve
        clear
        quietly svmat double `report', names(col)
        if "`replace'" != "" save `saving', replace
        else save `saving'
        restore
    }

    if "`scope'" == "targets" {
        return matrix components = `components'
        return matrix status = `status'
    }
    return matrix diagnostics = `out'
    return matrix compact = `compact'
    return scalar n_rows = `m'
    return local scope "`scope'"
    return local saved `saving'
end
