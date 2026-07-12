*! version 0.1.0  2026-07-12  Jakob Juergens
*! Postestimation diagnostics tables for np_dnn

program define np_dnn_diagnostics, rclass
    version 17

    syntax [, CANDIDATES(numlist integer >=1) ALLCandidates DETAIL ///
        SAVING(string asis) REPLACE]

    if "`e(cmd)'" != "np_dnn" {
        di as error "np_dnn_diagnostics requires active np_dnn estimates."
        exit 301
    }
    if "`candidates'" != "" & "`allcandidates'" != "" {
        di as error "Specify either candidates() or allcandidates, not both."
        exit 198
    }
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    tempname out compact raw targets candidate_meta candidate_diag
    local scope "targets"

    if "`candidates'`allcandidates'" == "" {
        capture matrix `raw' = e(diagnostics)
        if _rc {
            di as error "Active np_dnn estimates do not contain e(diagnostics); refit with the current package version."
            exit 498
        }
        matrix `targets' = e(targets)
        local m = rowsof(`raw')
        local p = colsof(`targets')
        mata: st_matrix("`out'", rdistnn_plot_target_diag(st_matrix("`raw'"), st_matrix("`targets'")))
        mata: st_matrix("`compact'", st_matrix("`out'")[,(1, `p' + 2, `p' + 3, `p' + 5, `p' + 14, `p' + 15)])
        local x0_names ""
        forvalues j = 1/`p' {
            local x0_names "`x0_names' x0_`j'"
        }
        matrix colnames `out' = target_id `x0_names' ///
            ess max_weight top1_weight_share top5_weight_share ///
            top10_weight_share n_nonzero_weight distance_available ///
            radius_rank_k radius_rank_s weighted_mean_distance ///
            weighted_median_distance outcome_available weighted_mean_y ///
            weighted_sd_y weighted_y_q10 weighted_y_q50 weighted_y_q90
        matrix rownames `out' = `: rownames `raw''
        matrix colnames `compact' = target_id ess max_weight ///
            top5_weight_share weighted_mean_y weighted_sd_y
        matrix rownames `compact' = `: rownames `raw''
    }
    else {
        if "`e(tuning_mode)'" != "grid" {
            di as error "Candidate diagnostics require np_dnn, tuningmode(grid)."
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
    }

    if "`detail'" != "" matrix list `out', format(%12.6g)
    else matrix list `compact', format(%12.6g)

    if `"`saving'"' != "" {
        preserve
        clear
        quietly svmat double `out', names(col)
        if "`replace'" != "" save `saving', replace
        else save `saving'
        restore
    }

    return matrix diagnostics = `out'
    return matrix compact = `compact'
    return scalar n_rows = `m'
    return local scope "`scope'"
    return local saved `saving'
end
