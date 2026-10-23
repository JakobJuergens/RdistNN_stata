*! version 0.1.0  2026-10-17  Jakob Juergens
*! Shared target/candidate selection and localization diagnostics

program define _rdistnn_diagnostics, rclass
    version 17
    syntax, COMMAND(name) [TARGETS(numlist integer >=1) ///
        CANDIDATES(numlist integer >=1) ALLCandidates]

    if "`e(cmd)'" != "`command'" {
        di as error "`command'_diagnostics requires active `command' estimates."
        exit 301
    }
    if "`candidates'" != "" & "`allcandidates'" != "" {
        di as error "Specify either candidates() or allcandidates, not both."
        exit 198
    }
    if "`targets'" != "" & "`candidates'`allcandidates'" != "" {
        di as error "targets() cannot be combined with candidate diagnostics."
        exit 198
    }

    _rdistnn_model_check
    tempname raw profiles out compact localization
    local stored : e(matrices)
    local scope "targets"
    if "`candidates'`allcandidates'" == "" {
        matrix `raw' = e(diagnostics)
        matrix `profiles' = e(targets)
        local total = rowsof(`raw')
        local p = colsof(`profiles')
        local selected "`targets'"
        if "`selected'" == "" {
            forvalues j = 1/`total' {
                local selected "`selected' `j'"
            }
        }
        local source_rows : rownames `raw'
        local selected_rows ""
        foreach j of local selected {
            if `j' > `total' {
                di as error "targets() entries must be in 1..`total'."
                exit 198
            }
            local rowname : word `j' of `source_rows'
            local selected_rows "`selected_rows' `rowname'"
        }
        mata: st_matrix("`out'", (strtoreal(tokens("`selected'"))', st_matrix("`profiles'")[strtoreal(tokens("`selected'"))',], st_matrix("`raw'")[strtoreal(tokens("`selected'"))',]))
        local x0_names ""
        forvalues j = 1/`p' {
            local x0_names "`x0_names' x0_`j'"
        }
        matrix colnames `out' = target_id `x0_names' `: colnames `raw''
        matrix rownames `out' = `selected_rows'
        matrix `localization' = `out'[1...,1..`=`p'+12']
        matrix colnames `localization' = target_id `x0_names' ///
            ess max_weight top1_weight_share top5_weight_share ///
            top10_weight_share n_nonzero_weight distance_available ///
            radius_rank_k radius_rank_s weighted_mean_distance ///
            weighted_median_distance
        return matrix localization = `localization'
        return scalar n_covariates = `p'
    }
    else {
        if "`e(tuning_mode)'" != "grid" {
            di as error "Candidate diagnostics require `command', tuningmode(grid)."
            exit 498
        }
        foreach required in tuning_candidates tuning_diagnostics {
            if !`: list required in stored' {
                di as error "Active estimates do not contain e(`required')."
                exit 498
            }
        }
        matrix `profiles' = e(tuning_candidates)
        matrix `raw' = e(tuning_diagnostics)
        local total = rowsof(`profiles')
        if rowsof(`raw') != `total' | colsof(`profiles') != 8 | colsof(`raw') != 5 {
            di as error "Active estimates contain incompatible candidate diagnostics."
            exit 498
        }
        local selected "`candidates'"
        if "`allcandidates'" != "" {
            forvalues j = 1/`total' {
                local selected "`selected' `j'"
            }
        }
        local selected_rows ""
        foreach j of local selected {
            if `j' > `total' {
                di as error "candidates() entries must be in 1..`total'."
                exit 198
            }
            local selected_rows "`selected_rows' candidate`j'"
        }
        mata: st_matrix("`out'", (st_matrix("`profiles'")[strtoreal(tokens("`selected'"))',], st_matrix("`raw'")[strtoreal(tokens("`selected'"))',]))
        matrix colnames `out' = candidate_id computed_id is_default ///
            feasible reason status lambda s ess max_weight top_n ///
            top5_weight_share degenerate
        matrix rownames `out' = `selected_rows'
        mata: st_matrix("`compact'", st_matrix("`out'")[,(1,3,4,6,7,8,9,10,12,13)])
        matrix colnames `compact' = candidate_id is_default feasible ///
            status lambda s ess max_weight top5_weight_share degenerate
        matrix rownames `compact' = `selected_rows'
        return matrix compact = `compact'
        local scope "candidates"
    }
    local selected : list retokenize selected
    return matrix diagnostics = `out'
    return scalar n_rows = `: word count `selected''
    return local selected "`selected'"
    return local scope "`scope'"
end
