*! version 0.1.2  2026-10-17  Jakob Juergens
*! Postestimation diagnostics tables for cate_dnn

program define cate_dnn_diagnostics, rclass
    version 17
    syntax [, TARGETS(numlist integer >=1) CANDIDATES(numlist integer >=1) ///
        ALLCandidates DETAIL SAVING(string asis) REPLACE]

    local targetopt ""
    local candidateopt ""
    if "`targets'" != "" local targetopt "targets(`targets')"
    if "`candidates'" != "" local candidateopt "candidates(`candidates')"
    _rdistnn_diagnostics, command(cate_dnn) `targetopt' `candidateopt' `allcandidates'
    tempname out compact localization report
    matrix `out' = r(diagnostics)
    local scope "`r(scope)'"
    local selected "`r(selected)'"
    local m = r(n_rows)
    local p = r(n_covariates)
    if "`scope'" == "targets" {
        matrix `localization' = r(localization)
        local x0_names ""
        forvalues j = 1/`p' {
            local x0_names "`x0_names' x0_`j'"
        }
        matrix colnames `out' = target_id `x0_names' ///
            ess max_weight top1_weight_share top5_weight_share ///
            top10_weight_share n_nonzero_weight distance_available ///
            radius_rank_k radius_rank_s weighted_mean_distance ///
            weighted_median_distance overlap_available treated_mass ///
            control_mass treated_share propensity_available ///
            weighted_mean_propensity weighted_propensity_q10 ///
            weighted_propensity_q90 physical_radius_90 ///
            arm_ess_available treated_ess control_ess ///
            score_contribution_available top5_score_share
        mata: st_matrix("`compact'", st_matrix("`out'")[,(1, `p'+2, `p'+5, `p'+21, `p'+23, `p'+24, `p'+26)])
        matrix colnames `compact' = target_id ess top5_weight_share ///
            physical_radius_90 treated_ess control_ess top5_score_share
        matrix rownames `compact' = `: rownames `out''
    }
    else matrix `compact' = r(compact)

    if "`detail'" != "" matrix list `out', format(%12.6g)
    else matrix list `compact', format(%12.6g)

    matrix `report' = `out'
    if `"`saving'"' != "" ///
        _rdistnn_diagnostics_save, matrix(`report') saving(`saving') `replace'

    if "`scope'" == "targets" {
        return matrix localization = `localization'
    }
    return matrix diagnostics = `out'
    return matrix compact = `compact'
    return scalar n_rows = `m'
    return local selected "`selected'"
    return local scope "`scope'"
    return local saved `saving'
end
