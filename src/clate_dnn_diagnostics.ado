*! version 0.1.5  2026-10-25  Jakob Juergens
*! Postestimation diagnostics tables for clate_dnn

program define clate_dnn_diagnostics, rclass
    version 17
    syntax [, TARGETS(numlist integer >=1) CANDIDATES(numlist integer >=1) ///
        ALLCandidates DETAIL SAVING(string asis) REPLACE]

    local targetopt ""
    local candidateopt ""
    if "`targets'" != "" local targetopt "targets(`targets')"
    if "`candidates'" != "" local candidateopt "candidates(`candidates')"
    _rdistnn_diagnostics, command(clate_dnn) `targetopt' `candidateopt' `allcandidates'
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
            weighted_median_distance overlap_available instrument_1_mass ///
            instrument_0_mass instrument_1_share propensity_available ///
            weighted_mean_propensity weighted_propensity_q10 ///
            weighted_propensity_q90 physical_radius_90 ///
            arm_ess_available instrument_1_ess instrument_0_ess ///
            score_contribution_available top5_score_share
        mata: st_matrix("`compact'", st_matrix("`out'")[,(1, `p'+2, `p'+5, `p'+21, `p'+23, `p'+24, `p'+26)])
        matrix colnames `compact' = target_id ess top5_weight_share ///
            physical_radius_90 instrument_1_ess instrument_0_ess top5_score_share
        matrix rownames `compact' = `: rownames `out''
        tempname components status component_raw component_vcov estimates status_raw
        matrix `component_raw' = e(components)
        matrix `component_vcov' = e(component_vcov)
        matrix `estimates' = e(estimates)
        matrix `status_raw' = e(status)
        matrix `components' = J(`m',9,.)
        local i = 0
        foreach j of local selected {
            local ++i
            matrix `components'[`i',1] = `j'
            matrix `components'[`i',2] = el(`component_raw',`j',1)
            matrix `components'[`i',3] = el(`component_raw',`j',2)
            matrix `components'[`i',4] = sqrt(el(`component_vcov',2*`j'-1,2*`j'-1))
            matrix `components'[`i',5] = sqrt(el(`component_vcov',2*`j',2*`j'))
            matrix `components'[`i',6] = el(`component_vcov',2*`j'-1,2*`j')
            if el(`components',`i',4)>0 & el(`components',`i',4)<. & ///
                el(`components',`i',5)>0 & el(`components',`i',5)<. ///
                matrix `components'[`i',7] = el(`components',`i',6) / ///
                    (el(`components',`i',4)*el(`components',`i',5))
            if el(`components',`i',5)>0 & el(`components',`i',5)<. ///
                matrix `components'[`i',8] = el(`components',`i',3)/el(`components',`i',5)
            matrix `components'[`i',9] = el(`estimates',1,`j')
        }
        matrix colnames `components' = target_id reduced_form first_stage ///
            numerator_se denominator_se component_covariance component_correlation ///
            denominator_stat estimate
        matrix rownames `components' = `: rownames `out''
        mata: st_matrix("`status'", st_matrix("`status_raw'")[strtoreal(tokens("`selected'"))',])
        matrix colnames `status' = `: colnames `status_raw''
        matrix rownames `status' = `: rownames `out''
    }
    else matrix `compact' = r(compact)

    if "`detail'" != "" matrix list `out', format(%12.6g)
    else matrix list `compact', format(%12.6g)
    if "`scope'" == "targets" {
        matrix list `components', format(%12.6g)
        matrix list `status', format(%9.0g)
        di as text "Arms refer to the instrument A; these are pointwise Wald diagnostics, without weak-IV-uniform guarantees."
    }

    if "`scope'" == "targets" matrix `report' = `out', `components'[1...,2...], `status'
    else matrix `report' = `out'
    if `"`saving'"' != "" ///
        _rdistnn_diagnostics_save, matrix(`report') saving(`saving') `replace'

    if "`scope'" == "targets" {
        return matrix localization = `localization'
        return matrix components = `components'
        return matrix status = `status'
    }
    return matrix diagnostics = `out'
    return matrix compact = `compact'
    return scalar n_rows = `m'
    return local selected "`selected'"
    return local scope "`scope'"
    return local saved `saving'
end
