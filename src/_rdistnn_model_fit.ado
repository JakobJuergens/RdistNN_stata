*! version 0.1.0  2026-10-18
// Post a derived fit from frozen numeric state; callers control e() restoration.
program define _rdistnn_model_fit, eclass
    version 17
    syntax, TARGETS(name) S(integer) [COVariance(string) CANDIDATE(integer 0)]
    if !inlist("`e(cmd)'","np_dnn","np_tdnn","cate_dnn","clate_dnn","plqr_dnn") error 301
    if e(model_version)!=1 {
        di as error "Active estimates have no supported saved model; refit with the current package version."
        exit 498
    }
    local matrices : e(matrices)
    foreach field in model_X model_signals model_inputs model_ranks model_distances model_targets model_row_ids {
        if !`: list field in matrices' {
            di as error "Saved model is missing e(`field'); refit with the current package version."
            exit 498
        }
    }
    if "`covariance'"=="" local covariance "compute"
    if !inlist("`covariance'","compute","none") error 198
    local active=("`covariance'"=="compute")
    local cmd "`e(cmd)'"
    local lambda=e(lambda)
    local parent_cmdline `"`e(cmdline)'"'
    if `candidate'>0 local lambda=el(e(tuning_candidates),`candidate',7)
    local p=colsof(`targets')
    local m=rowsof(`targets')
    if `p'!=colsof(e(model_X)) {
        di as error "Targets must have one column per fitted localization covariate."
        exit 198
    }
    foreach file in _rdistnn_mata.ado _rdistnn_clate_mata.ado _rdistnn_model_mata.ado {
        findfile "`file'"
        quietly run "`r(fn)'"
    }
    tempname b V diagnostics components componentV status neff ranks DX W fit ses
    mata: rdistnn_model_evaluate("`targets'",`s',`active',tokens("`b' `V' `diagnostics' `components' `componentV' `status' `neff' `ranks' `DX' `W' `fit'"))
    mata: st_matrix("`ses'",sqrt(diagonal(st_matrix("`V'")))')
    local scalars : e(scalars)
    local macros : e(macros)
    local ns 0
    foreach field of local scalars {
        local ++ns
        local scalar`ns'=e(`field')
    }
    local nm 0
    foreach field of local macros {
        local ++nm
        local macro`nm' `"`e(`field')'"'
    }
    local nt 0
    foreach field of local matrices {
        local ++nt
        tempname saved
        matrix `saved'=e(`field')
        local matrix`nt' "`saved'"
    }
    tempvar sample
    quietly gen byte `sample'=0
    capture quietly _rdistnn_sample_signature `e(signature_vars)' if e(sample)
    if !_rc {
        if `"`r(signature)'"'==`"`e(ordered_signature)'"' quietly replace `sample'=e(sample)
    }
    local labels ""
    local eqs ""
    local q=1
    if "`cmd'"=="plqr_dnn" local q=e(kd)
    forvalues j=1/`m' {
        if `q'==1 & "`cmd'"!="plqr_dnn" local labels "`labels' target`j'"
        else {
            local labels "`labels' `e(dvars)'"
            forvalues c=1/`q' {
                local eqs "`eqs' target`j'"
            }
        }
    }
    matrix colnames `b'=`labels'
    matrix colnames `ses'=`labels'
    matrix colnames `V'=`labels'
    matrix rownames `V'=`labels'
    if "`eqs'"!="" {
        matrix coleq `b'=`eqs'
        matrix coleq `ses'=`eqs'
        matrix coleq `V'=`eqs'
        matrix roweq `V'=`eqs'
    }
    tempname postb postV
    if "`cmd'"=="clate_dnn" {
        mata: st_matrix("`postb'",select(st_matrix("`b'"),st_matrix("`b'"):<.))
        // Undefined targets remain in canonical matrices, outside e(b).
        local valid ""
        forvalues j=1/`m' {
            if el(`b',1,`j')<. local valid "`valid' `j'"
        }
        if "`valid'"!="" {
            mata: st_matrix("`postV'",st_matrix("`V'")[strtoreal(tokens("`valid'")),strtoreal(tokens("`valid'"))])
            local posted_labels ""
            foreach j of local valid {
                local posted_labels "`posted_labels' target`j'"
            }
            matrix colnames `postb'=`posted_labels'
            matrix colnames `postV'=`posted_labels'
            matrix rownames `postV'=`posted_labels'
            mata: st_numscalar("r(_model_postcov)",!hasmissing(st_matrix("`postV'")))
            if `active' & r(_model_postcov) ereturn post `postb' `postV', esample(`sample')
            else ereturn post `postb', esample(`sample')
        }
        else ereturn post, esample(`sample')
    }
    else {
        matrix `postb'=`b'
        matrix `postV'=`V'
        if `active' ereturn post `postb' `postV', esample(`sample')
        else ereturn post `postb', esample(`sample')
    }
    local i 0
    foreach field of local scalars {
        local ++i
        if substr("`field'",1,7)!="tuning_" ereturn scalar `field'=`scalar`i''
    }
    local i 0
    foreach field of local macros {
        local ++i
        if substr("`field'",1,7)!="tuning_" & "`field'"!="properties" ereturn local `field' `"`macro`i''"'
    }
    local i 0
    foreach field of local matrices {
        local ++i
        if substr("`field'",1,7)!="tuning_" & !inlist("`field'","b","V") ereturn matrix `field'=`matrix`i'', copy
    }
    ereturn scalar m=`m'
    ereturn scalar p=`p'
    ereturn scalar kx=`p'
    ereturn scalar s=`s'
    ereturn scalar lambda=`lambda'
    if `candidate'>0 ereturn scalar candidate_id=`candidate'
    mata: st_numscalar("r(_model_min)",min(st_matrix("`neff'"))); st_numscalar("r(_model_max)",max(st_matrix("`neff'")))
    ereturn scalar neff_min=r(_model_min)
    ereturn scalar neff_max=r(_model_max)
    ereturn scalar neff=cond(`m'==1,el(`neff',1,1),.)
    ereturn scalar ess=cond(`m'==1,el(`neff',1,1),.)
    ereturn scalar ess_min=e(neff_min)
    ereturn scalar ess_max=e(neff_max)
    ereturn matrix estimates=`b', copy
    ereturn matrix se_targets=`ses', copy
    ereturn matrix estimate_vcov=`V', copy
    ereturn matrix targets=`targets', copy
    ereturn matrix x0_matrix=`targets', copy
    ereturn matrix neff_targets=`neff', copy
    ereturn matrix ess_targets=`neff', copy
    local cols "ess max_weight top1_weight_share top5_weight_share top10_weight_share n_nonzero_weight distance_available radius_rank_k radius_rank_s weighted_mean_distance weighted_median_distance"
    if "`cmd'"=="np_tdnn" local cols "signed_sum positive_mass negative_mass abs_mass inverse_l2_concentration max_abs_weight n_nonzero_weight abs_mean_distance"
    if "`cmd'"=="np_dnn" local cols "`cols' outcome_available weighted_mean_y weighted_sd_y weighted_y_q10 weighted_y_q50 weighted_y_q90"
    if inlist("`cmd'","cate_dnn","clate_dnn") local cols "`cols' overlap_available treated_mass control_mass treated_share propensity_available weighted_mean_propensity weighted_propensity_q10 weighted_propensity_q90 physical_radius_90 arm_ess_available treated_ess control_ess score_contribution_available top5_score_share"
    matrix colnames `diagnostics'=`cols'
    ereturn matrix diagnostics=`diagnostics', copy
    ereturn matrix target_diagnostics=`diagnostics', copy
    ereturn matrix model_ranks=`ranks', copy
    ereturn matrix model_distances=`DX', copy
    ereturn matrix model_targets=`targets', copy
    ereturn matrix fitted_values=`fit', copy
    if "`cmd'"=="clate_dnn" {
        matrix colnames `components'=reduced_form first_stage
        matrix colnames `status'=ratio_status negative_denominator near_zero_denominator covariance_status wald_available
        local component_labels ""
        forvalues j=1/`m' {
            local component_labels "`component_labels' y`j' d`j'"
        }
        matrix colnames `componentV'=`component_labels'
        matrix rownames `componentV'=`component_labels'
        ereturn matrix components=`components', copy
        ereturn matrix component_vcov=`componentV', copy
        ereturn matrix theta_vcov=`V', copy
        ereturn matrix status=`status', copy
        ereturn local covariance=cond(`active',"paired_frozen_score","none")
    }
    if "`cmd'"=="plqr_dnn" {
        ereturn matrix moment=`b', copy
        ereturn matrix moment_vcov=`V', copy
        ereturn matrix moment_by_target=`fit', copy
    }
    ereturn local x0 ""
    ereturn local x0matrix ""
    ereturn local target_source "saved_fit"
    ereturn local parent_cmdline `"`parent_cmdline'"'
    if `candidate'>0 ereturn local cmdline "rdistnn_candidate, candidate(`candidate') post"
    else {
        local seopt ""
        if `active' local seopt "sefit"
        ereturn local cmdline "rdistnn_predict, x0matrix(`targets') `seopt' post"
    }
    ereturn local tuning_mode "default"
    ereturn local tuning_source "saved_fit"
    local se_method "none"
    if `active' {
        if "`cmd'"=="clate_dnn" local se_method "paired_rank_shift"
        else if "`cmd'"=="np_tdnn" local se_method "signed_rank_shift"
        else local se_method "rank_shift"
    }
    ereturn local se_method "`se_method'"
    ereturn local jackknife=cond(`active',"fast","none")
end
