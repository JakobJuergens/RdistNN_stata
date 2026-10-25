*! version 0.1.5  2026-10-25  Jakob Juergens
*! Conditional LATE with supplied binary-IV nuisance predictions

program define clate_dnn, eclass
    version 17
    if replay() {
        if "`e(cmd)'" != "clate_dnn" error 301
        _clate_dnn_Display `0'
        exit
    }
    syntax varname(numeric) [if] [in], [ ///
        XCOLS(varlist numeric) ///
        X0(numlist) X0MATRIX(name) ///
        DCOL(varname numeric) ACOL(varname numeric) ECOL(varname numeric) ///
        MY0COL(varname numeric) MY1COL(varname numeric) ///
        MD0COL(varname numeric) MD1COL(varname numeric) ///
        PSEUDOOutcome DENOMinator(varname numeric) ONESIDED ///
        CLIP(real 1e-8) DENOMINATORTOL(real 1.490116119384766e-8) ///
        TRAINING(string) K(integer 1) S(integer 0) ///
        TUNINGMode(string) LAMBDA(real .4) LAMBDAGrid(numlist) ///
        SGRID(numlist) SDEFAULT(integer 0) ///
        DISTance(string) SCALERule(string) SCALEVec(numlist) SCALEEPS(real 1e-12) METHod(string) Level(cilevel) noSE]

    local original_cmdline `"clate_dnn `0'"'
    local input `varlist'
    local xvars `xcols'
    if "`xvars'"=="" {
        _rdistnn_messages error, rc(198) message(`"xcols() is required."')
    }
    if missing(`clip') | `clip'<=0 | `clip'>=.5 | 1-`clip'>=1 {
        _rdistnn_messages error, rc(198) message(`"clip() must be finite and strictly between zero and 0.5, with 1-clip() < 1."')
    }
    if missing(`denominatortol') | `denominatortol'<0 {
        _rdistnn_messages error, rc(198) message(`"denominatortol() must be finite and nonnegative."')
    }
    if "`training'"=="" local training "supplied"
    if !inlist("`training'","supplied","oracle","full_sample","oob","cross_fit") {
        _rdistnn_messages error, rc(198) message(`"training() must be supplied, oracle, full_sample, oob, or cross_fit."')
    }
    local route "augmented_iv"
    local nuisance_vars `dcol' `acol' `ecol' `my0col' `my1col' `md0col' `md1col'
    if "`pseudooutcome'"!="" {
        local route "paired_signals"
        if "`denominator'"=="" | "`nuisance_vars'`onesided'"!="" {
            _rdistnn_messages error, rc(198) message(`"pseudooutcome requires denominator() and cannot be combined with nuisance options or onesided."')
        }
    }
    else {
        if "`denominator'"!="" {
            _rdistnn_messages error, rc(198) message(`"denominator() requires pseudooutcome."')
        }
        foreach arg in dcol acol ecol my0col my1col md1col {
            if "``arg''"=="" {
                _rdistnn_messages error, rc(198) message(`"`arg'() is required for augmented IV signals."')
            }
        }
        if "`md0col'"=="" & "`onesided'"=="" {
            _rdistnn_messages error, rc(198) message(`"md0col() is required unless onesided is specified."')
        }
    }

    marksample touse, novarlist
    if "`route'"=="augmented_iv" {
        foreach prob in ecol md0col md1col {
            if "``prob''"!="" {
                quietly count if `touse' & (missing(``prob'') | ``prob''<0 | ``prob''>1)
                if r(N)>0 {
                    _rdistnn_messages error, rc(198) message(`"`prob'() must contain finite predictions in [0,1] before clipping."')
                }
            }
        }
        markout `touse' `input' `xvars' `nuisance_vars'
        foreach binary in dcol acol {
            quietly count if `touse' & !inlist(``binary'',0,1)
            if r(N)>0 {
                _rdistnn_messages error, rc(198) message(`"`binary'() must be coded 0/1."')
            }
        }
        if "`onesided'"!="" {
            quietly count if `touse' & `dcol'>`acol'
            if r(N)>0 {
                _rdistnn_messages error, rc(198) message(`"onesided requires D <= A on the estimation sample."')
            }
            if "`md0col'"!="" {
                quietly count if `touse' & `md0col'!=0
                if r(N)>0 {
                    _rdistnn_messages error, rc(198) message(`"md0col() must be zero when onesided is specified."')
                }
            }
        }
    }
    else markout `touse' `input' `denominator' `xvars'
    quietly count if `touse'
    if r(N)==0 error 2000

    tempvar phiy phid eclip
    if "`route'"=="augmented_iv" {
        local md0value `md0col'
        if "`md0value'"=="" local md0value "0"
        quietly gen double `eclip' = min(max(`ecol',`clip'),1-`clip') if `touse'
        quietly gen double `phiy' = `my1col'-`my0col' + ///
            `acol'*(`input'-`my1col')/`eclip' - ///
            (1-`acol')*(`input'-`my0col')/(1-`eclip') if `touse'
        quietly gen double `phid' = `md1col'-`md0value' + ///
            `acol'*(`dcol'-`md1col')/`eclip' - ///
            (1-`acol')*(`dcol'-`md0value')/(1-`eclip') if `touse'
    }
    else {
        quietly gen double `phiy' = `input' if `touse'
        quietly gen double `phid' = `denominator' if `touse'
    }
    quietly count if `touse' & (missing(`phiy') | missing(`phid'))
    if r(N)>0 {
        _rdistnn_messages error, rc(498) message(`"Inputs produce non-finite augmented signals; rescale outcome inputs."')
    }

    local kx : word count `xvars'
    if ("`x0'" == "") == ("`x0matrix'" == "") {
        _rdistnn_messages error, rc(198) message(`"Specify exactly one of x0() or x0matrix()."')
    }

    tempname targets_mat neff_mat
    local target_source "x0"
    if "`x0matrix'" != "" {
        confirm matrix `x0matrix'
        matrix `targets_mat' = `x0matrix'
        local target_source "X0"
    }
    else {
        local nx0 : word count `x0'
        if `nx0' != `kx' {
            _rdistnn_messages error, rc(198) message("x0() requires `kx' values; got `nx0'.") hint("Supply one value per covariate in xcols().")
        }
        local x0_matrix ""
        foreach v of local x0 {
            local x0_matrix `x0_matrix' `v',
        }
        local x0_matrix = substr("`x0_matrix'", 1, length("`x0_matrix'") - 1)
        matrix `targets_mat' = (`x0_matrix')
    }

    local ntargets = rowsof(`targets_mat')
    local ntarget_cols = colsof(`targets_mat')
    if `ntargets' < 1 {
        _rdistnn_messages error, rc(198) message(`"x0matrix() must contain at least one target row."')
    }
    if `ntarget_cols' != `kx' {
        _rdistnn_messages error, rc(503) message(`"Target matrix must have `kx' columns (one per covariate in xcols()); got `ntarget_cols'."')
    }
    mata: st_numscalar("r(target_missing)", hasmissing(st_matrix("`targets_mat'")))
    if scalar(r(target_missing)) {
        _rdistnn_messages error, rc(198) message(`"Evaluation targets may not contain missing values."')
    }

    local target_names ""
    forvalues j = 1/`ntargets' {
        local target_names `target_names' target`j'
    }
    matrix rownames `targets_mat' = `target_names'
    matrix colnames `targets_mat' = `xvars'

    quietly count if `touse'
    local n=r(N)
    if `k'<1 {
        _rdistnn_messages error, rc(198) message("k() must be a positive integer; got `k'.") hint("Set k() to at least 1.")
    }
    if "`distance'"=="" local distance "l2"
    if !inlist("`distance'","l2","mahalanobis","l1","l1_scaled","l2_scaled") {
        _rdistnn_messages error, rc(198) message(`"distance() must be l1, l2, l1_scaled, l2_scaled, or mahalanobis."')
    }
    if "`method'"=="" local method "exact"
    if !inlist("`method'","exact","asymp") {
        _rdistnn_messages error, rc(198) message(`"method() must be exact or asymp."')
    }
    local tuning_mode=lower("`tuningmode'")
    if "`tuning_mode'"=="" local tuning_mode "default"
    if !inlist("`tuning_mode'","default","grid") {
        _rdistnn_messages error, rc(198) message(`"tuningmode() must be default or grid."')
    }
    if "`distance'"=="mahalanobis" & `n'<2 {
        _rdistnn_messages error, rc(2001) message(`"distance(mahalanobis) requires at least two complete estimation observations."')
    }
    tempname candidate_mat grid_est grid_se grid_diag target_diag
    quietly _rdistnn_tuning_plan, n(`n') kx(`kx') k(`k') s(`s') lambda(`lambda') ///
        mode(`tuning_mode') sdefault(`sdefault') sgrid(`sgrid') lambdagrid(`lambdagrid') `se'
    matrix `candidate_mat'=r(candidates)
    local s=r(s)
    local default_lambda=r(default_lambda)
    local default_index=r(default_index)
    local ncandidates=r(ncandidates)
    local nfeasible=r(nfeasible)
    local nunique_feasible=r(nunique_feasible)
    local tuning_source "`r(tuning_source)'"
    local candidate_names "`r(candidate_names)'"
    local m=`ntargets'
    local resolved_s=`s'
    local resolved_method "`method'"
    local resolved_distance "`distance'"
    local isgrid=("`tuning_mode'"=="grid")
    local active=("`se'"=="")
    local ng=`ncandidates'
    local default=`default_index'
    tempname targets diag components theta vcov componentV status semat ///
        meta gridtheta gridse gridy gridd gridstatus neff
    matrix `targets'=`targets_mat'
    matrix `meta'=`candidate_mat'
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    findfile "_rdistnn_clate_mata.ado"
    quietly run "`r(fn)'"
    tempname model_groups model_signals model_ranks model_distances model_geometry
    mata: rdistnn_clate_fit_ado()

    // Post only the paired estimator's coefficients and uncertainty.
    local cnames ""
    local validnames ""
    local keep ""
    local compnames ""
    forvalues j=1/`m' {
        local coef "theta`j'"
        if `m'==1 local coef "theta"
        local cnames `cnames' `coef'
        local compnames `compnames' y`j' d`j'
        if el(`theta',1,`j')<. {
            local keep `keep' `j'
            local validnames `validnames' `coef'
        }
    }
    if "`keep'"!="" {
        tempname postb postV
        mata: rdistnn_clate_post_ado()
        matrix colnames `postb' = `validnames'
        matrix colnames `postV' = `validnames'
        matrix rownames `postV' = `validnames'
        if `active' & `postcov' ereturn post `postb' `postV', esample(`touse')
        else ereturn post `postb', esample(`touse')
    }
    else ereturn post, esample(`touse')
    ereturn scalar N=`n'
    ereturn scalar m=`m'
    ereturn scalar p=`kx'
    ereturn scalar kx=`kx'
    ereturn scalar k=`k'
    ereturn scalar s=`resolved_s'
    ereturn scalar lambda=`default_lambda'
    mata: st_numscalar("r(_cl_min)",min(st_matrix("`neff'"))); st_numscalar("r(_cl_max)",max(st_matrix("`neff'")))
    // Retain the scalar-copy boundary of the former localization adapter.
    // The canonical matrix values keep their full numeric precision.
    local localization_min=r(_cl_min)
    local localization_max=r(_cl_max)
    ereturn scalar neff_min=`localization_min'
    ereturn scalar neff_max=`localization_max'
    ereturn scalar ess_min=`localization_min'
    ereturn scalar ess_max=`localization_max'
    if `m'==1 {
        local localization_neff=el(`neff',1,1)
        ereturn scalar neff=`localization_neff'
        ereturn scalar ess=`localization_neff'
    }
    local neffnames ""
    forvalues j=1/`m' {
        local neffname "tau"
        if `m'>1 local neffname "tau`j'"
        local neffnames "`neffnames' `neffname'"
    }
    matrix colnames `neff'=`neffnames'
    ereturn matrix targets=`targets', copy
    ereturn matrix x0_matrix=`targets', copy
    ereturn matrix neff_targets=`neff', copy
    ereturn matrix ess_targets=`neff', copy
    if `isgrid' {
        matrix rownames `grid_diag'=`candidate_names'
        matrix colnames `grid_diag'=ess max_weight top_n top5_weight_share degenerate
        ereturn matrix tuning_candidates=`meta', copy
        ereturn matrix tuning_diagnostics=`grid_diag', copy
        ereturn scalar tuning_default_candidate=`default_index'
        ereturn scalar tuning_n_candidates=`ncandidates'
        ereturn scalar tuning_n_feasible=`nfeasible'
        ereturn scalar tuning_n_skipped=`ncandidates'-`nfeasible'
        ereturn scalar tuning_n_failed=0
        ereturn scalar tuning_n_unique=`nunique_feasible'
    }
    if "`scalerule'"=="" local scalerule "sd"
    ereturn local scalerule "`scalerule'"
    ereturn scalar scaleeps=`scaleeps'
    ereturn local x0 "`x0'"
    ereturn local x0matrix "`x0matrix'"
    ereturn local target_source "`target_source'"
    ereturn local xvars "`xvars'"
    ereturn local x_cols "`xvars'"
    ereturn local distance "`resolved_distance'"
    ereturn local method "`resolved_method'"
    ereturn local tuning_mode "`tuning_mode'"
    ereturn local tuning_source "`tuning_source'"
    ereturn local tuning_status_codes "0=infeasible 1=ok"
    ereturn local tuning_reason_codes "0=ok 1=k_gt_s 2=jackknife_requires_s_lt_n 3=s_gt_n"
    ereturn local jackknife=cond(`active',"fast","none")
    matrix colnames `theta' = `cnames'
    matrix colnames `semat' = `cnames'
    matrix colnames `vcov' = `cnames'
    matrix rownames `vcov' = `cnames'
    matrix colnames `componentV' = `compnames'
    matrix rownames `componentV' = `compnames'
    matrix rownames `components' = `: rownames `targets''
    matrix colnames `components' = reduced_form first_stage
    matrix rownames `status' = `: rownames `targets''
    matrix colnames `status' = ratio_status negative_denominator ///
        near_zero_denominator covariance_status wald_available
    matrix rownames `diag' = `: rownames `targets''
    matrix colnames `diag'=ess max_weight top1_weight_share top5_weight_share top10_weight_share ///
        n_nonzero_weight distance_available radius_rank_k radius_rank_s weighted_mean_distance ///
        weighted_median_distance overlap_available treated_mass control_mass treated_share ///
        propensity_available weighted_mean_propensity weighted_propensity_q10 weighted_propensity_q90 ///
        physical_radius_90 arm_ess_available treated_ess control_ess score_contribution_available top5_score_share
    ereturn matrix estimates = `theta', copy
    ereturn matrix se_targets = `semat', copy
    ereturn matrix theta_vcov = `vcov', copy
    ereturn matrix component_vcov = `componentV', copy
    ereturn matrix components = `components', copy
    ereturn matrix status = `status', copy
    ereturn matrix diagnostics = `diag', copy
    ereturn matrix target_diagnostics = `diag', copy
    if `isgrid' {
        foreach mat in gridtheta gridse gridy gridd {
            matrix colnames ``mat'' = `cnames'
            matrix rownames ``mat'' = `: rownames `meta''
        }
        matrix colnames `gridstatus' = candidate_id target_id ratio_status ///
            negative_denominator near_zero_denominator covariance_status wald_available
        ereturn matrix tuning_estimates = `gridtheta', copy
        ereturn matrix tuning_se = `gridse', copy
        ereturn matrix tuning_reduced_form = `gridy', copy
        ereturn matrix tuning_first_stage = `gridd', copy
        ereturn matrix tuning_target_status = `gridstatus', copy
    }
    ereturn scalar clip = `clip'
    ereturn scalar denominator_tol = `denominatortol'
    ereturn scalar one_sided = ("`onesided'"!="")
    ereturn local cmd "clate_dnn"
    ereturn local cmdline `"`original_cmdline'"'
    ereturn local score_method "`route'"
    ereturn local se_method = cond(`active',"paired_rank_shift","none")
    ereturn local vce = cond(`active',"paired_rank_shift","none")
    ereturn local vcetype = cond(`active',"Paired frozen-score jackknife","")
    ereturn local depvar "`input'"
    ereturn local covariance = cond(`active',"paired_frozen_score","none")
    ereturn local outcome "`input'"
    ereturn local signal "`input'"
    ereturn local denominator "`denominator'"
    ereturn local d_col "`dcol'"
    ereturn local a_col "`acol'"
    ereturn local e_col "`ecol'"
    ereturn local propensity_version = cond("`route'"=="augmented_iv","clipped","not_used")
    ereturn local my0_col "`my0col'"
    ereturn local my1_col "`my1col'"
    ereturn local md0_col "`md0col'"
    ereturn local md1_col "`md1col'"
    ereturn local training "`training'"
    ereturn local overlap_role "instrument"
    ereturn local ratio_status_codes "0=ok 1=zero_denominator 2=numerical_failure"
    ereturn local covariance_status_codes "0=ok 1=not_computed 2=zero_score_variance 3=numerical_failure"
    local signature_vars `xvars' `input' `nuisance_vars' `denominator'
    local signature_vars : list uniq signature_vars
    quietly _rdistnn_estimation_state `signature_vars'
    local geometry_option ""
    if inlist("`e(distance)'","mahalanobis","l1_scaled","l2_scaled") {
        local geometry_option "geometry(`model_geometry')"
    }
    quietly _rdistnn_model_store `signature_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances') groups(`model_groups') `geometry_option'
    _clate_dnn_Display, level(`level')
end

program define _clate_dnn_Display, rclass
    version 17
    syntax [, Level(cilevel)]
    _rdistnn_report, level(`level')
    return add
end
