*! version 0.1.0  2026-10-17  Jakob Juergens
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
        di as error "xcols() is required."
        exit 198
    }
    if missing(`clip') | `clip'<=0 | `clip'>=.5 | 1-`clip'>=1 {
        di as error "clip() must be finite and strictly between zero and 0.5, with 1-clip() < 1."
        exit 198
    }
    if missing(`denominatortol') | `denominatortol'<0 {
        di as error "denominatortol() must be finite and nonnegative."
        exit 198
    }
    if "`training'"=="" local training "supplied"
    if !inlist("`training'","supplied","oracle","full_sample","oob","cross_fit") {
        di as error "training() must be supplied, oracle, full_sample, oob, or cross_fit."
        exit 198
    }
    local route "augmented_iv"
    local nuisance_vars `dcol' `acol' `ecol' `my0col' `my1col' `md0col' `md1col'
    if "`pseudooutcome'"!="" {
        local route "paired_signals"
        if "`denominator'"=="" | "`nuisance_vars'`onesided'"!="" {
            di as error "pseudooutcome requires denominator() and cannot be combined with nuisance options or onesided."
            exit 198
        }
    }
    else {
        if "`denominator'"!="" {
            di as error "denominator() requires pseudooutcome."
            exit 198
        }
        foreach arg in dcol acol ecol my0col my1col md1col {
            if "``arg''"=="" {
                di as error "`arg'() is required for augmented IV signals."
                exit 198
            }
        }
        if "`md0col'"=="" & "`onesided'"=="" {
            di as error "md0col() is required unless onesided is specified."
            exit 198
        }
    }

    marksample touse, novarlist
    if "`route'"=="augmented_iv" {
        foreach prob in ecol md0col md1col {
            if "``prob''"!="" {
                quietly count if `touse' & (missing(``prob'') | ``prob''<0 | ``prob''>1)
                if r(N)>0 {
                    di as error "`prob'() must contain finite predictions in [0,1] before clipping."
                    exit 198
                }
            }
        }
        markout `touse' `input' `xvars' `nuisance_vars'
        foreach binary in dcol acol {
            quietly count if `touse' & !inlist(``binary'',0,1)
            if r(N)>0 {
                di as error "`binary'() must be coded 0/1."
                exit 198
            }
        }
        if "`onesided'"!="" {
            quietly count if `touse' & `dcol'>`acol'
            if r(N)>0 {
                di as error "onesided requires D <= A on the estimation sample."
                exit 198
            }
            if "`md0col'"!="" {
                quietly count if `touse' & `md0col'!=0
                if r(N)>0 {
                    di as error "md0col() must be zero when onesided is specified."
                    exit 198
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
        di as error "Inputs produce non-finite augmented signals; rescale outcome inputs."
        exit 498
    }

    local kx : word count `xvars'
    if ("`x0'" == "") == ("`x0matrix'" == "") {
        di as error "Specify exactly one of x0() or x0matrix()."
        exit 198
    }

    tempname targets_mat neff_mat
    local target_source "x0"
    if "`x0matrix'" != "" {
        capture confirm matrix `x0matrix'
        if _rc {
            di as error "x0matrix() must name an existing numeric Stata matrix."
            exit 198
        }
        matrix `targets_mat' = `x0matrix'
        local target_source "X0"
    }
    else {
        local nx0 : word count `x0'
        if `nx0' != `kx' {
            di as error "x0() must have `kx' elements (one per covariate in xcols())."
            exit 198
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
        di as error "x0matrix() must contain at least one target row."
        exit 198
    }
    if `ntarget_cols' != `kx' {
        di as error "Target matrix must have `kx' columns (one per covariate in xcols()); got `ntarget_cols'."
        exit 198
    }
    mata: st_numscalar("r(target_missing)", hasmissing(st_matrix("`targets_mat'")))
    if scalar(r(target_missing)) {
        di as error "Evaluation targets may not contain missing values."
        exit 198
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
        di as error "k() must be >= 1."
        exit 198
    }
    if "`distance'"=="" local distance "l2"
    if !inlist("`distance'","l2","mahalanobis","l1","l1_scaled","l2_scaled") {
        di as error "distance() must be l1, l2, l1_scaled, l2_scaled, or mahalanobis."
        exit 198
    }
    if "`method'"=="" local method "exact"
    if !inlist("`method'","exact","asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }
    local tuning_mode=lower("`tuningmode'")
    if "`tuning_mode'"=="" local tuning_mode "default"
    if !inlist("`tuning_mode'","default","grid") {
        di as error "tuningmode() must be default or grid."
        exit 198
    }
    if "`distance'"=="mahalanobis" & `n'<2 {
        di as error "distance(mahalanobis) requires at least two complete estimation observations."
        exit 2001
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
    mata {
        _cl_X = st_data(., "`xvars'", "`touse'")
        rdistnn_duplicate_warning(_cl_X)
        _cl_X0 = st_matrix("`targets'")
        _cl_geometry = _cl_DX = _cl_ranks = J(0,0,.)
        rdistnn_localization_geometry(_cl_X,_cl_X0,"`resolved_distance'","`scalerule'",
            strtoreal(tokens("`scalevec'")),`scaleeps',_cl_geometry,_cl_DX,_cl_ranks)
        _cl_signals = st_data(., "`phiy' `phid'", "`touse'")
        _cl_meta = st_matrix("`meta'")
        _cl_grid = _cl_gridse = _cl_gridy = _cl_gridd = J(`ng',`m',.)
        _cl_gridstatus = J(`ng'*`m',7,.)
        _cl_c = _cl_V = _cl_CV = _cl_status = J(0,0,.)
        _cl_b = J(1,0,.)
        for (_cl_g=1; _cl_g<=`ng'; _cl_g++) {
            if (_cl_meta[_cl_g,4]==0) continue
            _cl_prior=0
            for (_cl_h=1; _cl_h<_cl_g; _cl_h++) {
                if (_cl_meta[_cl_h,4] & _cl_meta[_cl_h,2]==_cl_meta[_cl_g,2]) {
                    _cl_prior=_cl_h
                    break
                }
            }
            if (_cl_prior) {
                _cl_grid[_cl_g,.]=_cl_grid[_cl_prior,.]
                _cl_gridse[_cl_g,.]=_cl_gridse[_cl_prior,.]
                _cl_gridy[_cl_g,.]=_cl_gridy[_cl_prior,.]
                _cl_gridd[_cl_g,.]=_cl_gridd[_cl_prior,.]
                _cl_gridstatus[((_cl_g-1)*`m'+1)..(_cl_g*`m'),.] =
                    (J(`m',1,_cl_g),_cl_gridstatus[((_cl_prior-1)*`m'+1)..(_cl_prior*`m'),2..7])
                continue
            }
            rdistnn_clate_fit(_cl_ranks,_cl_signals,`k',_cl_meta[_cl_g,8],
                "`resolved_method'",`active',`denominatortol',
                _cl_c,_cl_b,_cl_V,_cl_CV,_cl_status,_cl_X)
            _cl_grid[_cl_g,.] = _cl_b
            _cl_gridse[_cl_g,.] = sqrt(diagonal(_cl_V))'
            _cl_gridy[_cl_g,.] = _cl_c[,1]'
            _cl_gridd[_cl_g,.] = _cl_c[,2]'
            _cl_gridstatus[((_cl_g-1)*`m'+1)..(_cl_g*`m'),.] =
                (J(`m',1,_cl_g),(1::`m'),_cl_status)
            if (_cl_meta[_cl_g,2]==_cl_meta[`default',2]) {
                st_matrix("`components'",_cl_c)
                st_matrix("`theta'",_cl_b)
                st_matrix("`vcov'",_cl_V)
                st_matrix("`componentV'",_cl_CV)
                st_matrix("`status'",_cl_status)
            }
        }
        _cl_b = st_matrix("`theta'")
        _cl_W = rdistnn_weights_from_rmat(_cl_ranks,`k',`resolved_s',"`resolved_method'",_cl_X)
        _cl_A = _cl_e = J(0,1,.)
        if ("`route'"=="augmented_iv") {
            _cl_A = st_data(.,"`acol'","`touse'")
            _cl_e = st_data(.,"`eclip'","`touse'")
        }
        st_matrix("`diag'",rdistnn_diagnostics("clate_dnn",_cl_X,_cl_X0,
            _cl_signals,_cl_W,_cl_DX,`k',`resolved_s',_cl_b,_cl_A,_cl_e))
        st_matrix("`neff'",1:/colsum(_cl_W:^2))
        _cl_grid_diag=J(`ng',5,.)
        for(_cl_g=1;_cl_g<=`ng';_cl_g++) {
            _cl_grid_diag[_cl_g,3]=min((5,`n'))
            if(_cl_meta[_cl_g,4]==0) continue
            _cl_grid_W=rdistnn_weights_from_rmat(_cl_ranks,`k',_cl_meta[_cl_g,8],"`resolved_method'",_cl_X)
            _cl_grid_diag[_cl_g,.]=rdistnn_grid_weight_diagnostics(_cl_grid_W)
        }
        st_matrix("`grid_diag'",_cl_grid_diag)
        st_matrix("`semat'",sqrt(diagonal(st_matrix("`vcov'")))')
        st_matrix("`gridtheta'",_cl_grid)
        st_matrix("`gridse'",_cl_gridse)
        st_matrix("`gridy'",_cl_gridy)
        st_matrix("`gridd'",_cl_gridd)
        st_matrix("`gridstatus'",_cl_gridstatus)
    }

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
        mata: _cl_keep=strtoreal(tokens("`keep'")); st_matrix("`postb'",st_matrix("`theta'")[,_cl_keep]); st_matrix("`postV'",st_matrix("`vcov'")[_cl_keep,_cl_keep]); st_local("postcov",strofreal(!hasmissing(st_matrix("`postV'"))))
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
    tempname model_signals model_ranks model_distances model_geometry
    mata: st_matrix("`model_signals'",_cl_signals); st_matrix("`model_ranks'",_cl_ranks); st_matrix("`model_distances'",_cl_DX)
    local geometry_option ""
    if inlist("`e(distance)'","mahalanobis","l1_scaled","l2_scaled") {
        mata: st_matrix("`model_geometry'",_cl_geometry)
        local geometry_option "geometry(`model_geometry')"
    }
    quietly _rdistnn_model_store `signature_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances') `geometry_option'
    _clate_dnn_Display, level(`level')
end

program define _clate_dnn_Display, rclass
    version 17
    _rdistnn_report `0'
    return add
end
