*! version 0.1.2  2026-07-12  Jakob Juergens
*! Nonparametric regression via distributional nearest neighbors

program define np_dnn, eclass
    version 17

    if replay() {
        if "`e(cmd)'" != "np_dnn" error 301
        _np_dnn_Display `0'
        exit
    }

    syntax varname(numeric) [if] [in] , ///
        [                                ///
        XCOLS(varlist numeric)           ///  R-aligned x_cols convention
        x0(numlist)                      ///  single evaluation point
        X0MATRIX(name)                   ///  named m x kx target matrix (R X0)
        k(integer 1)                     ///  # selected within each s-subset
        s(integer 0)                     ///  subset/kernel order; larger = tighter
        TUNINGMode(string)               ///  default | grid (R tuning_mode)
        LAMBDA(real .4)                  ///  paper finite-sample multiplier
        LAMBDAGrid(numlist)              ///  paper-rule candidate multipliers
        SGRID(numlist)                   ///  direct Stata candidate orders
        SDEFAULT(integer 0)              ///  default order for sgrid()
        Distance(string)                 ///  localization geometry
        SCALERule(string) SCALEVec(numlist) SCALEEPS(real 1e-12) ///
        Method(string)                   ///  exact (default) | asymp
        Level(cilevel)                   ///  CI level (default 95)
        noSE                             ///  skip jackknife standard errors
        ]

    local outcome `varlist'
    local xvars `xcols'
    if "`xvars'" == "" {
        di as error "xcols() is required."
        exit 198
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

    marksample touse
    markout `touse' `xvars' `outcome'
    quietly count if `touse'
    if r(N) == 0 {
        di as error "No observations."
        exit 2000
    }
    local n = r(N)

    if `k' < 1 {
        _rdistnn_messages error, rc(198) message("k() must be a positive integer; got `k'.") hint("Set k() to at least 1.")
    }
    if "`distance'" == "" local distance "l2"
    if !inlist("`distance'", "l2", "mahalanobis", "l1", "l1_scaled", "l2_scaled") {
        di as error "distance() must be l1, l2, l1_scaled, l2_scaled, or mahalanobis."
        exit 198
    }
    if "`distance'" == "mahalanobis" & `n' < 2 {
        di as error "distance(mahalanobis) requires at least two complete estimation observations."
        exit 2001
    }
    if "`method'" == "" local method "exact"
    if !inlist("`method'", "exact", "asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }

    local tuning_mode = lower("`tuningmode'")
    if "`tuning_mode'" == "" local tuning_mode "default"
    if !inlist("`tuning_mode'", "default", "grid") {
        di as error "tuningmode() must be default or grid."
        exit 198
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

    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    tempname b V
    mata {
        _np_X = st_data(., "`xvars'", "`touse'")
        _np_X0 = st_matrix("`targets_mat'")
        _np_geometry = _np_DX = _np_ranks = J(0,0,.)
        rdistnn_localization_geometry(_np_X,_np_X0,"`distance'","`scalerule'",
            strtoreal(tokens("`scalevec'")),`scaleeps',_np_geometry,_np_DX,_np_ranks)
        _np_groups=rdistnn_distance_groups(_np_ranks,_np_DX,1)
        _np_y = st_data(., "`outcome'", "`touse'")
        _np_candidates = st_matrix("`candidate_mat'")
        _np_grid_est = J(`ncandidates', `ntargets', .)
        _np_grid_se = J(`ncandidates', `ntargets', .)
        _np_grid_diag = J(`ncandidates', 5, .)

        for (_np_g = 1; _np_g <= `ncandidates'; _np_g++) {
            _np_grid_diag[_np_g, 3] = min((5, `n'))
            if (_np_candidates[_np_g, 4] == 0) continue
            _np_prior = 0
            if (_np_g > 1) {
                for (_np_h = 1; _np_h < _np_g; _np_h++) {
                    if (_np_candidates[_np_h, 4] == 1 &
                        _np_candidates[_np_h, 2] == _np_candidates[_np_g, 2]) {
                        _np_prior = _np_h
                        break
                    }
                }
            }
            if (_np_prior > 0) {
                _np_grid_est[_np_g,.] = _np_grid_est[_np_prior,.]
                _np_grid_se[_np_g,.] = _np_grid_se[_np_prior,.]
                _np_grid_diag[_np_g,.] = _np_grid_diag[_np_prior,.]
            }
            else {
                _np_sg = _np_candidates[_np_g, 8]
                _np_Wg = rdistnn_weights_from_rmat(_np_ranks, `k', _np_sg, "`method'",_np_DX,_np_groups)
                _np_grid_est[_np_g,.] = rdistnn_weighted_signals(_np_y, _np_Wg)
                _np_grid_diag[_np_g,.]=rdistnn_grid_weight_diagnostics(_np_Wg)
                if ("`se'" == "") {
                    _np_Vg = rdistnn_jk_cov_multi_from_ranks(
                        _np_ranks, _np_y, `k', _np_sg, "`method'",_np_DX,_np_groups)
                    _np_grid_se[_np_g,.] = sqrt(diagonal(_np_Vg))'
                }
                // Save the default's first representative, including when the
                // default row duplicates an earlier candidate after other orders.
                if (_np_candidates[_np_g, 2] == _np_candidates[`default_index', 2]) {
                    _np_W = _np_Wg
                    if ("`se'" == "") _np_V = _np_Vg
                }
            }
        }

        _np_mu = _np_grid_est[`default_index',.]
        _np_neff = 1 :/ colsum(_np_W :^ 2)
        _np_target_diag = rdistnn_diagnostics("np_dnn",_np_X,_np_X0,
            _np_y,_np_W,_np_DX,`k',`s',_np_mu,J(0,1,.),J(0,1,.))
        st_matrix("`b'", _np_mu)
        st_matrix("`neff_mat'", _np_neff)
        st_matrix("`grid_est'", _np_grid_est)
        st_matrix("`grid_se'", _np_grid_se)
        st_matrix("`grid_diag'", _np_grid_diag)
        st_matrix("`target_diag'", _np_target_diag)
        st_numscalar("r(neff_min)", min(_np_neff))
        st_numscalar("r(neff_max)", max(_np_neff))
    }
    if "`se'" == "" mata: st_matrix("`V'", _np_V)

    local coef_names "mu"
    if `ntargets' > 1 {
        local coef_names ""
        forvalues j = 1/`ntargets' {
            local coef_names `coef_names' mu`j'
        }
    }
    matrix colnames `b' = `coef_names'
    matrix colnames `neff_mat' = `coef_names'
    matrix colnames `grid_est' = `coef_names'
    matrix colnames `grid_se' = `coef_names'
    matrix rownames `grid_est' = `candidate_names'
    matrix rownames `grid_se' = `candidate_names'
    matrix rownames `grid_diag' = `candidate_names'
    matrix colnames `grid_diag' = ess max_weight top_n top5_weight_share degenerate
    matrix rownames `target_diag' = `target_names'
    matrix colnames `target_diag' = ess max_weight top1_weight_share ///
        top5_weight_share top10_weight_share n_nonzero_weight ///
        distance_available radius_rank_k radius_rank_s ///
        weighted_mean_distance weighted_median_distance outcome_available ///
        weighted_mean_y weighted_sd_y weighted_y_q10 weighted_y_q50 ///
        weighted_y_q90
    if "`se'" == "" {
        matrix rownames `V' = `coef_names'
        matrix colnames `V' = `coef_names'
        ereturn post `b' `V', esample(`touse')
    }
    else {
        ereturn post `b', esample(`touse')
    }
    ereturn scalar N = `n'
    ereturn scalar m = `ntargets'
    ereturn scalar neff_min = scalar(r(neff_min))
    ereturn scalar neff_max = scalar(r(neff_max))
    ereturn scalar ess_min = scalar(r(neff_min))
    ereturn scalar ess_max = scalar(r(neff_max))
    if `ntargets' == 1 {
        ereturn scalar neff = el(`neff_mat', 1, 1)
        ereturn scalar ess = el(`neff_mat', 1, 1)
    }
    ereturn scalar k = `k'
    ereturn scalar s = `s'
    ereturn scalar lambda = `default_lambda'
    ereturn scalar kx = `kx'
    ereturn scalar p = `kx'
    ereturn matrix targets = `targets_mat', copy
    ereturn matrix x0_matrix = `targets_mat', copy
    ereturn matrix neff_targets = `neff_mat', copy
    ereturn matrix ess_targets = `neff_mat', copy
    ereturn matrix diagnostics = `target_diag', copy
    ereturn matrix target_diagnostics = `target_diag', copy
    if "`tuning_mode'" == "grid" {
        ereturn matrix tuning_candidates = `candidate_mat', copy
        ereturn matrix tuning_estimates = `grid_est', copy
        ereturn matrix tuning_se = `grid_se', copy
        ereturn matrix tuning_diagnostics = `grid_diag', copy
        ereturn scalar tuning_default_candidate = `default_index'
        ereturn scalar tuning_n_candidates = `ncandidates'
        ereturn scalar tuning_n_feasible = `nfeasible'
        ereturn scalar tuning_n_skipped = `ncandidates' - `nfeasible'
        ereturn scalar tuning_n_failed = 0
        ereturn scalar tuning_n_unique = `nunique_feasible'
    }
    if "`scalerule'"=="" local scalerule "sd"
    ereturn local scalerule "`scalerule'"
    ereturn scalar scaleeps=`scaleeps'
    ereturn local x0 "`x0'"
    ereturn local x0matrix "`x0matrix'"
    ereturn local target_source "`target_source'"
    ereturn local xvars "`xvars'"
    ereturn local x_cols "`xvars'"
    ereturn local distance "`distance'"
    ereturn local method "`method'"
    ereturn local tuning_mode "`tuning_mode'"
    ereturn local tuning_source "`tuning_source'"
    ereturn local tuning_status_codes "0=infeasible 1=ok"
    ereturn local tuning_reason_codes "0=ok 1=k_gt_s 2=jackknife_requires_s_lt_n 3=s_gt_n"
    ereturn local se_method = cond("`se'" == "", "rank_shift", "none")
    ereturn local jackknife = cond("`se'" == "", "fast", "none")
    ereturn local vce = cond("`se'" == "", "rank_shift", "none")
    ereturn local vcetype = cond("`se'" == "", "Jackknife", "")
    ereturn local depvar "`outcome'"
    ereturn local outcome "`outcome'"
    ereturn local signal "`outcome'"
    ereturn local cmd "np_dnn"
    ereturn local cmdline "np_dnn `0'"

    local signature_vars `xvars' `outcome'
    local signature_vars : list uniq signature_vars
    quietly _rdistnn_estimation_state `signature_vars'
    tempname model_groups model_signals model_ranks model_distances model_geometry
    mata: st_matrix("`model_groups'",_np_groups); st_matrix("`model_signals'",_np_y); st_matrix("`model_ranks'",_np_ranks); st_matrix("`model_distances'",_np_DX)
    local geometry_option ""
    if inlist("`e(distance)'","mahalanobis","l1_scaled","l2_scaled") {
        mata: st_matrix("`model_geometry'",_np_geometry)
        local geometry_option "geometry(`model_geometry')"
    }
    quietly _rdistnn_model_store `signature_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances') groups(`model_groups') `geometry_option'
    _np_dnn_Display, level(`level')
end

program define _np_dnn_Display, rclass
    version 17
    _rdistnn_report `0'
    return add
end
