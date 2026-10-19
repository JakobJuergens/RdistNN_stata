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
        x(varlist numeric)               ///  legacy covariate option
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
    if "`x'" != "" & "`xcols'" != "" {
        di as error "Specify either x() or xcols(), not both."
        exit 198
    }
    local xvars `x'
    local x_source "x"
    if "`xcols'" != "" {
        local xvars `xcols'
        local x_source "xcols"
    }
    if "`xvars'" == "" {
        di as error "xcols() is required (legacy alias: x())."
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

    marksample touse
    markout `touse' `xvars' `outcome'
    quietly count if `touse'
    if r(N) == 0 {
        di as error "No observations."
        exit 2000
    }
    local n = r(N)

    if `k' < 1 {
        di as error "k() must be >= 1."
        exit 198
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
    if `lambda' <= 0 | missing(`lambda') {
        di as error "lambda() must be a positive finite number."
        exit 198
    }
    if "`lambdagrid'" != "" & "`sgrid'" != "" {
        di as error "Specify either lambdagrid() or sgrid(), not both."
        exit 198
    }
    if "`tuning_mode'" == "default" & ///
       ("`lambdagrid'" != "" | "`sgrid'" != "" | `sdefault' != 0) {
        di as error "lambdagrid(), sgrid(), and sdefault() require tuningmode(grid)."
        exit 198
    }
    if "`tuning_mode'" == "grid" & `s' != 0 {
        di as error "s() cannot be combined with tuningmode(grid)."
        exit 198
    }
    if `sdefault' < 0 {
        di as error "sdefault() must be a positive integer when specified."
        exit 198
    }
    if `sdefault' > 0 & "`sgrid'" == "" {
        di as error "sdefault() requires sgrid()."
        exit 198
    }

    local paper_base = `n'^((`kx' + 1) / (`kx' + 2))
    local tuning_source "paper"
    local candidate_s ""
    local candidate_lambda ""
    local default_index 1

    if "`tuning_mode'" == "default" {
        if `s' == 0 local s = floor(`lambda' * `paper_base')
        else local tuning_source "s"
        local default_s = `s'
        local default_lambda = `lambda'
        local candidate_s "`s'"
        local candidate_lambda "`lambda'"
    }
    else if "`sgrid'" != "" {
        local tuning_source "sgrid"
        foreach sg of local sgrid {
            if `sg' < 1 | `sg' != floor(`sg') | missing(`sg') {
                di as error "sgrid() must contain positive integers."
                exit 198
            }
            local candidate_s "`candidate_s' `sg'"
            local candidate_lambda "`candidate_lambda' ."
        }
        local candidate_s : list retokenize candidate_s
        local candidate_lambda : list retokenize candidate_lambda
        local first_s : word 1 of `candidate_s'
        if `sdefault' == 0 local default_s = `first_s'
        else local default_s = `sdefault'
        local default_lambda = .
        local default_index 0
        local g = 0
        foreach sg of local candidate_s {
            local ++g
            if `default_index' == 0 & `sg' == `default_s' local default_index = `g'
        }
        if `default_index' == 0 {
            local candidate_s "`default_s' `candidate_s'"
            local candidate_lambda ". `candidate_lambda'"
            local default_index 1
        }
    }
    else {
        local tuning_source "lambda_grid"
        if "`lambdagrid'" == "" local lambdagrid ".3 .4 .5 .75"
        foreach lg of local lambdagrid {
            if `lg' <= 0 | missing(`lg') {
                di as error "lambdagrid() must contain positive finite numbers."
                exit 198
            }
            local sg = floor(`lg' * `paper_base')
            local candidate_s "`candidate_s' `sg'"
            local candidate_lambda "`candidate_lambda' `lg'"
        }
        local candidate_s : list retokenize candidate_s
        local candidate_lambda : list retokenize candidate_lambda
        local default_s = floor(`lambda' * `paper_base')
        local default_lambda = `lambda'
        local default_index 0
        local g = 0
        foreach sg of local candidate_s {
            local ++g
            local lg : word `g' of `candidate_lambda'
            if `default_index' == 0 & `sg' == `default_s' & ///
               abs(`lg' - `default_lambda') <= sqrt(c(epsdouble)) * max(1, abs(`default_lambda')) {
                local default_index = `g'
            }
        }
        if `default_index' == 0 {
            local g = 0
            foreach sg of local candidate_s {
                local ++g
                if `default_index' == 0 & `sg' == `default_s' local default_index = `g'
            }
        }
        if `default_index' == 0 {
            local candidate_s "`default_s' `candidate_s'"
            local candidate_lambda ". `candidate_lambda'"
            local default_index 1
        }
    }

    local s = `default_s'
    local ncandidates : word count `candidate_s'
    tempname candidate_mat grid_est grid_se grid_diag target_diag
    matrix `candidate_mat' = J(`ncandidates', 8, .)
    local candidate_names ""
    local seen_s ""
    local nfeasible 0
    local nunique 0
    local nunique_feasible 0
    forvalues g = 1/`ncandidates' {
        local sg : word `g' of `candidate_s'
        local lg : word `g' of `candidate_lambda'
        local feasible 1
        local reason 0
        if `k' > `sg' {
            local feasible 0
            local reason 1
        }
        else if "`se'" == "" & `sg' >= `n' {
            local feasible 0
            local reason 2
        }
        else if "`se'" != "" & `sg' > `n' {
            local feasible 0
            local reason 3
        }
        local computed_id 0
        local h = 0
        foreach old_s of local seen_s {
            local ++h
            if `computed_id' == 0 & `old_s' == `sg' local computed_id = `h'
        }
        if `computed_id' == 0 {
            local ++nunique
            local computed_id = `nunique'
            local seen_s "`seen_s' `sg'"
            if `feasible' local ++nunique_feasible
        }
        local nfeasible = `nfeasible' + `feasible'
        matrix `candidate_mat'[`g', 1] = `g'
        matrix `candidate_mat'[`g', 2] = `computed_id'
        matrix `candidate_mat'[`g', 3] = (`g' == `default_index')
        matrix `candidate_mat'[`g', 4] = `feasible'
        matrix `candidate_mat'[`g', 5] = `reason'
        matrix `candidate_mat'[`g', 6] = `feasible'
        matrix `candidate_mat'[`g', 7] = `lg'
        matrix `candidate_mat'[`g', 8] = `sg'
        local candidate_names "`candidate_names' candidate`g'"
    }
    local candidate_names : list retokenize candidate_names
    matrix rownames `candidate_mat' = `candidate_names'
    matrix colnames `candidate_mat' = candidate_id computed_id is_default feasible reason status lambda s

    if el(`candidate_mat', `default_index', 4) == 0 {
        local default_reason = el(`candidate_mat', `default_index', 5)
        if `default_reason' == 1 {
            di as error "Default s(`s') < k(`k'). Increase n or lambda(), or decrease k()."
        }
        else if `default_reason' == 2 {
            di as error "Default s(`s') must be less than n(`n') for jackknife SE; specify nose or reduce lambda()."
        }
        else {
            di as error "Default s(`s') must be no greater than n(`n') when nose is specified."
        }
        exit 198
    }

    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    tempname b V
    mata {
        _np_X = st_data(., "`xvars'", "`touse'")
        _np_X0 = st_matrix("`targets_mat'")
        _np_geometry = rdistnn_fit_geometry(_np_X,"`distance'","`scalerule'",strtoreal(tokens("`scalevec'")),`scaleeps')
        _np_DX = rdistnn_distance_matrix(_np_X, _np_X0, "`distance'", _np_geometry)
        _np_ranks = rdistnn_ranks_from_dmat(_np_DX)
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
                _np_Wg = rdistnn_weights_from_rmat(_np_ranks, `k', _np_sg, "`method'")
                _np_grid_est[_np_g,.] = rdistnn_weighted_signals(_np_y, _np_Wg)
                _np_wr = rdistnn_weights_from_ranks((1::`n'), `k', _np_sg, "`method'")
                _np_topn = min((5, `n'))
                _np_ess = 1 / quadsum(_np_wr :^ 2)
                _np_grid_diag[_np_g,.] = (_np_ess,
                    max(abs(_np_wr)), _np_topn,
                    quadsum(abs(_np_wr[1.._np_topn])) / quadsum(abs(_np_wr)),
                    _np_ess >= .9 * `n')
                if ("`se'" == "") {
                    _np_Vg = rdistnn_jk_cov_multi_from_ranks(
                        _np_ranks, _np_y, `k', _np_sg, "`method'")
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
        _np_target_diag = J(`ntargets', 17, .)
        for (_np_j = 1; _np_j <= `ntargets'; _np_j++) {
            _np_wj = _np_W[,_np_j]
            _np_dxj = _np_DX[,_np_j]
            _np_word = order((_np_wj, (1::`n')), (-1, 2))
            _np_dxsort = sort(_np_dxj, 1)
            _np_wsum = quadsum(_np_wj)
            _np_ybar = rdistnn_weighted_signals(_np_y, _np_wj) / _np_wsum
            _np_target_diag[_np_j, 1] = 1 / quadsum(_np_wj :^ 2)
            _np_target_diag[_np_j, 2] = max(_np_wj)
            _np_target_diag[_np_j, 3] = quadsum(_np_wj[_np_word[1]])
            _np_target_diag[_np_j, 4] = quadsum(
                _np_wj[_np_word[1..min((5, `n'))]])
            _np_target_diag[_np_j, 5] = quadsum(
                _np_wj[_np_word[1..min((10, `n'))]])
            _np_target_diag[_np_j, 6] = sum(
                abs(_np_wj) :> sqrt(epsilon(1)))
            _np_target_diag[_np_j, 7] = 1
            _np_target_diag[_np_j, 8] = _np_dxsort[min((`k', `n'))]
            _np_target_diag[_np_j, 9] = _np_dxsort[min((`s', `n'))]
            _np_target_diag[_np_j, 10] = cross(_np_dxj, _np_wj) / _np_wsum
            _np_target_diag[_np_j, 11] = rdistnn_weighted_quantile(
                _np_dxj, _np_wj, .5)
            _np_target_diag[_np_j, 12] = 1
            _np_target_diag[_np_j, 13] = _np_ybar
            _np_target_diag[_np_j, 14] = rdistnn_weighted_sd(_np_y, _np_wj)
            _np_target_diag[_np_j, 15] = rdistnn_weighted_quantile(
                _np_y, _np_wj, .1)
            _np_target_diag[_np_j, 16] = rdistnn_weighted_quantile(
                _np_y, _np_wj, .5)
            _np_target_diag[_np_j, 17] = rdistnn_weighted_quantile(
                _np_y, _np_wj, .9)
        }
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
    ereturn local x_source "`x_source'"
    ereturn local distance "`distance'"
    ereturn local method "`method'"
    ereturn local tuning_mode "`tuning_mode'"
    ereturn local tuning_source "`tuning_source'"
    ereturn local tuning_status_codes "0=infeasible 1=ok"
    ereturn local tuning_reason_codes "0=ok 1=k_gt_s 2=jackknife_requires_s_lt_n 3=s_gt_n"
    ereturn local se_method = cond("`se'" == "", "rank_shift", "none")
    ereturn local jackknife = cond("`se'" == "", "fast", "none")
    ereturn local outcome "`outcome'"
    ereturn local signal "`outcome'"
    ereturn local cmd "np_dnn"
    ereturn local cmdline "np_dnn `0'"

    local signature_vars `xvars' `outcome'
    local signature_vars : list uniq signature_vars
    quietly _rdistnn_estimation_state `signature_vars'
    tempname model_signals model_ranks model_distances model_geometry
    mata: st_matrix("`model_signals'",_np_y); st_matrix("`model_ranks'",_np_ranks); st_matrix("`model_distances'",_np_DX)
    local geometry_option ""
    if inlist("`e(distance)'","mahalanobis","l1_scaled","l2_scaled") {
        mata: st_matrix("`model_geometry'",_np_geometry)
        local geometry_option "geometry(`model_geometry')"
    }
    quietly _rdistnn_model_store `signature_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances') `geometry_option'
    _np_dnn_Display, level(`level')
end

program define _np_dnn_Display, rclass
    version 17
    _rdistnn_report `0'
    return add
end
