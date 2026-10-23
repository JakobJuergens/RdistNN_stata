*! version 0.0.9  2026-07-12  Jakob Juergens
*! Conditional average treatment effect via distributional nearest neighbors

program define cate_dnn, eclass
    version 17

    // ── replay ────────────────────────────────────────────────────────
    if replay() {
        if "`e(cmd)'" != "cate_dnn" error 301
        _cate_dnn_Display `0'
        exit
    }

    syntax varname(numeric) [if] [in] , ///
        [                                ///
        XCOLS(varlist numeric)           ///  R-aligned x_cols convention
        x0(numlist)                      ///  single evaluation point
        X0MATRIX(name)                   ///  named m x kx target matrix (R X0)
        TCOL(varname numeric)            ///  R-aligned t_col convention
        ECOL(varname numeric)            ///  R-aligned e_col convention
        M0COL(varname numeric)           ///  R-aligned m0_col convention
        M1COL(varname numeric)           ///  R-aligned m1_col convention
        PSEUDOOutcome                    ///  positional variable is CATE signal
        CLIP(real 1e-8)                  ///  propensity clipping constant
        k(integer 1)                     ///  # selected within each s-subset
        s(integer 0)                     ///  subset/kernel order; larger = tighter
        TUNINGMode(string)               ///  default | grid (R tuning_mode)
        LAMBDA(real .4)                  ///  paper finite-sample multiplier
        LAMBDAGrid(numlist)              ///  paper-rule candidate multipliers
        SGRID(numlist)                   ///  direct Stata candidate orders
        SDEFAULT(integer 0)              ///  default order for sgrid()
        Distance(string)                 ///  localization geometry
        SCALERule(string) SCALEVec(numlist) SCALEEPS(real 1e-12) ///
        METHod(string)                   ///  exact (default) | asymp
        Level(cilevel)                   ///  CI level (default 95)
        noSE                             ///  skip jackknife standard error
    ]

    local input_var `varlist'
    local outcome `input_var'
    local xvars `xcols'
    if "`xvars'" == "" {
        di as error "xcols() is required."
        exit 198
    }
    local kx : word count `xvars'
    local treat_var `tcol'
    local propensity_var `ecol'
    local m0_var `m0col'
    local m1_var `m1col'

    local pseudo_flag "`pseudooutcome'"

    // ── resolve evaluation targets ────────────────────────────────────
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

    // ── resolve score interface ───────────────────────────────────────
    // CATE supports precomputed pseudo-outcomes and full AIPW construction.
    local score_method ""
    if "`pseudo_flag'" != "" {
        if "`treat_var'`propensity_var'`m0_var'`m1_var'" != "" {
            di as error "pseudo_outcome cannot be combined with nuisance options."
            exit 198
        }
        local score_method "pseudo_outcome"
        local outcome ""
    }
    else {
        if "`treat_var'" == "" | "`propensity_var'" == "" {
            di as error "tcol() and ecol() are required when constructing a pseudo-outcome."
            exit 198
        }
        if "`m0_var'" == "" | "`m1_var'" == "" {
            di as error "Supply both m0col() and m1col() for AIPW, or use pseudooutcome."
            exit 198
        }
        local score_method "aipw"
    }

    if missing(`clip') | `clip' <= 0 | `clip' >= .5 | 1-`clip'>=1 {
        di as error "clip() must be a finite number strictly between 0 and 0.5, with 1-clip() < 1."
        exit 198
    }

    // ── mark sample ──────────────────────────────────────────────────
    marksample touse, novarlist
    if "`score_method'" != "pseudo_outcome" {
        quietly count if `touse' & missing(`propensity_var')
        if r(N) > 0 {
            di as error "ecol() must contain finite values in [0,1]; missing values are not allowed."
            exit 198
        }
        quietly count if `touse' & (`propensity_var' < 0 | `propensity_var' > 1)
        if r(N) > 0 {
            di as error "ecol() must lie in [0,1] before clipping."
            exit 198
        }
    }
    if "`score_method'" == "pseudo_outcome" {
        markout `touse' `xvars' `input_var'
    }
    else {
        markout `touse' `outcome' `treat_var' `xvars' `propensity_var' `m0_var' `m1_var'
    }

    // ── validate ──────────────────────────────────────────────────────
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

    if "`method'" == "" local method "exact"
    if !inlist("`method'", "exact", "asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }

    // ── resolve scalar or grid tuning ────────────────────────────────
    local tuning_mode = lower("`tuningmode'")
    if "`tuning_mode'" == "" local tuning_mode "default"
    if !inlist("`tuning_mode'", "default", "grid") {
        di as error "tuningmode() must be default or grid."
        exit 198
    }
    if "`distance'" == "mahalanobis" & `n' < 2 {
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

    // ── pseudo-outcome ────────────────────────────────────────────────
    local rho_use `input_var'
    if "`score_method'" != "pseudo_outcome" {
        quietly count if `touse' & !inlist(`treat_var', 0, 1)
        if r(N) > 0 {
            di as error "tcol() must be coded 0/1."
            exit 198
        }

        tempvar eclip rho_var
        quietly gen double `eclip' = min(max(`propensity_var', `clip'), 1 - `clip') if `touse'
        // phi = m1-m0 + D*(Y-m1)/e - (1-D)*(Y-m0)/(1-e)
        quietly gen double `rho_var' = ///
            `treat_var' * (`outcome' - `m1_var') / `eclip' ///
            - (1 - `treat_var') * (`outcome' - `m0_var') / (1 - `eclip') ///
            + (`m1_var' - `m0_var') if `touse'
        local rho_use `rho_var'
    }

    // ── load Mata ────────────────────────────────────────────────────
    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    // ── compute ranks once, then evaluate unique candidate orders ─────
    tempname b V
    mata {
        _ct_X     = st_data(., "`xvars'", "`touse'")
        _ct_X0    = st_matrix("`targets_mat'")
        _ct_geometry = _ct_DX = _ct_ranks = J(0,0,.)
        rdistnn_localization_geometry(_ct_X,_ct_X0,"`distance'","`scalerule'",
            strtoreal(tokens("`scalevec'")),`scaleeps',_ct_geometry,_ct_DX,_ct_ranks)
        _ct_groups=rdistnn_distance_groups(_ct_ranks,_ct_DX,1)
        _ct_rho   = st_data(., "`rho_use'", "`touse'")
        _ct_candidates = st_matrix("`candidate_mat'")
        _ct_grid_est = J(`ncandidates', `ntargets', .)
        _ct_grid_se = J(`ncandidates', `ntargets', .)
        _ct_grid_diag = J(`ncandidates', 5, .)

        for (_ct_g = 1; _ct_g <= `ncandidates'; _ct_g++) {
            _ct_grid_diag[_ct_g, 3] = min((5, `n'))
            if (_ct_candidates[_ct_g, 4] == 0) continue
            _ct_prior = 0
            if (_ct_g > 1) {
                for (_ct_h = 1; _ct_h < _ct_g; _ct_h++) {
                    if (_ct_candidates[_ct_h, 4] == 1 &
                        _ct_candidates[_ct_h, 2] == _ct_candidates[_ct_g, 2]) {
                        _ct_prior = _ct_h
                        break
                    }
                }
            }
            if (_ct_prior > 0) {
                _ct_grid_est[_ct_g,.] = _ct_grid_est[_ct_prior,.]
                _ct_grid_se[_ct_g,.] = _ct_grid_se[_ct_prior,.]
                _ct_grid_diag[_ct_g,.] = _ct_grid_diag[_ct_prior,.]
            }
            else {
                _ct_sg = _ct_candidates[_ct_g, 8]
                _ct_Wg = rdistnn_weights_from_rmat(_ct_ranks, `k', _ct_sg, "`method'",_ct_DX,_ct_groups)
                _ct_grid_est[_ct_g,.] = rdistnn_weighted_signals(_ct_rho, _ct_Wg)
                _ct_grid_diag[_ct_g,.]=rdistnn_grid_weight_diagnostics(_ct_Wg)
                if ("`se'" == "") {
                    _ct_Vg = rdistnn_jk_cov_multi_from_ranks(
                        _ct_ranks, _ct_rho, `k', _ct_sg, "`method'",_ct_DX,_ct_groups)
                    _ct_grid_se[_ct_g,.] = sqrt(diagonal(_ct_Vg))'
                }
                // Save the default's first representative, including when the
                // default row duplicates an earlier candidate after other orders.
                if (_ct_candidates[_ct_g, 2] == _ct_candidates[`default_index', 2]) {
                    _ct_W = _ct_Wg
                    if ("`se'" == "") _ct_V = _ct_Vg
                }
            }
        }

        _ct_tau = _ct_grid_est[`default_index',.]
        _ct_neff = 1 :/ colsum(_ct_W :^ 2)
        _ct_T = _ct_E = J(0,1,.)
        if ("`score_method'" != "pseudo_outcome") {
            _ct_T = st_data(., "`treat_var'", "`touse'")
            _ct_E = st_data(., "`eclip'", "`touse'")
        }
        _ct_target_diag = rdistnn_diagnostics("cate_dnn",_ct_X,_ct_X0,
            _ct_rho,_ct_W,_ct_DX,`k',`s',_ct_tau,_ct_T,_ct_E)
        st_matrix("`b'", _ct_tau)
        st_matrix("`neff_mat'", _ct_neff)
        st_matrix("`grid_est'", _ct_grid_est)
        st_matrix("`grid_se'", _ct_grid_se)
        st_matrix("`grid_diag'", _ct_grid_diag)
        st_matrix("`target_diag'", _ct_target_diag)
        st_numscalar("r(neff_min)", min(_ct_neff))
        st_numscalar("r(neff_max)", max(_ct_neff))
    }

    if "`se'" == "" mata: st_matrix("`V'", _ct_V)

    // ── store eclass results ─────────────────────────────────────────
    local coef_names "tau"
    if `ntargets' > 1 {
        local coef_names ""
        forvalues j = 1/`ntargets' {
            local coef_names `coef_names' tau`j'
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
        weighted_mean_distance weighted_median_distance overlap_available ///
        treated_mass control_mass treated_share propensity_available ///
        weighted_mean_propensity weighted_propensity_q10 ///
        weighted_propensity_q90 physical_radius_90 arm_ess_available ///
        treated_ess control_ess score_contribution_available ///
        top5_score_share
    if "`se'" == "" {
        matrix rownames `V' = `coef_names'
        matrix colnames `V' = `coef_names'
        ereturn post `b' `V', esample(`touse')
    }
    else {
        ereturn post `b', esample(`touse')
    }
    ereturn scalar N      = `n'
    ereturn scalar m      = `ntargets'
    ereturn scalar neff_min = scalar(r(neff_min))
    ereturn scalar neff_max = scalar(r(neff_max))
    ereturn scalar ess_min = scalar(r(neff_min))
    ereturn scalar ess_max = scalar(r(neff_max))
    if `ntargets' == 1 {
        ereturn scalar neff = el(`neff_mat', 1, 1)
        ereturn scalar ess  = el(`neff_mat', 1, 1)
    }
    ereturn scalar k      = `k'
    ereturn scalar s      = `s'
    ereturn scalar lambda = `default_lambda'
    ereturn scalar kx     = `kx'
    ereturn scalar p      = `kx'
    ereturn scalar clip   = `clip'
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
    ereturn local  x0       "`x0'"
    ereturn local  x0matrix "`x0matrix'"
    ereturn local  target_source "`target_source'"
    ereturn local  xvars    "`xvars'"
    ereturn local  x_cols   "`xvars'"
    ereturn local  distance "`distance'"
    ereturn local  method   "`method'"
    ereturn local  tuning_mode "`tuning_mode'"
    ereturn local  tuning_source "`tuning_source'"
    ereturn local  tuning_status_codes "0=infeasible 1=ok"
    ereturn local  tuning_reason_codes "0=ok 1=k_gt_s 2=jackknife_requires_s_lt_n 3=s_gt_n"
    ereturn local  score_method "`score_method'"
    ereturn local  propensity_version = cond("`score_method'" == "pseudo_outcome", "not_used", "clipped")
    ereturn local  se_method = cond("`se'" == "", "rank_shift", "none")
    ereturn local  jackknife = cond("`se'" == "", "fast", "none")
    ereturn local  vce = cond("`se'" == "", "rank_shift", "none")
    ereturn local  vcetype = cond("`se'" == "", "Frozen-score jackknife", "")
    ereturn local  depvar "`input_var'"
    ereturn local  outcome  "`outcome'"
    ereturn local  signal   "`input_var'"
    ereturn local  t_col    "`treat_var'"
    ereturn local  e_col    "`propensity_var'"
    ereturn local  m0_col   "`m0_var'"
    ereturn local  m1_col   "`m1_var'"
    ereturn local  cmd      "cate_dnn"
    ereturn local  cmdline  "cate_dnn `0'"

    local signature_vars `xvars' `input_var' `treat_var' `propensity_var' `m0_var' `m1_var'
    local signature_vars : list uniq signature_vars
    quietly _rdistnn_estimation_state `signature_vars'
    tempname model_groups model_signals model_ranks model_distances model_geometry
    mata: st_matrix("`model_groups'",_ct_groups); st_matrix("`model_signals'",_ct_rho); st_matrix("`model_ranks'",_ct_ranks); st_matrix("`model_distances'",_ct_DX)
    local geometry_option ""
    if inlist("`e(distance)'","mahalanobis","l1_scaled","l2_scaled") {
        mata: st_matrix("`model_geometry'",_ct_geometry)
        local geometry_option "geometry(`model_geometry')"
    }
    quietly _rdistnn_model_store `signature_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances') groups(`model_groups') `geometry_option'
    _cate_dnn_Display, level(`level')
end

// ── display subroutine ────────────────────────────────────────────────
program define _cate_dnn_Display, rclass
    version 17
    _rdistnn_report `0'
    return add
end
