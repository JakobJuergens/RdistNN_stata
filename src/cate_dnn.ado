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
        x(varlist numeric)               ///  legacy covariate option
        XCOLS(varlist numeric)           ///  R-aligned x_cols convention
        x0(numlist)                      ///  single evaluation point
        X0MATRIX(name)                   ///  named m x kx target matrix (R X0)
        TReat(varname numeric)           ///  binary treatment D
        TCOL(varname numeric)            ///  R-aligned t_col convention
        Propensity(varname numeric)      ///  estimated/oracle e(X)
        ECOL(varname numeric)            ///  R-aligned e_col convention
        M0(varname numeric)              ///  estimated/oracle E[Y|X,D=0]
        M0COL(varname numeric)           ///  R-aligned m0_col convention
        M1(varname numeric)              ///  estimated/oracle E[Y|X,D=1]
        M1COL(varname numeric)           ///  R-aligned m1_col convention
        MHAT(varname numeric)            ///  estimated/oracle E[Y|X]
        MCOL(varname numeric)            ///  R-aligned m_col convention
        PSEUDOOutcome                    ///  positional variable is CATE signal
        CLIP(real 1e-8)                  ///  propensity clipping constant
        k(integer 1)                     ///  # selected within each s-subset
        s(integer 0)                     ///  subset/kernel order; larger = tighter
        TUNINGMode(string)               ///  default | grid (R tuning_mode)
        LAMBDA(real .4)                  ///  paper finite-sample multiplier
        LAMBDAGrid(numlist)              ///  paper-rule candidate multipliers
        SGRID(numlist)                   ///  direct Stata candidate orders
        SDEFAULT(integer 0)              ///  default order for sgrid()
        Distance(string)                 ///  l2 (default) | mahalanobis
        METHod(string)                   ///  exact (default) | asymp
        Level(cilevel)                   ///  CI level (default 95)
        noSE                             ///  skip jackknife standard error
        noCONStant                       ///  (reserved; absorbed silently)
    ]

    local input_var `varlist'
    local outcome `input_var'
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

    if "`treat'" != "" & "`tcol'" != "" {
        di as error "Specify either treat() or tcol(), not both."
        exit 198
    }
    local treat_var `treat'
    if "`tcol'" != "" local treat_var `tcol'

    if "`propensity'" != "" & "`ecol'" != "" {
        di as error "Specify either propensity() or ecol(), not both."
        exit 198
    }
    local propensity_var `propensity'
    if "`ecol'" != "" local propensity_var `ecol'

    if "`m0'" != "" & "`m0col'" != "" {
        di as error "Specify either m0() or m0col(), not both."
        exit 198
    }
    local m0_var `m0'
    if "`m0col'" != "" local m0_var `m0col'

    if "`m1'" != "" & "`m1col'" != "" {
        di as error "Specify either m1() or m1col(), not both."
        exit 198
    }
    local m1_var `m1'
    if "`m1col'" != "" local m1_var `m1col'

    if "`mhat'" != "" & "`mcol'" != "" {
        di as error "Specify either mhat() or mcol(), not both."
        exit 198
    }
    local m_var `mhat'
    if "`mcol'" != "" local m_var `mcol'

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

    // ── resolve score interface ───────────────────────────────────────
    // The Stata interface mirrors the three R cate_dnn() score paths:
    // precomputed pseudo-outcome, full AIPW, and Robinson/R-learner.
    local score_method ""
    if "`pseudo_flag'" != "" {
        if "`treat_var'`propensity_var'`m0_var'`m1_var'`m_var'" != "" {
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
        if "`m_var'" != "" & ("`m0_var'" != "" | "`m1_var'" != "") {
            di as error "Specify either mcol() (Robinson) or m0col() and m1col() (AIPW), not both."
            exit 198
        }
        if "`m_var'" != "" {
            local score_method "robinson"
        }
        else if "`m0_var'" != "" & "`m1_var'" != "" {
            local score_method "aipw"
        }
        else {
            di as error "Supply mcol() for Robinson, or both m0col() and m1col() for AIPW."
            exit 198
        }
    }

    if missing(`clip') | `clip' <= 0 | `clip' >= .5 {
        di as error "clip() must be a finite number strictly between 0 and 0.5."
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
    else if "`score_method'" == "robinson" {
        markout `touse' `outcome' `treat_var' `xvars' `propensity_var' `m_var'
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
        di as error "k() must be >= 1."
        exit 198
    }

    if "`distance'" == "" local distance "l2"
    if !inlist("`distance'", "l2", "mahalanobis") {
        di as error "distance() must be l2 or mahalanobis."
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
        if "`score_method'" == "aipw" {
            // phi = m1-m0 + D*(Y-m1)/e - (1-D)*(Y-m0)/(1-e)
            quietly gen double `rho_var' = ///
                `treat_var' * (`outcome' - `m1_var') / `eclip' ///
                - (1 - `treat_var') * (`outcome' - `m0_var') / (1 - `eclip') ///
                + (`m1_var' - `m0_var') ///
                if `touse'
        }
        else {
            // Robinson/R-learner: (D-e)*(Y-m)/(e*(1-e))
            quietly gen double `rho_var' = ///
                (`treat_var' - `eclip') * (`outcome' - `m_var') ///
                / (`eclip' * (1 - `eclip')) if `touse'
        }
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
        _ct_DX    = rdistnn_distance_matrix(_ct_X, _ct_X0, "`distance'")
        _ct_ranks = rdistnn_ranks_from_dmat(_ct_DX)
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
                _ct_Wg = rdistnn_weights_from_rmat(_ct_ranks, `k', _ct_sg, "`method'")
                _ct_grid_est[_ct_g,.] = cross(_ct_rho, _ct_Wg)
                _ct_wr = rdistnn_weights_from_ranks((1::`n'), `k', _ct_sg, "`method'")
                _ct_topn = min((5, `n'))
                _ct_ess = 1 / quadsum(_ct_wr :^ 2)
                _ct_grid_diag[_ct_g,.] = (_ct_ess,
                    max(abs(_ct_wr)), _ct_topn,
                    quadsum(abs(_ct_wr[1.._ct_topn])) / quadsum(abs(_ct_wr)),
                    _ct_ess >= .9 * `n')
                if ("`se'" == "") {
                    _ct_Vg = rdistnn_jk_cov_multi_from_ranks(
                        _ct_ranks, _ct_rho, `k', _ct_sg, "`method'")
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
        _ct_target_diag = J(`ntargets', 25, .)
        if ("`score_method'" != "pseudo_outcome") {
            _ct_T = st_data(., "`treat_var'", "`touse'")
            _ct_E = st_data(., "`eclip'", "`touse'")
        }
        for (_ct_j = 1; _ct_j <= `ntargets'; _ct_j++) {
            _ct_wj = _ct_W[,_ct_j]
            _ct_dxj = _ct_DX[,_ct_j]
            _ct_word = order((_ct_wj, (1::`n')), (-1, 2))
            _ct_dxsort = sort(_ct_dxj, 1)
            _ct_target_diag[_ct_j, 1] = 1 / quadsum(_ct_wj :^ 2)
            _ct_target_diag[_ct_j, 2] = max(_ct_wj)
            _ct_target_diag[_ct_j, 3] = quadsum(_ct_wj[_ct_word[1]])
            _ct_target_diag[_ct_j, 4] = quadsum(
                _ct_wj[_ct_word[1..min((5, `n'))]])
            _ct_target_diag[_ct_j, 5] = quadsum(
                _ct_wj[_ct_word[1..min((10, `n'))]])
            _ct_target_diag[_ct_j, 6] = sum(
                abs(_ct_wj) :> sqrt(epsilon(1)))
            _ct_target_diag[_ct_j, 7] = 1
            _ct_target_diag[_ct_j, 8] = _ct_dxsort[min((`k', `n'))]
            _ct_target_diag[_ct_j, 9] = _ct_dxsort[min((`s', `n'))]
            _ct_target_diag[_ct_j, 10] = cross(_ct_dxj, _ct_wj) /
                quadsum(_ct_wj)
            _ct_target_diag[_ct_j, 11] = rdistnn_weighted_quantile(
                _ct_dxj, _ct_wj, .5)
            // Paper diagnostic: ordinary Euclidean radius in original X
            // coordinates, regardless of the metric used to form ranks.
            _ct_physical_dx = sqrt(rowsum((_ct_X :- _ct_X0[_ct_j,.]) :^ 2))
            _ct_target_diag[_ct_j, 20] = rdistnn_weighted_quantile(
                _ct_physical_dx, _ct_wj, .9)
            if ("`score_method'" == "pseudo_outcome") {
                _ct_target_diag[_ct_j, 12] = 0
                _ct_target_diag[_ct_j, 16] = 0
                _ct_target_diag[_ct_j, 21] = 0
            }
            else {
                _ct_tmass = cross(_ct_T, _ct_wj)
                _ct_cmass = cross(1 :- _ct_T, _ct_wj)
                _ct_target_diag[_ct_j, 12] = 1
                _ct_target_diag[_ct_j, 13] = _ct_tmass
                _ct_target_diag[_ct_j, 14] = _ct_cmass
                _ct_target_diag[_ct_j, 15] = _ct_tmass /
                    (_ct_tmass + _ct_cmass)
                _ct_target_diag[_ct_j, 16] = 1
                _ct_target_diag[_ct_j, 17] = cross(_ct_E, _ct_wj) /
                    quadsum(_ct_wj)
                _ct_target_diag[_ct_j, 18] = rdistnn_weighted_quantile(
                    _ct_E, _ct_wj, .1)
                _ct_target_diag[_ct_j, 19] = rdistnn_weighted_quantile(
                    _ct_E, _ct_wj, .9)
                _ct_target_diag[_ct_j, 21] = 1
                _ct_tw = select(_ct_wj, _ct_T :== 1)
                _ct_cw = select(_ct_wj, _ct_T :== 0)
                if (rows(_ct_tw) > 0 & quadsum(_ct_tw :^ 2) > 0) {
                    _ct_target_diag[_ct_j, 22] =
                        quadsum(_ct_tw)^2 / quadsum(_ct_tw :^ 2)
                }
                if (rows(_ct_cw) > 0 & quadsum(_ct_cw :^ 2) > 0) {
                    _ct_target_diag[_ct_j, 23] =
                        quadsum(_ct_cw)^2 / quadsum(_ct_cw :^ 2)
                }
            }
            _ct_target_diag[_ct_j, 24] = 1
            _ct_abs_contrib = abs(_ct_wj :* (_ct_rho :- _ct_tau[_ct_j]))
            _ct_contrib_total = quadsum(_ct_abs_contrib)
            if (_ct_contrib_total > 0 & _ct_contrib_total < .) {
                _ct_contrib_order = order((_ct_abs_contrib, (1::`n')), (-1, 2))
                _ct_target_diag[_ct_j, 25] = quadsum(
                    _ct_abs_contrib[_ct_contrib_order[1..min((5, `n'))]]) /
                    _ct_contrib_total
            }
        }
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
    ereturn local  x0       "`x0'"
    ereturn local  x0matrix "`x0matrix'"
    ereturn local  target_source "`target_source'"
    ereturn local  xvars    "`xvars'"
    ereturn local  x_cols   "`xvars'"
    ereturn local  x_source "`x_source'"
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
    ereturn local  outcome  "`outcome'"
    ereturn local  signal   "`input_var'"
    ereturn local  treat    "`treat_var'"
    ereturn local  t_col    "`treat_var'"
    ereturn local  propensity "`propensity_var'"
    ereturn local  e_col    "`propensity_var'"
    ereturn local  m0_col   "`m0_var'"
    ereturn local  m1_col   "`m1_var'"
    ereturn local  m_col    "`m_var'"
    ereturn local  cmd      "cate_dnn"
    ereturn local  cmdline  "cate_dnn `0'"

    local signature_vars `xvars' `input_var' `treat_var' `propensity_var' `m0_var' `m1_var' `m_var'
    local signature_vars : list uniq signature_vars
    quietly _rdistnn_estimation_state `signature_vars'
    tempname model_signals model_ranks model_distances
    mata: st_matrix("`model_signals'",_ct_rho); st_matrix("`model_ranks'",_ct_ranks); st_matrix("`model_distances'",_ct_DX)
    quietly _rdistnn_model_store `signature_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances')
    _cate_dnn_Display, level(`level')
end

// ── display subroutine ────────────────────────────────────────────────
program define _cate_dnn_Display, rclass
    version 17
    _rdistnn_report `0'
    return add
end
