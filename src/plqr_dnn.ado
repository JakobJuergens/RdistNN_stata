*! version 0.1.0  2026-07-12  Jakob Juergens
*! Experimental PLQR localized-moment evaluation at supplied theta

program define plqr_dnn, eclass
    version 17

    if replay() {
        if "`e(cmd)'" != "plqr_dnn" error 301
        _plqr_dnn_Display `0'
        exit
    }

    syntax varname(numeric) [if] [in] , ///
        x(varlist numeric)               ///  covariates
        x0(numlist)                      ///  evaluation point
        D(varlist numeric)               ///  partially linear regressors
        MD(varlist numeric)              ///  supplied m_D(X), same order as d()
        [                                ///
        G(varname numeric)               ///  canonical supplied g_tau(X)
        Q(varname numeric)               ///  compatibility alias for g()
        THETAEVAL(numlist)                ///  canonical supplied theta vector
        beta(numlist)                    ///  compatibility alias for thetaeval()
        tau(real 0.5)                    ///  quantile index
        k(integer 1)                     ///  # selected within each s-subset
        s(integer 0)                     ///  subset/kernel order; larger = tighter
        Distance(string)                 ///  localization geometry
        SCALERule(string) SCALEVec(numlist) SCALEEPS(real 1e-12) ///
        Method(string)                   ///  exact (default) | asymp
        Level(cilevel)                   ///
        ]

    if "`g'" != "" & "`q'" != "" {
        di as error "Specify either g() or compatibility alias q(), not both."
        exit 198
    }
    if "`g'`q'" == "" {
        di as error "g() is required (compatibility alias: q())."
        exit 198
    }
    if "`thetaeval'" != "" & "`beta'" != "" {
        di as error "Specify either thetaeval() or compatibility alias beta(), not both."
        exit 198
    }
    if "`thetaeval'`beta'" == "" {
        di as error "thetaeval() is required (compatibility alias: beta())."
        exit 198
    }

    local gvar `g'
    local g_option "g"
    if "`gvar'" == "" {
        local gvar `q'
        local g_option "q"
    }
    local theta `thetaeval'
    local theta_option "thetaeval"
    if "`theta'" == "" {
        local theta `beta'
        local theta_option "beta"
    }

    marksample touse
    markout `touse' `varlist' `x' `d' `gvar' `md'

    local outcome `varlist'
    local xvars   `x'
    local dvars   `d'
    local mdvars  `md'
    local gvar    `gvar'
    local kx : word count `xvars'
    local kd : word count `dvars'
    local kmd : word count `mdvars'
    local ntheta : word count `theta'

    quietly count if `touse'
    if r(N) == 0 {
        di as error "No observations."
        exit 2000
    }
    local n = r(N)

    local nx0 : word count `x0'
    if `nx0' != `kx' {
        di as error "x0() must have `kx' elements (one per covariate in x())."
        exit 198
    }

    if `kd' != `kmd' {
        di as error "d() and md() must contain the same number of variables."
        exit 198
    }

    if `ntheta' != `kd' {
        di as error "thetaeval() must have `kd' elements (one per variable in d())."
        exit 198
    }

    if `tau' <= 0 | `tau' >= 1 {
        di as error "tau() must lie in (0,1)."
        exit 198
    }

    if `k' < 1 {
        di as error "k() must be >= 1."
        exit 198
    }

    if "`distance'" == "" local distance "l2"
    if !inlist("`distance'", "l2", "mahalanobis", "l1", "l1_scaled", "l2_scaled") {
        di as error "distance() must be l1, l2, l1_scaled, l2_scaled, or mahalanobis."
        exit 198
    }
    if "`distance'" == "mahalanobis" & `kx' < 2 {
        di as error "distance(mahalanobis) requires at least 2 covariates."
        exit 198
    }

    if "`method'" == "" local method "exact"
    if !inlist("`method'", "exact", "asymp") {
        di as error "method() must be exact or asymp."
        exit 198
    }

    if `s' == 0 {
        local s = floor(.4 * `n'^((`kx' + 1) / (`kx' + 2)))
    }
    if `s' < `k' {
        di as error "s(`s') < k(`k')."
        exit 198
    }
    if `s' >= `n' {
        di as error "s(`s') must be less than n(`n') for jackknife SE."
        exit 198
    }

    findfile "_rdistnn_mata.ado"
    quietly run "`r(fn)'"

    local x0_mata ""
    foreach v of local x0 {
        local x0_mata `x0_mata' `v',
    }
    local x0_mata = substr("`x0_mata'", 1, length("`x0_mata'") - 1)

    local theta_mata ""
    foreach v of local theta {
        local theta_mata `theta_mata' `v',
    }
    local theta_mata = substr("`theta_mata'", 1, length("`theta_mata'") - 1)

    tempname b V theta_mat moment_mat moment_vcov_mat
    mata {
        _pl_X      = st_data(., "`xvars'", "`touse'")
        _pl_y      = st_data(., "`outcome'", "`touse'")
        _pl_D      = st_data(., "`dvars'", "`touse'")
        _pl_g0     = st_data(., "`gvar'", "`touse'")
        _pl_md     = st_data(., "`mdvars'", "`touse'")
        _pl_x0     = (`x0_mata')
        _pl_theta  = (`theta_mata')'
        // Match the R reference and the other estimators: rank the reported
        // distance matrix, with original observation order resolving ties.
        _pl_geometry = rdistnn_fit_geometry(_pl_X,"`distance'","`scalerule'",strtoreal(tokens("`scalevec'")),`scaleeps')
        _pl_DX     = rdistnn_distance_matrix(_pl_X, _pl_x0, "`distance'", _pl_geometry)
        _pl_ranks  = rdistnn_ranks_from_dmat(_pl_DX)
        _pl_w      = rdistnn_weights_from_ranks(_pl_ranks, `k', `s', "`method'")
        _pl_resid  = `tau' :- (_pl_y :<= (_pl_D * _pl_theta :+ _pl_g0))
        _pl_score  = (_pl_D :- _pl_md) :* (_pl_resid * J(1, cols(_pl_D), 1))
        _pl_moment = rdistnn_weighted_signals(_pl_w, _pl_score)
        _pl_V      = rdistnn_jk_cov_from_ranks(_pl_ranks, _pl_score, `k', `s', "`method'")
        _pl_neff   = 1 / sum(_pl_w :^2)
        st_matrix("`b'", _pl_moment)
        st_matrix("`V'", _pl_V)
        st_matrix("`theta_mat'", _pl_theta')
        st_numscalar("r(neff)", _pl_neff)
    }

    matrix rownames `b' = moment
    matrix colnames `b' = `dvars'
    matrix rownames `V' = `dvars'
    matrix colnames `V' = `dvars'
    matrix colnames `theta_mat' = `dvars'
    matrix `moment_mat' = `b'
    matrix `moment_vcov_mat' = `V'

    ereturn post `b' `V', esample(`touse')
    ereturn matrix moment = `moment_mat', copy
    ereturn matrix moment_vcov = `moment_vcov_mat', copy
    ereturn matrix theta_eval = `theta_mat', copy
    ereturn matrix beta = `theta_mat', copy
    ereturn scalar N      = `n'
    ereturn scalar neff   = scalar(r(neff))
    ereturn scalar k      = `k'
    ereturn scalar s      = `s'
    ereturn scalar kx     = `kx'
    ereturn scalar kd     = `kd'
    ereturn scalar m      = 1
    ereturn scalar p      = `kx'
    ereturn scalar tau    = `tau'
    if "`scalerule'"=="" local scalerule "sd"
    ereturn local scalerule "`scalerule'"
    ereturn scalar scaleeps=`scaleeps'
    ereturn local  x0       "`x0'"
    ereturn local  xvars    "`xvars'"
    ereturn local  dvars    "`dvars'"
    ereturn local  mdvars   "`mdvars'"
    ereturn local  gvar     "`gvar'"
    ereturn local  qvar     "`gvar'"
    ereturn local  g_option "`g_option'"
    ereturn local  theta_option "`theta_option'"
    ereturn local  theta_source "supplied"
    ereturn local  inference_target "moment"
    ereturn local  root_solver "none"
    ereturn local  yvar     "`outcome'"
    ereturn local  distance "`distance'"
    ereturn local  method   "`method'"
    ereturn local  se_method "rank_shift"
    ereturn local  jackknife "fast"
    ereturn local  cmd      "plqr_dnn"
    ereturn local  cmdline  "plqr_dnn `0'"

    tempname model_signals model_ranks model_distances model_geometry
    mata: st_matrix("`model_signals'",_pl_score); st_matrix("`model_ranks'",_pl_ranks); st_matrix("`model_distances'",_pl_DX)
    local snapshot_vars `xvars' `outcome' `dvars' `gvar' `mdvars'
    local snapshot_vars : list uniq snapshot_vars
    local geometry_option ""
    if inlist("`e(distance)'","mahalanobis","l1_scaled","l2_scaled") {
        mata: st_matrix("`model_geometry'",_pl_geometry)
        local geometry_option "geometry(`model_geometry')"
    }
    quietly _rdistnn_model_store `snapshot_vars', signals(`model_signals') ranks(`model_ranks') distances(`model_distances') `geometry_option'
    tempname ses targets
    mata: st_matrix("`ses'",sqrt(diagonal(st_matrix("e(V)")))'); st_matrix("`targets'",_pl_x0)
    matrix colnames `ses'=`dvars'
    matrix colnames `targets'=`xvars'
    ereturn matrix estimates=`moment_mat', copy
    ereturn matrix se_targets=`ses', copy
    ereturn matrix estimate_vcov=`moment_vcov_mat', copy
    ereturn matrix targets=`targets', copy
    tempname bytarget diagnostics
    matrix `bytarget'=`moment_mat''
    ereturn matrix moment_by_target=`bytarget', copy
    foreach file in _rdistnn_clate_mata.ado _rdistnn_model_mata.ado {
        findfile "`file'"
        quietly run "`r(fn)'"
    }
    mata: st_matrix("`diagnostics'",rdistnn_model_diagnostics("plqr_dnn",_pl_X,_pl_x0,_pl_score,_pl_w,st_matrix("`model_distances'"),`k',`s',_pl_moment,J(0,1,.),J(0,1,.)))
    matrix colnames `diagnostics'=ess max_weight top1_weight_share top5_weight_share top10_weight_share n_nonzero_weight distance_available radius_rank_k radius_rank_s weighted_mean_distance weighted_median_distance
    ereturn matrix diagnostics=`diagnostics', copy
    ereturn matrix target_diagnostics=`diagnostics', copy
    quietly signestimationsample `snapshot_vars'
    quietly _rdistnn_sample_signature `snapshot_vars' if e(sample)
    ereturn local ordered_signature "`r(signature)'"
    ereturn local signature_vars "`snapshot_vars'"
    _plqr_dnn_Display, level(`level')
end

program define _plqr_dnn_Display, rclass
    syntax [, Level(cilevel)]
    _rdistnn_report, level(`level')
    return add
end
